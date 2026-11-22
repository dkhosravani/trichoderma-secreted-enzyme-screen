"""Unit tests use synthetic records and do not claim experimental validation."""
from src.parse_annotations import (extract_ec_numbers, parse_record,
                                   evidence_group, signal_cleavage_position)
from src.calculate_properties import calculate
from src.select_candidates import select_shortlist, valid_ec_numbers
from src.audit_pipeline import percentile


def record(ec=None, signal=False, sequence='ACDEFGHIKLMN'):
    name = {'fullName': {'value': 'Example enzyme'}}
    if ec is not None:
        name['ecNumbers'] = [{'value': ec}]
    return {
        'primaryAccession': 'TEST001',
        'entryType': 'UniProtKB reviewed (Swiss-Prot)',
        'proteinDescription': {'recommendedName': name},
        'genes': [{'geneName': {'value': 'example_gene'}}],
        'sequence': {'value': sequence},
        'features': ([{'type': 'Signal', 'location': {'start': {'value': 1},
                    'end': {'value': 3}},
                    'evidences': [{'evidenceCode': 'ECO:0000256',
                                   'source': 'SAM', 'id': 'SignalP'}]}]
                     if signal else []),
    }


def test_signal_true():
    assert parse_record(record(signal=True))['signal_peptide_annotated']


def test_signal_missing():
    assert not parse_record(record())['signal_peptide_annotated']


def test_transmembrane_not_signal():
    item = record()
    item['features'] = [{'type': 'Transmembrane'}]
    assert not parse_record(item)['signal_peptide_annotated']


def test_ec_true():
    assert parse_record(record(ec='3.2.1.4'))['glycosidase_ec_3_2_1']


def test_ec_false():
    assert not parse_record(record(ec='1.1.1.1'))['glycosidase_ec_3_2_1']


def test_ec_incomplete():
    parsed = parse_record(record(ec='3.2.1.-'))
    assert not parsed['glycosidase_ec_3_2_1']
    assert parsed['incomplete_ec_3_2_1']


def test_ec_missing():
    assert parse_record(record())['ec_numbers'] == ''


def test_ec_nested():
    description = {'alternativeNames': [{'fullName': {'value': 'Alt'},
                    'ecNumbers': [{'value': '3.2.1.8'}]}]}
    assert extract_ec_numbers(description) == ['3.2.1.8']


def test_length_from_sequence():
    assert parse_record(record(sequence='ACDEFG'))['length_aa'] == 6


def test_names():
    result = parse_record(record())
    assert result['protein_name'] == 'Example enzyme'
    assert result['gene_name'] == 'example_gene'


def test_evidence_codes():
    result = parse_record(record(signal=True))
    assert result['signal_evidence_codes'] == 'ECO:0000256'
    assert result['signal_evidence_group'] == 'ECO:0000256_only'


def test_manual_model_evidence():
    assert evidence_group(['ECO:0000255']) == 'ECO:0000255_only'


def test_unrecognized_evidence():
    assert evidence_group(['ECO:0000269']) == 'other_or_mixed_evidence'


def test_cleavage_coordinates():
    value = parse_record(record(signal=True))
    assert value['signal_end_aa'] == 3
    assert value['signal_location_status'] == 'valid_annotated_coordinates'


def test_bad_cleavage_coordinates():
    features = [{'type': 'Signal', 'location': {'start': {'value': 2},
                 'end': {'value': 4}}}]
    assert signal_cleavage_position(features, 12)[0] == ''


def test_no_cleavage_coordinates():
    assert signal_cleavage_position([], 12)[1] == 'no_signal_annotation'


def test_valid_properties():
    mw, pi, status = calculate('ACDEFGHIK')
    assert status == 'ok' and mw > 0 and 0 < pi < 14


def test_ambiguous_sequence():
    assert calculate('ACDEX')[2] == 'unsupported_residues:X'


def test_empty_sequence():
    assert calculate('') == ('', '', 'missing_sequence')


def test_valid_ec_selection():
    row = {'ec_numbers': '3.2.1.4; 3.2.1.-; 1.1.1.1'}
    assert valid_ec_numbers(row) == ['3.2.1.4']


def shortlist_row(accession, ec, reviewed='False'):
    return {'uniprot_accession': accession, 'gene_name': 'x',
            'protein_name': 'Example', 'reviewed': reviewed,
            'signal_peptide_annotated': 'True',
            'glycosidase_ec_3_2_1': 'True', 'property_status': 'ok',
            'ec_numbers': ec}


def test_shortlist_prefers_reviewed():
    rows = [shortlist_row('A', '3.2.1.4'),
            shortlist_row('B', '3.2.1.4', 'True')]
    selected, eligible, categories = select_shortlist(rows, target=1)
    assert (selected[0]['uniprot_accession'], eligible, categories) == ('B', 2, 1)


def test_shortlist_distinct_categories():
    rows = [shortlist_row('A', '3.2.1.4'), shortlist_row('B', '3.2.1.8')]
    selected, _, _ = select_shortlist(rows, target=2)
    assert len({r['representative_ec'] for r in selected}) == 2


def test_percentile():
    assert percentile([1, 2, 3, 4, 5], .25) == 2
