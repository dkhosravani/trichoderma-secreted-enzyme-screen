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

# Added after the scientific review: distribution tests and reference recovery.
from src.audit_pipeline import compare_distributions, group_stats, fmt_p
from src.check_reference_enzymes import reference_results


def test_mann_whitney_effect_sign():
    result = compare_distributions([1, 2, 3], [8, 9, 10])
    assert result['rank_biserial'] == -1
    assert 0 <= result['p_value'] <= 1
    assert result['n_signal'] == result['n_no_signal'] == 3


def test_mann_whitney_reverse_direction():
    assert compare_distributions([8, 9], [1, 2])['rank_biserial'] == 1


def test_mann_whitney_empty_group():
    result = compare_distributions([], [1, 2])
    assert result['p_value'] is None and result['rank_biserial'] is None


def test_iqr_with_missing_values():
    result = group_stats([{'p': '1'}, {'p': '2'}, {'p': ''},
                          {'p': '3'}, {'p': '4'}, {'p': '5'}], 'p')
    assert result == {'n_valid': 5, 'median': 3, 'q1': 2, 'q3': 4}


def test_uncertain_signal_coordinates_not_clipped():
    feature = {'type': 'Signal', 'location': {
        'start': {'value': 1}, 'end': {'value': 3, 'modifier': 'UNSURE'}}}
    assert signal_cleavage_position([feature], 10)[0] == ''


def test_reference_check_known_ec_incomplete():
    row = shortlist_row('G0RVK1', '3.2.1.-')
    row['glycosidase_ec_3_2_1'] = 'False'
    row['gene_name'] = 'cel7a'
    refs = [{'reference': 'CBH I', 'accession': 'G0RVK1', 'gene': 'cel7a',
             'role': 'positive-reference'}]
    results = reference_results([row], references=refs)
    assert results[0]['match_status'] == 'matched'
    assert results[0]['primary_screen_pass'] is False


def test_reference_check_gene_search_and_missing():
    row = shortlist_row('KNOWN', '3.2.1.4')
    row['gene_name'] = 'cel6a'
    refs = [
        {'reference': 'CBH II', 'accession': '', 'gene': 'cel6a',
         'role': 'positive-reference'},
        {'reference': 'Other', 'accession': 'MISSING', 'gene': '',
         'role': 'comparison'},
    ]
    checks = reference_results([row], refs)
    assert checks[0]['accession'] == 'KNOWN' and checks[0]['primary_screen_pass']
    assert checks[1]['match_status'] == 'not_found'


def test_deterministic_shortlist_order():
    rows = [shortlist_row('Z', '3.2.1.4'), shortlist_row('A', '3.2.1.4'),
            shortlist_row('C', '3.2.1.8')]
    first = select_shortlist(rows, target=2)[0]
    second = select_shortlist(list(reversed(rows)), target=2)[0]
    assert [(r['uniprot_accession'], r['representative_ec']) for r in first] == [
        (r['uniprot_accession'], r['representative_ec']) for r in second]


def test_p_value_formatter():
    assert 'underflow' in fmt_p(0)
    assert fmt_p(None) == 'not computed'


def test_readme_auto_results_no_double_marker(tmp_path, monkeypatch):
    from src import audit_pipeline
    (tmp_path / 'README.md').write_text('Header\n' + audit_pipeline.START + '\nold\n'
                                        + audit_pipeline.END + '\nFooter\n')
    monkeypatch.setattr(audit_pipeline, 'ROOT', tmp_path)
    audit_pipeline.update_readme_generated_section('new finding')
    text = (tmp_path / 'README.md').read_text()
    assert 'new finding' in text and 'old' not in text
    assert text.count(audit_pipeline.START) == text.count(audit_pipeline.END) == 1


def test_transmembrane_outside_signal_flag():
    entry = record(signal=True)
    entry['features'].append({'type': 'Transmembrane',
                             'location': {'start': {'value': 6}, 'end': {'value': 9}}})
    assert parse_record(entry)['tm_outside_signal_status'] == 'tm_annotation_outside_signal_region'


def test_transmembrane_inside_signal_flag():
    entry = record(signal=True)
    entry['features'].append({'type': 'Transmembrane',
                             'location': {'start': {'value': 1}, 'end': {'value': 2}}})
    assert parse_record(entry)['tm_outside_signal_status'] == 'tm_features_within_signal_region'


def test_qm6a_cbhii_reference_uses_accession_not_only_gene_name():
    from src.check_reference_enzymes import REFERENCES
    ref = next(r for r in REFERENCES if r['reference'].startswith('CBH II'))
    assert ref['accession'] == 'G0RB58'
    example = shortlist_row('G0RB58', '3.2.1.91')
    example['gene_name'] = ''
    checks = reference_results([example], references=[ref])
    assert checks[0]['match_status'] == 'matched'
    assert checks[0]['accession'] == 'G0RB58'
