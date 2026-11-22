"""Extract UniProt feature/EC annotations without inferring true secretion."""
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / 'results'
EC_GLYCOSIDASE = re.compile(r'3\.2\.1\.\d+\Z')


def extract_ec_numbers(description):
    numbers = set()

    def visit(value):
        if isinstance(value, dict):
            for key, item in value.items():
                if key == 'ecNumbers' and isinstance(item, list):
                    for ec in item:
                        if isinstance(ec, dict) and ec.get('value'):
                            numbers.add(ec['value'])
                else:
                    visit(item)
        elif isinstance(value, list):
            for item in value:
                visit(item)

    visit(description)
    return sorted(numbers)


def get_protein_name(description):
    names = [description.get('recommendedName') or {}]
    names += description.get('submissionNames') or []
    names += description.get('alternativeNames') or []
    for name in names:
        value = (name.get('fullName') or {}).get('value')
        if value:
            return value
    return ''


def evidence_group(codes):
    codes = set(codes)
    if not codes:
        return 'no_evidence_code'
    if codes == {'ECO:0000256'}:
        return 'ECO:0000256_only'
    if codes == {'ECO:0000255'}:
        return 'ECO:0000255_only'
    if codes == {'ECO:0000255', 'ECO:0000256'}:
        return 'both_model_codes'
    return 'other_or_mixed_evidence'


def signal_cleavage_position(features, sequence_length):
    """Return end of ONE annotated N-terminal signal region, or a reason."""
    signals = [f for f in features if f.get('type') == 'Signal']
    if not signals:
        return '', 'no_signal_annotation'
    if len(signals) != 1:
        return '', 'multiple_signal_features'
    location = signals[0].get('location') or {}
    start = (location.get('start') or {}).get('value')
    end = (location.get('end') or {}).get('value')
    if not isinstance(start, int) or not isinstance(end, int):
        return '', 'missing_or_uncertain_coordinates'
    if start != 1 or not 1 <= end < sequence_length:
        return '', 'non_nt_or_out_of_range'
    return end, 'valid_annotated_coordinates'


def parse_record(record):
    description = record.get('proteinDescription') or {}
    features = record.get('features') or []
    sequence = (record.get('sequence') or {}).get('value') or ''
    signals = [f for f in features if f.get('type') == 'Signal']
    ec = extract_ec_numbers(description)

    evidences = [
        e for feature in signals for e in (feature.get('evidences') or [])
        if isinstance(e, dict)
    ]
    codes = sorted({e['evidenceCode'] for e in evidences if e.get('evidenceCode')})
    sources = sorted({
        f"{e['source']}:{e.get('id', '')}"
        for e in evidences if e.get('source')
    })
    cleavage_end, cleavage_status = signal_cleavage_position(features, len(sequence))

    gene = ''
    for gene_info in record.get('genes') or []:
        gene = (gene_info.get('geneName') or {}).get('value', '')
        if gene:
            break

    return {
        'uniprot_accession': record['primaryAccession'],
        'gene_name': gene,
        'protein_name': get_protein_name(description),
        'reviewed': record.get('entryType') == 'UniProtKB reviewed (Swiss-Prot)',
        'signal_peptide_annotated': bool(signals),
        'signal_evidence_codes': '; '.join(codes),
        'signal_evidence_sources': '; '.join(sources),
        'signal_evidence_group': evidence_group(codes) if signals else 'no_signal_annotation',
        'signal_end_aa': cleavage_end,
        'signal_location_status': cleavage_status,
        'transmembrane_annotated': any(f.get('type') == 'Transmembrane' for f in features),
        'ec_numbers': '; '.join(ec),
        'glycosidase_ec_3_2_1': any(EC_GLYCOSIDASE.fullmatch(number) for number in ec),
        'incomplete_ec_3_2_1': '3.2.1.-' in ec,
        'length_aa': len(sequence) if sequence else '',
    }


def write_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    with (ROOT / 'data/proteome_records.json').open(encoding='utf-8') as handle:
        records = json.load(handle)
    with (RESULTS / 'metadata.json').open(encoding='utf-8') as handle:
        metadata = json.load(handle)

    accessions = [r['primaryAccession'] for r in records]
    if len(records) != metadata['downloaded_count'] or len(set(accessions)) != len(records):
        raise ValueError('Downloaded count or accession uniqueness failed')

    rows = sorted((parse_record(r) for r in records), key=lambda x: x['uniprot_accession'])
    signals = [r for r in rows if r['signal_peptide_annotated']]
    gh = [r for r in rows if r['glycosidase_ec_3_2_1']]
    candidates = [r for r in signals if r['glycosidase_ec_3_2_1']]

    summary = {
        'total_proteins': len(rows),
        'annotated_signal_peptide': len(signals),
        'no_signal_peptide_annotation': len(rows) - len(signals),
        'glycosidase_ec_3_2_1': len(gh),
        'signal_and_glycosidase': len(candidates),
        'glycosidase_without_signal_annotation': len(gh) - len(candidates),
        'signal_annotated_without_ec': sum(not r['ec_numbers'] for r in signals),
        'incomplete_ec_3_2_1_entries': sum(r['incomplete_ec_3_2_1'] for r in rows),
        'missing_sequences': sum(r['length_aa'] == '' for r in rows),
    }
    RESULTS.mkdir(exist_ok=True)
    write_csv(RESULTS / 'protein_annotations.csv', rows)
    with (RESULTS / 'annotation_summary.json').open('w', encoding='utf-8') as handle:
        json.dump(summary, handle, indent=2)
    print('ANNOTATION SUMMARY')
    for key, value in summary.items():
        print(f'{key}: {value}')


if __name__ == '__main__':
    main()
