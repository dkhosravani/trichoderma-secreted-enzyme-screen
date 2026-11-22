"""Calculate theoretical full-length and post-signal-segment properties."""
import csv
import json
from collections import Counter
from pathlib import Path
from Bio.SeqUtils.ProtParam import ProteinAnalysis

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / 'results'
VALID_AA = set('ACDEFGHIKLMNPQRSTVWY')


def calculate(sequence):
    sequence = (sequence or '').strip().upper()
    if not sequence:
        return '', '', 'missing_sequence'
    unsupported = set(sequence) - VALID_AA
    if unsupported:
        return '', '', 'unsupported_residues:' + ''.join(sorted(unsupported))
    protein = ProteinAnalysis(sequence)
    return round(protein.molecular_weight() / 1000, 3), round(protein.isoelectric_point(), 3), 'ok'


def write_csv(path, rows):
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    with (ROOT / 'data/proteome_records.json').open(encoding='utf-8') as handle:
        records = json.load(handle)
    with (RESULTS / 'protein_annotations.csv').open(encoding='utf-8', newline='') as handle:
        rows = list(csv.DictReader(handle))

    sequences = {r['primaryAccession']: (r.get('sequence') or {}).get('value', '') for r in records}
    if len(sequences) != len(records) or len(rows) != len(sequences):
        raise ValueError('Duplicate or missing source records')

    counts = Counter()
    candidates = []
    for row in rows:
        accession = row['uniprot_accession']
        if accession not in sequences:
            raise ValueError(f'Unknown accession: {accession}')
        sequence = sequences[accession]
        if sequence and int(row['length_aa']) != len(sequence):
            raise ValueError(f'Length mismatch: {accession}')

        mw, pi, full_status = calculate(sequence)
        row['calculated_mw_kda'] = mw
        row['calculated_pi'] = pi
        row['property_status'] = full_status
        counts[full_status] += 1

        # A UniProt signal feature may be predicted. Removing the annotated
        # signal region DOES NOT establish the experimentally mature form.
        if row['signal_location_status'] == 'valid_annotated_coordinates':
            end = int(row['signal_end_aa'])
            clipped = sequence[end:]
            clipped_mw, clipped_pi, clipped_status = calculate(clipped)
            row['post_signal_length_aa'] = len(clipped)
        else:
            clipped_mw, clipped_pi = '', ''
            clipped_status = 'not_calculated:' + row['signal_location_status']
            row['post_signal_length_aa'] = ''
        row['post_signal_mw_kda'] = clipped_mw
        row['post_signal_pi'] = clipped_pi
        row['post_signal_status'] = clipped_status

        if row['signal_peptide_annotated'] == 'True' and row['glycosidase_ec_3_2_1'] == 'True':
            candidates.append(row)

    write_csv(RESULTS / 'protein_properties.csv', rows)
    write_csv(RESULTS / 'candidates.csv', candidates)

    summary = {
        'total_proteins': len(rows),
        'properties_calculated': counts['ok'],
        'properties_skipped': len(rows) - counts['ok'],
        'skip_reasons': {k: v for k, v in counts.items() if k != 'ok'},
        'candidate_count': len(candidates),
        'candidates_with_properties': sum(r['property_status'] == 'ok' for r in candidates),
        'post_signal_properties_calculated': sum(r['post_signal_status'] == 'ok' for r in rows),
        'candidate_post_signal_properties_calculated': sum(r['post_signal_status'] == 'ok' for r in candidates),
    }
    with (RESULTS / 'property_summary.json').open('w', encoding='utf-8') as handle:
        json.dump(summary, handle, indent=2)
    print('PROPERTY SUMMARY')
    for key, value in summary.items():
        print(f'{key}: {value}')


if __name__ == '__main__':
    main()
