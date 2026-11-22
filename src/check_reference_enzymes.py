"""Exploratory positive/negative controls with documented identifiers.

These are *not* the list used to define the screen. The check happens
AFTER classification to avoid selection bias. A protein absent from our
strict EC-defined set is not automatically a false secretion prediction.
"""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / 'results'

# Three accessions came from inspection of the QM6a UniProt records.
# CBH II (Cel6A) can be represented in QM6a as G0RB58, annotated as
# 'Glucanase' rather than under the cel6a gene name. The candidate mapping
# was identified through the UniProt UniRef100_A0A024SH76 cluster, containing
# G0RB58 (QM6a) together with Exoglucanase 2 records. This is a reference
# mapping for manual inspection, not independent experimental validation.
REFERENCES = [
    {'reference': 'CBH I (Cel7A)', 'accession': 'G0RVK1', 'gene': 'cel7a', 'role': 'positive-reference'},
    {'reference': 'CBH II (Cel6A, UniRef-linked)', 'accession': 'G0RB58', 'gene': 'cel6a', 'role': 'positive-reference'},
    {'reference': 'Endoglucanase II (Cel5A)', 'accession': 'G0RB67', 'gene': 'cel5a', 'role': 'positive-reference'},
    {'reference': 'Xylanase II (Xyn2)', 'accession': 'G0RUP7', 'gene': 'xyn2', 'role': 'positive-reference'},
    {'reference': 'eIF5A comparison', 'accession': 'G0RKY3', 'gene': '', 'role': 'comparison'},
]


def reference_results(rows, references=REFERENCES):
    by_accession = {row['uniprot_accession']: row for row in rows}
    by_gene = {}
    for row in rows:
        # A single canonical gene name is parsed by the upstream module.
        gene = row['gene_name'].strip().lower()
        if gene:
            by_gene.setdefault(gene, []).append(row)

    output = []
    for ref in references:
        accession = ref['accession']
        gene = ref['gene']
        match_method = 'accession' if accession else 'gene_name'
        matches = [by_accession[accession]] if accession and accession in by_accession else []
        if not matches and gene:
            matches = by_gene.get(gene.lower(), [])
            match_method = 'gene_name'
        if len(matches) != 1:
            output.append({
                'reference': ref['reference'], 'expected_role': ref['role'],
                'accession': accession, 'gene_name': gene,
                'match_method': match_method,
                'match_status': 'not_found' if not matches else 'ambiguous_gene_match',
                'signal_annotated': '', 'complete_ec_3_2_1': '',
                'primary_screen_pass': '', 'ec_numbers': '',
                'interpretation': 'Not resolved from these annotations; no recovery claim',
            })
            continue
        row = matches[0]
        signal = row['signal_peptide_annotated'] == 'True'
        ec = row['glycosidase_ec_3_2_1'] == 'True'
        output.append({
            'reference': ref['reference'], 'expected_role': ref['role'],
            'accession': row['uniprot_accession'], 'gene_name': row['gene_name'],
            'match_method': match_method, 'match_status': 'matched',
            'signal_annotated': signal, 'complete_ec_3_2_1': ec,
            'primary_screen_pass': signal and ec,
            'ec_numbers': row['ec_numbers'],
            'interpretation': (
                'Passes both annotation criteria' if signal and ec else
                'Not selected by strict EC/signal annotation criteria; does not disprove secretion'
            ),
        })
    return output


def main():
    with (RESULTS / 'protein_properties.csv').open(encoding='utf-8', newline='') as file:
        rows = list(csv.DictReader(file))
    checks = reference_results(rows)
    with (RESULTS / 'reference_check.csv').open('w', newline='', encoding='utf-8') as file:
        writer = csv.DictWriter(file, fieldnames=list(checks[0]))
        writer.writeheader()
        writer.writerows(checks)
    summary = {
        'reference_count': len(checks),
        'resolved': sum(r['match_status'] == 'matched' for r in checks),
        'passed_primary_screen': sum(r['primary_screen_pass'] is True for r in checks),
        'not_found_or_ambiguous': sum(r['match_status'] != 'matched' for r in checks),
        'results': checks,
        'note': 'Targeted illustrative inspection, not independent sensitivity validation.',
    }
    (RESULTS / 'reference_summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print('REFERENCE CHECK')
    for row in checks:
        print(row['reference'], '|', row['match_status'], '|', row['accession'],
              '| primary screen:', row['primary_screen_pass'])
    print('Saved results/reference_check.csv')


if __name__ == '__main__':
    main()
