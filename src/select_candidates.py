"""Select a reproducible EC-diverse shortlist, not a performance ranking."""
import csv
import json
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / 'results'
TARGET = 10
EC_PATTERN = re.compile(r'3\.2\.1\.\d+\Z')


def valid_ec_numbers(row):
    return sorted({part.strip() for part in row['ec_numbers'].split(';')
                   if EC_PATTERN.fullmatch(part.strip())})


def annotation_priority(row):
    return (
        row['reviewed'] != 'True',
        not bool(row['gene_name'].strip()),
        not bool(row['protein_name'].strip()),
        row['uniprot_accession'],
    )


def select_shortlist(rows, target=TARGET):
    eligible = [r for r in rows if r['property_status'] == 'ok'
                and r['signal_peptide_annotated'] == 'True'
                and r['glycosidase_ec_3_2_1'] == 'True'
                and valid_ec_numbers(r)]
    ec_counts = Counter(ec for row in eligible for ec in valid_ec_numbers(row))
    ec_order = sorted(ec_counts, key=lambda ec: (-ec_counts[ec], ec))
    selected, used = [], set()

    for ec in ec_order:
        choices = [r for r in eligible if ec in valid_ec_numbers(r)
                   and r['uniprot_accession'] not in used]
        if choices:
            chosen = min(choices, key=annotation_priority)
            selected.append({**chosen, 'representative_ec': ec,
                             'selection_reason': 'EC category diversity'})
            used.add(chosen['uniprot_accession'])
        if len(selected) == target:
            break

    for row in sorted(eligible, key=annotation_priority):
        if len(selected) == target:
            break
        if row['uniprot_accession'] not in used:
            selected.append({**row, 'representative_ec': valid_ec_numbers(row)[0],
                             'selection_reason': 'Additional representative'})
            used.add(row['uniprot_accession'])
    return selected, len(eligible), len(ec_counts)


def main():
    with (RESULTS / 'candidates.csv').open(encoding='utf-8', newline='') as handle:
        candidates = list(csv.DictReader(handle))
    with (RESULTS / 'property_summary.json').open(encoding='utf-8') as handle:
        properties = json.load(handle)
    if len(candidates) != properties['candidate_count']:
        raise ValueError('Candidate count mismatch')

    selected, eligible_count, ec_categories = select_shortlist(candidates)
    if not selected:
        raise ValueError('No eligible candidates')
    with (RESULTS / 'candidate_shortlist.csv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(selected[0]))
        writer.writeheader()
        writer.writerows(selected)

    metadata = {
        'eligible_candidates': eligible_count,
        'distinct_ec_categories': ec_categories,
        'shortlist_count': len(selected),
        'target_count': TARGET,
        'selection_method': 'EC diversity by category frequency, then review/name completeness and accession',
        'not_an_activity_or_yield_ranking': True,
    }
    with (RESULTS / 'shortlist_metadata.json').open('w', encoding='utf-8') as handle:
        json.dump(metadata, handle, indent=2)
    print(f'SHORTLIST: {len(selected)} from {eligible_count}; EC categories: {ec_categories}')
    for row in selected:
        print(row['uniprot_accession'], '|', row['protein_name'], '|', row['representative_ec'])


if __name__ == '__main__':
    main()
