"""Audit all observed results and record interpretive summary for the README."""
import csv
import json
import re
from collections import Counter
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / 'results'
EC_PATTERN = re.compile(r'3\.2\.1\.\d+\Z')


def csv_rows(name):
    with (RESULTS / name).open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle))


def read_json(name):
    with (RESULTS / name).open(encoding='utf-8') as handle:
        return json.load(handle)


def percentile(values, percentage):
    """Linear interpolation, matching pandas' default percentile convention."""
    items = sorted(values)
    if not items:
        return None
    i = (len(items) - 1) * percentage
    low, high = int(i), min(int(i) + 1, len(items) - 1)
    return round(items[low] + (items[high] - items[low]) * (i - low), 3)


def group_stats(rows, field):
    values = [float(row[field]) for row in rows if row[field] not in ('', None)]
    return {'n_valid': len(values), 'median': round(median(values), 3) if values else None,
            'q1': percentile(values, .25), 'q3': percentile(values, .75)}


def main():
    all_rows = csv_rows('protein_properties.csv')
    candidates = csv_rows('candidates.csv')
    shortlist = csv_rows('candidate_shortlist.csv')
    metadata = read_json('metadata.json')
    annotation = read_json('annotation_summary.json')
    properties = read_json('property_summary.json')
    selection = read_json('shortlist_metadata.json')

    def check(condition, description):
        if not condition:
            raise AssertionError(description)

    ids = [r['uniprot_accession'] for r in all_rows]
    candidate_ids = [r['uniprot_accession'] for r in candidates]
    shortlist_ids = [r['uniprot_accession'] for r in shortlist]
    check(len(ids) == len(set(ids)) == metadata['downloaded_count'], 'Proteome duplicate/count')
    check(len(all_rows) == annotation['total_proteins'], 'Total count differs')
    check(len(candidates) == annotation['signal_and_glycosidase'] == properties['candidate_count'],
          'Candidate count differs')
    check(len(shortlist) == selection['shortlist_count'], 'Shortlist count differs')
    check(len(shortlist_ids) == len(set(shortlist_ids)), 'Duplicate shortlist accession')
    check(set(shortlist_ids) <= set(candidate_ids), 'Shortlist not subset of candidates')

    signal = [r for r in all_rows if r['signal_peptide_annotated'] == 'True']
    gh = [r for r in all_rows if r['glycosidase_ec_3_2_1'] == 'True']
    signal_gh = [r for r in gh if r['signal_peptide_annotated'] == 'True']
    no_signal_gh = [r for r in gh if r['signal_peptide_annotated'] == 'False']
    check(len(signal) == annotation['annotated_signal_peptide'], 'Signal count differs')
    check(len(gh) == annotation['glycosidase_ec_3_2_1'], 'EC count differs')
    check(len(signal_gh) == len(candidates), 'Candidate classification differs')
    check(len(no_signal_gh) == annotation['glycosidase_without_signal_annotation'],
          'EC without signal count differs')
    check(len(gh) == len(signal_gh) + len(no_signal_gh), 'EC partition invalid')
    check(len([r for r in all_rows if r['property_status'] == 'ok'])
          == properties['properties_calculated'], 'Property status count differs')
    check(len([r for r in all_rows if r['post_signal_status'] == 'ok'])
          == properties['post_signal_properties_calculated'], 'Post-signal count differs')
    check(len([r for r in all_rows if r['incomplete_ec_3_2_1'] == 'True'])
          == annotation['incomplete_ec_3_2_1_entries'], 'Incomplete EC count differs')

    for r in candidates:
        check(r['signal_peptide_annotated'] == 'True' and r['glycosidase_ec_3_2_1'] == 'True',
              'Incorrect candidate filtering')
        check(r['property_status'] == 'ok', 'Missing candidate properties')
    for r in shortlist:
        check(EC_PATTERN.fullmatch(r['representative_ec']) is not None,
              'Incomplete EC in shortlist')
        check(r['representative_ec'] in r['ec_numbers'].split('; '), 'Representative EC not in record')
    check(len({r['representative_ec'] for r in shortlist}) == len(shortlist),
          'Shortlist EC categories not unique')

    by_group = Counter(r['signal_evidence_group'] for r in candidates)
    source_codes = Counter()
    for r in candidates:
        for code in filter(None, r['signal_evidence_codes'].split('; ')):
            source_codes[code] += 1
    check(sum(by_group.values()) == len(candidates), 'Evidence group count differs')

    no_signal = [r for r in all_rows if r['signal_peptide_annotated'] == 'False']
    valid_signal = [r for r in signal if r['post_signal_status'] == 'ok']
    clipped_candidates = [r for r in candidates if r['post_signal_status'] == 'ok']
    shifts = [round(float(r['post_signal_pi']) - float(r['calculated_pi']), 3)
              for r in clipped_candidates]
    longest = max(all_rows, key=lambda r: int(r['length_aa']) if r['length_aa'] else -1)

    details = {
        'total_entries': len(all_rows),
        'signal_annotated': len(signal),
        'ec_3_2_1_entries': len(gh),
        'candidate_entries': len(candidates),
        'candidate_reviewed': sum(r['reviewed'] == 'True' for r in candidates),
        'candidate_signal_evidence_groups': dict(sorted(by_group.items())),
        'candidate_signal_evidence_codes': dict(sorted(source_codes.items())),
        'ec_3_2_1_without_signal_annotation': len(no_signal_gh),
        'ec_3_2_1_without_signal_examples': [
            {'accession': r['uniprot_accession'], 'protein_name': r['protein_name']}
            for r in sorted(no_signal_gh, key=lambda x: x['uniprot_accession'])[:5]
        ],
        'incomplete_ec_3_2_1_entries': annotation['incomplete_ec_3_2_1_entries'],
        'signal_annotated_post_signal_calculated': len(valid_signal),
        'candidate_post_signal_calculated': len(clipped_candidates),
        'candidate_median_post_signal_pi_shift': round(median(shifts), 3) if shifts else None,
        'signal_group_full_length_pi': group_stats(signal, 'calculated_pi'),
        'no_signal_group_full_length_pi': group_stats(no_signal, 'calculated_pi'),
        'signal_group_length_aa': group_stats(signal, 'length_aa'),
        'no_signal_group_length_aa': group_stats(no_signal, 'length_aa'),
        'longest_record': {'accession': longest['uniprot_accession'],
                           'length_aa': longest['length_aa'],
                           'protein_name': longest['protein_name']},
        'interpretation': 'Annotation screen only; signal feature and EC do not verify secretion or activity.',
    }
    with (RESULTS / 'analysis_details.json').open('w', encoding='utf-8') as handle:
        json.dump(details, handle, indent=2)

    # Generate a human-readable report from measured data, not guessed values.
    report = [
        '# Annotation evidence and follow-up analysis',
        '',
        'These counts come from the downloaded UniProtKB snapshot. They are',
        'annotation observations, not experimental proof of extracellular secretion.',
        '',
        f"- Proteome entries: **{len(all_rows):,}**",
        f"- Signal-feature-annotated entries: **{len(signal):,}**",
        f"- Complete EC 3.2.1.* entries: **{len(gh):,}**",
        f"- Entries meeting both criteria: **{len(candidates):,}**",
        f"- EC 3.2.1.* entries without a Signal annotation: **{len(no_signal_gh):,}**",
        f"- Candidates in reviewed UniProtKB entries: **{details['candidate_reviewed']:,}**",
        '',
        '## Signal-feature evidence among candidates',
        '',
        '| Exclusive evidence-code group | Candidate proteins |',
        '|---|---:|',
    ]
    report += [f'| {key} | {count} |' for key, count in sorted(by_group.items())]
    report += [
        '',
        'ECO:0000256 denotes a sequence-model match used in an automatic',
        'assertion; ECO:0000255 is a sequence-model match used in a manual',
        'assertion. Neither code demonstrates experimental secretion.',
        'Other or missing evidence codes must be inspected individually.',
        '',
        '## Full-length protein property summaries',
        '',
        '| Annotation group | Valid pI | Median pI | IQR (Q1–Q3) |',
        '|---|---:|---:|---|',
    ]
    for name, stat in [
        ('Signal annotated', details['signal_group_full_length_pi']),
        ('No Signal annotation', details['no_signal_group_full_length_pi']),
    ]:
        report.append(f"| {name} | {stat['n_valid']} | {stat['median']} | {stat['q1']}–{stat['q3']} |")
    report += [
        '',
        '## Post-signal-segment analysis',
        '',
        f"- Candidates with usable annotated signal coordinates: **{len(clipped_candidates)}**",
        f"- Median (post-signal pI minus full-length pI): **{details['candidate_median_post_signal_pi_shift']}**",
        '',
        'These are theoretical properties of sequences after removal of the',
        'annotated N-terminal Signal region. This is not a verified mature protein:',
        'other processing steps and uncertain cleavage positions may exist.',
        '',
        '## Follow-up observations',
        '',
        f"- Entries with incomplete EC 3.2.1.-: **{annotation['incomplete_ec_3_2_1_entries']}**",
        f"- Longest annotated protein: **{longest['uniprot_accession']}** ({longest['length_aa']} aa; {longest['protein_name']}).",
        '',
        'Absence of a Signal feature is not evidence of intracellular localization.',
        'The long-sequence outlier is reported but its gene-model correctness',
        'has not been independently validated.',
        '',
        '### Example EC 3.2.1.* entries without signal annotation',
        '',
        '| Accession | Protein name |',
        '|---|---|',
    ]
    report += [f"| {r['accession']} | {r['protein_name']} |"
               for r in details['ec_3_2_1_without_signal_examples']]
    (RESULTS / 'analysis_report.md').write_text('\n'.join(report) + '\n', encoding='utf-8')

    print('FINAL AUDIT PASSED')
    print('Evidence groups among candidates:', dict(by_group))
    print('Reviewed candidates:', details['candidate_reviewed'])
    print('EC 3.2.1.* without signal annotation:', len(no_signal_gh))
    print('Candidates with post-signal calculations:', len(clipped_candidates))
    print('Median candidate post-signal pI shift:', details['candidate_median_post_signal_pi_shift'])
    print('Longest sequence:', longest['uniprot_accession'], longest['length_aa'], 'aa')
    print('Saved results/analysis_details.json')


if __name__ == '__main__':
    main()
