"""Re-audit downloaded results; report evidence, statistics and caveats.

This script updates the generated key findings in README.md from actual CSV
files. It never assumes counts from a previous UniProt release.
"""
import csv
import json
import re
from collections import Counter
from pathlib import Path
from statistics import median

from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / 'results'
EC_PATTERN = re.compile(r'3\.2\.1\.\d+\Z')
START = '<!-- AUTO_RESULTS_BEGIN -->'
END = '<!-- AUTO_RESULTS_END -->'


def csv_rows(name):
    with (RESULTS / name).open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle))


def read_json(name):
    with (RESULTS / name).open(encoding='utf-8') as handle:
        return json.load(handle)


def percentile(values, fraction):
    """Linear interpolation equivalent to pandas/NumPy defaults."""
    numbers = sorted(float(x) for x in values)
    if not numbers:
        return None
    position = (len(numbers) - 1) * fraction
    low = int(position)
    high = min(low + 1, len(numbers) - 1)
    return round(numbers[low] + (numbers[high] - numbers[low]) * (position - low), 3)


def numeric_values(rows, field):
    return [float(row[field]) for row in rows if row[field] not in ('', None)]


def group_stats(rows, field):
    values = numeric_values(rows, field)
    return {
        'n_valid': len(values),
        'median': round(median(values), 3) if values else None,
        'q1': percentile(values, .25),
        'q3': percentile(values, .75),
    }


def compare_distributions(signal_values, other_values):
    """Two-sided Mann–Whitney and rank-biserial effect (signal minus other).

    Positive effect means higher values/ranks tend to occur in the signal
    group. Does not assume equal shape, causality, or any secretion proof.
    """
    x = list(signal_values)
    y = list(other_values)
    if not x or not y:
        return {'method': 'two-sided_asymptotic', 'n_signal': len(x),
                'n_no_signal': len(y), 'u': None, 'p_value': None,
                'rank_biserial': None}
    result = mannwhitneyu(x, y, alternative='two-sided', method='asymptotic')
    u = float(result.statistic)
    effect = 2 * u / (len(x) * len(y)) - 1
    return {
        'method': 'two-sided_asymptotic', 'n_signal': len(x), 'n_no_signal': len(y),
        'u': round(u, 3), 'p_value': float(result.pvalue),
        'rank_biserial': round(effect, 4),
    }


def fmt_p(p):
    if p is None:
        return 'not computed'
    if p == 0:
        return '<1e-300 (underflow)'
    return f'{p:.3g}'


def fmt_number(value):
    return 'NA' if value is None else str(value)


def ensure(condition, message):
    if not condition:
        raise AssertionError(message)


def update_readme_generated_section(section):
    path = ROOT / 'README.md'
    text = path.read_text(encoding='utf-8')
    ensure(text.count(START) == 1 and text.count(END) == 1,
           'README results markers are missing or duplicated')
    before, rest = text.split(START, 1)
    _, after = rest.split(END, 1)
    path.write_text(before + START + '\n' + section.rstrip() + '\n' + END + after,
                    encoding='utf-8')


def md_safe(value):
    return str(value).replace('|', '\\|').replace('\n', ' ')


def run_audit():
    all_rows = csv_rows('protein_properties.csv')
    candidates = csv_rows('candidates.csv')
    shortlist = csv_rows('candidate_shortlist.csv')
    metadata = read_json('metadata.json')
    annotation = read_json('annotation_summary.json')
    properties = read_json('property_summary.json')
    selection = read_json('shortlist_metadata.json')

    ids = [r['uniprot_accession'] for r in all_rows]
    candidate_ids = [r['uniprot_accession'] for r in candidates]
    shortlist_ids = [r['uniprot_accession'] for r in shortlist]
    ensure(len(ids) == len(set(ids)) == metadata['downloaded_count'], 'Record uniqueness/count')
    ensure(len(all_rows) == annotation['total_proteins'], 'Total annotation count')
    ensure(len(candidate_ids) == len(set(candidate_ids)), 'Candidate duplicates')
    ensure(len(candidates) == annotation['signal_and_glycosidase'] == properties['candidate_count'],
           'Candidate counts disagree')
    ensure(len(shortlist) == selection['shortlist_count'], 'Shortlist count')
    ensure(len(shortlist_ids) == len(set(shortlist_ids)), 'Shortlist duplicates')
    ensure(set(shortlist_ids) <= set(candidate_ids), 'Shortlist contains noncandidate')

    signal = [r for r in all_rows if r['signal_peptide_annotated'] == 'True']
    no_signal = [r for r in all_rows if r['signal_peptide_annotated'] == 'False']
    ec = [r for r in all_rows if r['glycosidase_ec_3_2_1'] == 'True']
    ec_no_signal = [r for r in ec if r['signal_peptide_annotated'] == 'False']
    ec_signal = [r for r in ec if r['signal_peptide_annotated'] == 'True']
    ensure(len(signal) == annotation['annotated_signal_peptide'], 'Signal count')
    ensure(len(no_signal) == annotation['no_signal_peptide_annotation'], 'No signal count')
    ensure(len(ec) == annotation['glycosidase_ec_3_2_1'], 'EC count')
    ensure(len(ec_signal) == len(candidates), 'Combined selection count')
    ensure(len(ec_no_signal) == annotation['glycosidase_without_signal_annotation'], 'EC sans signal')
    ensure(len(ec) == len(ec_signal) + len(ec_no_signal), 'EC partition')
    ensure(sum(r['incomplete_ec_3_2_1'] == 'True' for r in all_rows)
           == annotation['incomplete_ec_3_2_1_entries'], 'Incomplete EC count')
    ensure(sum(r['property_status'] == 'ok' for r in all_rows)
           == properties['properties_calculated'], 'Full property count')
    ensure(sum(r['post_signal_status'] == 'ok' for r in all_rows)
           == properties['post_signal_properties_calculated'], 'Post signal count')
    ensure(sum(r['property_status'] == 'ok' for r in candidates)
           == properties['candidates_with_properties'], 'Candidate property count')
    for r in candidates:
        ensure(r['signal_peptide_annotated'] == 'True' and
               r['glycosidase_ec_3_2_1'] == 'True', 'Candidate inclusion')
    for r in shortlist:
        ensure(EC_PATTERN.fullmatch(r['representative_ec']), 'Invalid shortlist EC')
        ensure(r['representative_ec'] in r['ec_numbers'].split('; '), 'Representative EC mismatch')
    ensure(len({r['representative_ec'] for r in shortlist}) >= min(
           len(shortlist), selection['distinct_ec_categories']),
           'Shortlist failed to cover eligible EC categories')

    groups = Counter(r['signal_evidence_group'] for r in candidates)
    codes = Counter(
        code for r in candidates for code in r['signal_evidence_codes'].split('; ')
        if code
    )
    ensure(sum(groups.values()) == len(candidates), 'Signal evidence partition')
    post = [r for r in candidates if r['post_signal_status'] == 'ok']
    shifts = [float(r['post_signal_pi']) - float(r['calculated_pi']) for r in post]
    longest = max(all_rows, key=lambda r: int(r['length_aa']) if r['length_aa'] else -1)

    stats_pi_signal = group_stats(signal, 'calculated_pi')
    stats_pi_other = group_stats(no_signal, 'calculated_pi')
    stats_length_signal = group_stats(signal, 'length_aa')
    stats_length_other = group_stats(no_signal, 'length_aa')
    mw_pi = compare_distributions(numeric_values(signal, 'calculated_pi'),
                                  numeric_values(no_signal, 'calculated_pi'))
    mw_length = compare_distributions(numeric_values(signal, 'length_aa'),
                                      numeric_values(no_signal, 'length_aa'))

    details = {
        'total_entries': len(all_rows),
        'retrieval_timestamp_utc': metadata.get('downloaded_at_utc'),
        'uniprot_release': metadata.get('uniprot_release'),
        'signal_annotated': len(signal),
        'complete_ec_3_2_1': len(ec),
        'candidate_entries': len(candidates),
        'candidate_reviewed': sum(r['reviewed'] == 'True' for r in candidates),
        'candidate_signal_evidence_groups': dict(sorted(groups.items())),
        'candidate_signal_evidence_codes': dict(sorted(codes.items())),
        'ec_3_2_1_without_signal_annotation': len(ec_no_signal),
        'ec_3_2_1_without_signal_examples': [
            {'accession': r['uniprot_accession'], 'protein_name': r['protein_name']}
            for r in sorted(ec_no_signal, key=lambda r: r['uniprot_accession'])[:5]
        ],
        'ec_no_signal_transmembrane_feature_count': sum(
            r.get('transmembrane_annotated') == 'True' for r in ec_no_signal
        ),
        'candidate_tm_outside_signal': sum(
            r.get('tm_outside_signal_status') == 'tm_annotation_outside_signal_region'
            for r in candidates
        ),
        'incomplete_ec_3_2_1_entries': annotation['incomplete_ec_3_2_1_entries'],
        'full_length_properties_calculated': properties['properties_calculated'],
        'candidates_with_post_signal_properties': len(post),
        'candidate_post_signal_pi_shift_median': round(median(shifts), 3) if shifts else None,
        'signal_full_length_pi': stats_pi_signal,
        'no_signal_full_length_pi': stats_pi_other,
        'signal_length_aa': stats_length_signal,
        'no_signal_length_aa': stats_length_other,
        'mann_whitney_full_length_pi': mw_pi,
        'mann_whitney_length': mw_length,
        'longest_record': {
            'accession': longest['uniprot_accession'],
            'length_aa': int(longest['length_aa']),
            'protein_name': longest['protein_name'],
        },
    }
    (RESULTS / 'analysis_details.json').write_text(json.dumps(details, indent=2) + '\n',
                                                   encoding='utf-8')
    with (RESULTS / 'ec_without_signal.csv').open('w', newline='', encoding='utf-8') as file:
        writer = csv.DictWriter(file, fieldnames=list(ec_no_signal[0]) if ec_no_signal else list(all_rows[0]))
        writer.writeheader()
        writer.writerows(ec_no_signal)

    references = read_json('reference_summary.json') if (RESULTS / 'reference_summary.json').exists() else None
    reference_lines = []
    if references:
        for item in references['results']:
            passed = item['primary_screen_pass']
            pass_str = ('yes' if passed is True else 'no' if passed is False else 'not resolved')
            reference_lines.append(
                f"| {md_safe(item['reference'])} | {item['accession'] or '—'} | "
                f"{md_safe(item['ec_numbers']) or '—'} | {pass_str} |"
            )

    # This block is refreshed automatically using actual observations.
    section = '\n'.join([
        'The table below is generated from the CSV outputs by `python src/audit_pipeline.py`.',
        f"**Snapshot:** {details['retrieval_timestamp_utc'] or 'see metadata.json'}; "
        f"UniProt release: {details['uniprot_release'] or 'not reported'}.  ",
        '',
        '| Observation | Count |',
        '|---|---:|',
        f"| UniProt entries | {len(all_rows):,} |",
        f"| Signal annotations | {len(signal):,} ({100 * len(signal) / len(all_rows):.2f}%) |",
        f"| No Signal annotation | {len(no_signal):,} |",
        f"| Complete EC `3.2.1.*` | {len(ec):,} |",
        f"| Signal + complete EC candidates | {len(candidates):,} |",
        f"| Complete EC but no Signal annotation | {len(ec_no_signal):,} |",
        f"| Entries with incomplete `3.2.1.-` annotation | {details['incomplete_ec_3_2_1_entries']:,} |",
        f"| Full-length MW/pI calculated | {properties['properties_calculated']:,} |",
        '',
        '**Signal-feature evidence among the EC-selected candidates:**',
        '',
        '| Evidence code category | Candidates |',
        '|---|---:|',
        *[f'| `{md_safe(key)}` | {value} |' for key, value in sorted(groups.items())],
        f"| **Reviewed UniProt entries among candidates** | **{details['candidate_reviewed']}** |",
        '',
        'The ECO:0000255 and ECO:0000256 groups describe sequence-model evidence',
        '(manual vs automatic assertions). When these are the only observed groups,',
        '**no experimental evidence codes were found in the extracted Signal features**;',
        'this does not assess all published experimental evidence.',
        '',
        '**Full-length sequence summaries (median and IQR):**',
        '',
        '| Annotation group | Valid pI | Median pI | pI IQR (Q1–Q3) | Median length (aa) | Length IQR (Q1–Q3) |',
        '|---|---:|---:|---|---:|---|',
        f"| Signal annotated | {stats_pi_signal['n_valid']} | {fmt_number(stats_pi_signal['median'])} | "
        f"{fmt_number(stats_pi_signal['q1'])}–{fmt_number(stats_pi_signal['q3'])} | "
        f"{fmt_number(stats_length_signal['median'])} | {fmt_number(stats_length_signal['q1'])}–{fmt_number(stats_length_signal['q3'])} |",
        f"| No Signal annotation | {stats_pi_other['n_valid']} | {fmt_number(stats_pi_other['median'])} | "
        f"{fmt_number(stats_pi_other['q1'])}–{fmt_number(stats_pi_other['q3'])} | "
        f"{fmt_number(stats_length_other['median'])} | {fmt_number(stats_length_other['q1'])}–{fmt_number(stats_length_other['q3'])} |",
        '',
        f"**Exploratory two-sided Mann–Whitney U (full-length pI):** "
        f"`U={fmt_number(mw_pi['u'])}`, `p={fmt_p(mw_pi['p_value'])}`, "
        f"rank-biserial effect `r={fmt_number(mw_pi['rank_biserial'])}` "
        '(positive values indicate higher ranks in the Signal-annotated group).',
        'The test compares distributions; it does not prove secretion or a causal difference.',
        '',
        f"**After removal of the annotated signal segment:** "
        f"{len(post)} candidates had calculable theoretical pI; "
        f"median post-signal minus full-length pI = "
        f"**{fmt_number(details['candidate_post_signal_pi_shift_median'])}**. "
        'This is not a verified mature-protein measurement.',
        '',
        f"**TM check:** {details['candidate_tm_outside_signal']} candidate(s) have an "
        'annotated Transmembrane region extending beyond the Signal region. '
        'This is a possible membrane-association flag, not a secretion assay.',
        '',
        f"**Sequence outlier:** `{longest['uniprot_accession']}` has "
        f"{int(longest['length_aa']):,} residues; the predicted gene model was not independently validated.",
    ])
    update_readme_generated_section(section)

    report = [
        '# Detailed annotation audit (generated)',
        '',
        'Re-run `python src/audit_pipeline.py` to recalculate these observations.',
        '',
        section,
        '',
        '## EC matches lacking a Signal feature',
        '',
        f"There were **{len(ec_no_signal)}** such records. "
        'This is an annotation gap/question, not proof of intracellular localization.',
        f"Of these, {details['ec_no_signal_transmembrane_feature_count']} have any UniProt Transmembrane feature.",
        'The complete group is in [`ec_without_signal.csv`](ec_without_signal.csv).',
        '',
        '| UniProt accession | Protein name |',
        '|---|---|',
        *[f"| {item['accession']} | {md_safe(item['protein_name'])} |"
          for item in details['ec_3_2_1_without_signal_examples']],
        '',
        '## Known reference proteins (exploratory check)',
        '',
        'The following accessions/gene labels were chosen before inspecting pipeline recovery.',
        'No statement of sensitivity is made from this small reference set.',
        '| Reference | Accession | EC annotation | Meets both strict criteria? |',
        '|---|---|---|---|',
        *(reference_lines or ['| Reference check not run yet | — | — | — |']),
        '',
        'The known Cel7A example may be excluded by the EC rule if its accession',
        'contains only an incomplete EC number; this is not a failed secretion experiment.',
        '',
        '## Descriptive statistics and effect size',
        '',
        f"Full-length pI: U={fmt_number(mw_pi['u'])}, p={fmt_p(mw_pi['p_value'])}, "
        f"rank-biserial={fmt_number(mw_pi['rank_biserial'])}.",
        f"Length: U={fmt_number(mw_length['u'])}, p={fmt_p(mw_length['p_value'])}, "
        f"rank-biserial={fmt_number(mw_length['rank_biserial'])}.",
        'Tests are exploratory and compare distributions, not causation;',
        'with very large unequal groups, small effects may have very small p-values.',
        '',
        '## Boundaries',
        '',
        'EC 3.2.1.* refers to an EC-number glycosidase category, not a CAZy GH-family prediction.',
        'The signal feature is not confirmed experimental secretion.',
        'The post-signal sequence is not necessarily the experimentally processed mature protein.',
        'No activity or industrial productivity was measured.',
    ]
    (RESULTS / 'analysis_report.md').write_text('\n'.join(report) + '\n', encoding='utf-8')
    print('FINAL AUDIT PASSED')
    print('Candidate evidence groups:', dict(groups))
    print('EC glycosidases without Signal:', len(ec_no_signal))
    print('pI Mann–Whitney:', fmt_p(mw_pi['p_value']), 'rank-biserial:', mw_pi['rank_biserial'])
    print('Median post-signal pI change:', details['candidate_post_signal_pi_shift_median'])
    print('README key findings generated from real data')
    return details


def main():
    run_audit()


if __name__ == '__main__':
    main()
