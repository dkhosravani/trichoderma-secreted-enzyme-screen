"""Rebuild figures from measured annotations and the audit report."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # Works in CI and on machines without a display.
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / 'results'
FIGURES = ROOT / 'figures'
FIGURES.mkdir(exist_ok=True)
TEAL = '#187F83'
GRAY = '#9AA6B2'


def label_for_group(signal, other):
    return [f'Signal feature (n={signal:,})', f'No Signal feature (n={other:,})']


def main():
    data = pd.read_csv(RESULTS / 'protein_properties.csv')
    candidates = pd.read_csv(RESULTS / 'candidates.csv')
    with (RESULTS / 'analysis_details.json').open(encoding='utf-8') as file:
        audit = json.load(file)
    signal = data['signal_peptide_annotated'].astype(str).str.lower().eq('true')
    count_signal, count_other = int(signal.sum()), int((~signal).sum())
    labels = label_for_group(count_signal, count_other)
    length = pd.to_numeric(data['length_aa'], errors='coerce')
    pi = pd.to_numeric(data['calculated_pi'], errors='coerce')

    # Simple annotation count as a supplementary plot.
    fig, ax = plt.subplots(figsize=(9, 5))
    counts = [count_signal, count_other]
    bars = ax.bar(labels, counts, color=[TEAL, GRAY])
    ax.bar_label(bars, labels=[f'{n:,} ({n / len(data):.2%})' for n in counts], padding=5)
    ax.set_ylim(0, max(counts) * 1.18)
    ax.set_title('UniProt Signal feature annotation counts')
    ax.set_ylabel('Protein records')
    fig.tight_layout()
    fig.savefig(FIGURES / 'signal_annotation_fraction.png', dpi=250)
    plt.close(fig)

    # Supplementary boxplots; outliers are not removed from the statistics.
    for values, ylabel, title, name in [
        (length, 'Amino acids (full-length)', 'Protein length by Signal annotation',
         'protein_length_comparison.png'),
        (pi, 'Theoretical pI (full-length)', 'Full-length protein pI by Signal annotation',
         'theoretical_pi_comparison.png'),
    ]:
        groups = [values[signal].dropna(), values[~signal].dropna()]
        fig, ax = plt.subplots(figsize=(9, 5))
        boxes = ax.boxplot(groups, showfliers=False, patch_artist=True)
        ax.set_xticks([1, 2], labels)
        for patch, color in zip(boxes['boxes'], [TEAL, GRAY]):
            patch.set_facecolor(color)
            patch.set_alpha(.5)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        fig.text(.5, .015, 'Outlier markers hidden; whiskers extend to 1.5 × IQR',
                 ha='center', fontsize=9)
        fig.tight_layout(rect=[0, .04, 1, 1])
        fig.savefig(FIGURES / name, dpi=250)
        plt.close(fig)

    # Compare distributions, not raw counts, because groups are very unequal.
    x = pi[signal].dropna().to_numpy(dtype=float)
    y = pi[~signal].dropna().to_numpy(dtype=float)
    if len(x) and len(y):
        edges = np.linspace(min(x.min(), y.min()), max(x.max(), y.max()), 45)
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.hist(y, bins=edges, density=True, alpha=.45, color=GRAY,
                label=labels[1])
        ax.hist(x, bins=edges, density=True, alpha=.55, color=TEAL,
                label=labels[0])
        ax.axvline(float(np.median(x)), linestyle='--', linewidth=1, color=TEAL)
        ax.axvline(float(np.median(y)), linestyle='--', linewidth=1, color='#596778')
        ax.set_title('Distribution of theoretical pI across annotation groups')
        ax.set_ylabel('Density (normalized within each group)')
        ax.set_xlabel('Theoretical pI of full-length protein')
        stats = audit['mann_whitney_full_length_pi']
        p = stats['p_value']
        p_text = ('p<1e-300' if p == 0 else f'p={p:.3g}' if p is not None else 'p=NA')
        ax.text(.99, .97,
                f"Exploratory Mann–Whitney, {p_text}\nrank-biserial r={stats['rank_biserial']}",
                ha='right', va='top', transform=ax.transAxes, fontsize=9,
                bbox={'facecolor': 'white', 'alpha': .75, 'edgecolor': 'none'})
        ax.legend(loc='upper left')
        fig.tight_layout()
        fig.savefig(FIGURES / 'full_length_pi_distribution.png', dpi=250)
        plt.close(fig)

    # Per-candidate evidence categories; mutually exclusive by parser definition.
    evidence = candidates['signal_evidence_group'].fillna('no_evidence_code').value_counts()
    display_names = {
        'ECO:0000256_only': 'ECO:0000256 (automatic model)',
        'ECO:0000255_only': 'ECO:0000255 (manual model)',
        'both_model_codes': 'Both model codes',
        'other_or_mixed_evidence': 'Other / mixed code',
        'no_evidence_code': 'Missing evidence code',
    }
    fig, ax = plt.subplots(figsize=(9, 5))
    bars = ax.barh([display_names.get(x, x) for x in evidence.index],
                   evidence.to_list(), color=TEAL)
    ax.bar_label(bars, padding=4)
    ax.invert_yaxis()
    ax.set_xlim(0, max(evidence) * 1.16)
    ax.set_xlabel('Candidate proteins')
    ax.set_title('Signal-feature evidence among EC-selected candidates')
    fig.text(.5, .02, 'Sequence-model evidence is not an experimental secretion measurement',
             ha='center', fontsize=9)
    fig.tight_layout(rect=[0, .045, 1, 1])
    fig.savefig(FIGURES / 'candidate_signal_evidence.png', dpi=250)
    plt.close(fig)

    # Paired difference after removing the *annotated* signal region.
    post = candidates.loc[candidates['post_signal_status'].eq('ok')]
    if len(post):
        differences = (pd.to_numeric(post['post_signal_pi'])
                       - pd.to_numeric(post['calculated_pi'])).dropna()
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.hist(differences, bins=min(18, max(6, len(differences) // 4)),
                color=TEAL, alpha=.8)
        ax.axvline(0, color='black', linestyle='--', linewidth=1, label='No change')
        median_shift = differences.median()
        ax.axvline(median_shift, color='#B55425', linestyle=':', linewidth=2,
                   label=f'Median shift: {median_shift:.3f}')
        ax.set_xlabel('Post-signal theoretical pI − full-length theoretical pI')
        ax.set_ylabel('Candidate proteins')
        ax.set_title(f'pI change after removing annotated signal region (n={len(differences)})')
        ax.legend()
        fig.text(.5, .015, 'Post-signal is an annotation-based segment, not confirmed mature protein',
                 ha='center', fontsize=9)
        fig.tight_layout(rect=[0, .04, 1, 1])
        fig.savefig(FIGURES / 'post_signal_pi_shift.png', dpi=250)
        plt.close(fig)

    print('Figures regenerated:')
    for path in sorted(FIGURES.glob('*.png')):
        print(' -', path.name)


if __name__ == '__main__':
    main()
