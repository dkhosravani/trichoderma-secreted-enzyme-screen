"""Create reproducible descriptive figures from validated CSV outputs."""
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / 'results'
FIGURES = ROOT / 'figures'
FIGURES.mkdir(exist_ok=True)


def main():
    data = pd.read_csv(RESULTS / 'protein_properties.csv')
    candidates = pd.read_csv(RESULTS / 'candidates.csv')
    signal = data['signal_peptide_annotated'].eq(True)
    length = pd.to_numeric(data['length_aa'], errors='coerce')
    pi = pd.to_numeric(data['calculated_pi'], errors='coerce')
    labels = [f'Signal annotated (n={signal.sum():,})',
              f'No signal annotation (n={(~signal).sum():,})']

    # Annotation frequencies: annotations are NOT measured secretion.
    counts = [int(signal.sum()), int((~signal).sum())]
    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(labels, counts)
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2, count + max(counts) * .02,
                f'{count:,} ({100 * count / len(data):.2f}%)', ha='center')
    ax.set_ylim(0, max(counts) * 1.16)
    ax.set_ylabel('UniProt entries')
    ax.set_title('Signal peptide feature annotations')
    fig.tight_layout()
    fig.savefig(FIGURES / 'signal_annotation_fraction.png', dpi=250)
    plt.close(fig)

    for values, ylabel, title, filename in [
        (length, 'Full-length protein (aa)', 'Protein length by annotation group',
         'protein_length_comparison.png'),
        (pi, 'Theoretical full-length pI', 'Protein pI by annotation group',
         'theoretical_pi_comparison.png'),
    ]:
        groups = [values[signal].dropna(), values[~signal].dropna()]
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.boxplot(groups, showfliers=False)
        ax.set_xticks([1, 2], labels)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        fig.text(.5, .01, 'Outliers hidden; whiskers = 1.5 × IQR', ha='center', fontsize=9)
        fig.tight_layout(rect=[0, .04, 1, 1])
        fig.savefig(FIGURES / filename, dpi=250)
        plt.close(fig)

    # Evidence codes in the final candidate pool: mutually exclusive categories.
    groups = candidates['signal_evidence_group'].fillna('no_evidence_code').value_counts()
    fig, ax = plt.subplots(figsize=(9, 5))
    bars = ax.barh(list(groups.index), list(groups.values))
    ax.bar_label(bars, padding=3)
    ax.invert_yaxis()
    ax.set_xlim(0, max(groups.values) * 1.2)
    ax.set_xlabel('Candidate proteins')
    ax.set_title('Evidence types for signal annotations in EC-selected candidates')
    fig.tight_layout()
    fig.savefig(FIGURES / 'candidate_signal_evidence.png', dpi=250)
    plt.close(fig)

    # Paired difference for protein segments after the annotated signal region.
    clipped = candidates.loc[candidates['post_signal_status'].eq('ok')].copy()
    if not clipped.empty:
        difference = (pd.to_numeric(clipped['post_signal_pi'])
                      - pd.to_numeric(clipped['calculated_pi']))
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(difference, bins=min(15, max(5, len(difference) // 4)))
        ax.axvline(0, linestyle='--', color='black', linewidth=1)
        ax.set_xlabel('Post-signal pI minus full-length pI')
        ax.set_ylabel('Candidate count')
        ax.set_title(f'pI shift after removal of annotated signal region (n={len(difference)})')
        fig.tight_layout()
        fig.savefig(FIGURES / 'post_signal_pi_shift.png', dpi=250)
        plt.close(fig)

    print('Figures saved in figures/:')
    for file in sorted(FIGURES.glob('*.png')):
        print(' -', file.name)


if __name__ == '__main__':
    main()
