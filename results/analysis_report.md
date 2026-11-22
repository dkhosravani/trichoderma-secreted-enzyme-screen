# Annotation evidence and follow-up analysis

These counts come from the downloaded UniProtKB snapshot. They are
annotation observations, not experimental proof of extracellular secretion.

- Proteome entries: **9,115**
- Signal-feature-annotated entries: **653**
- Complete EC 3.2.1.* entries: **108**
- Entries meeting both criteria: **69**
- EC 3.2.1.* entries without a Signal annotation: **39**
- Candidates in reviewed UniProtKB entries: **4**

## Signal-feature evidence among candidates

| Exclusive evidence-code group | Candidate proteins |
|---|---:|
| ECO:0000255_only | 4 |
| ECO:0000256_only | 65 |

ECO:0000256 denotes a sequence-model match used in an automatic
assertion; ECO:0000255 is a sequence-model match used in a manual
assertion. Neither code demonstrates experimental secretion.
Other or missing evidence codes must be inspected individually.

## Full-length protein property summaries

| Annotation group | Valid pI | Median pI | IQR (Q1–Q3) |
|---|---:|---:|---|
| Signal annotated | 653 | 5.403 | 4.836–6.353 |
| No Signal annotation | 8461 | 6.366 | 5.585–8.719 |

## Post-signal-segment analysis

- Candidates with usable annotated signal coordinates: **69**
- Median (post-signal pI minus full-length pI): **-0.074**

These are theoretical properties of sequences after removal of the
annotated N-terminal Signal region. This is not a verified mature protein:
other processing steps and uncertain cleavage positions may exist.

## Follow-up observations

- Entries with incomplete EC 3.2.1.-: **11**
- Longest annotated protein: **G0RTM6** (20873 aa; Carrier domain-containing protein).

Absence of a Signal feature is not evidence of intracellular localization.
The long-sequence outlier is reported but its gene-model correctness
has not been independently validated.

### Example EC 3.2.1.* entries without signal annotation

| Accession | Protein name |
|---|---|
| G0R6Z3 | mannan endo-1,6-alpha-mannosidase |
| G0R9W7 | Mannan endo-1,6-alpha-mannosidase |
| G0RA95 | Mannosyl-oligosaccharide glucosidase |
| G0RBU1 | Glycogen debranching enzyme |
| G0RC60 | chitinase |
