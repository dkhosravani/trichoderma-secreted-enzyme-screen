# Detailed annotation audit (generated)

Re-run `python src/audit_pipeline.py` to recalculate these observations.

The table below is generated from the CSV outputs by `python src/audit_pipeline.py`.
**Snapshot:** 2026-11-22T22:12:38.267306+00:00; UniProt release: 2026_04.  

| Observation | Count |
|---|---:|
| UniProt entries | 9,115 |
| Signal annotations | 653 (7.16%) |
| No Signal annotation | 8,462 |
| Complete EC `3.2.1.*` | 108 |
| Signal + complete EC candidates | 69 |
| Complete EC but no Signal annotation | 39 |
| Entries with incomplete `3.2.1.-` annotation | 11 |
| Full-length MW/pI calculated | 9,114 |

**Signal-feature evidence among the EC-selected candidates:**

| Evidence code category | Candidates |
|---|---:|
| `ECO:0000255_only` | 4 |
| `ECO:0000256_only` | 65 |
| **Reviewed UniProt entries among candidates** | **4** |

The ECO:0000255 and ECO:0000256 groups describe sequence-model evidence
(manual vs automatic assertions). When these are the only observed groups,
**no experimental evidence codes were found in the extracted Signal features**;
this does not assess all published experimental evidence.

**Full-length sequence summaries (median and IQR):**

| Annotation group | Valid pI | Median pI | pI IQR (Q1–Q3) | Median length (aa) | Length IQR (Q1–Q3) |
|---|---:|---:|---|---:|---|
| Signal annotated | 653 | 5.403 | 4.836–6.353 | 391.0 | 241.0–545.0 |
| No Signal annotation | 8461 | 6.366 | 5.585–8.719 | 412.0 | 264.0–607.0 |

**Exploratory two-sided Mann–Whitney U (full-length pI):** `U=1583906.5`, `p=5.82e-74`, rank-biserial effect `r=-0.4266` (positive values indicate higher ranks in the Signal-annotated group).
The test compares distributions; it does not prove secretion or a causal difference.

**After removal of the annotated signal segment:** 69 candidates had calculable theoretical pI; median post-signal minus full-length pI = **-0.074**. This is not a verified mature-protein measurement.

**TM check:** 4 candidate(s) have an annotated Transmembrane region extending beyond the Signal region. This is a possible membrane-association flag, not a secretion assay.

**Sequence outlier:** `G0RTM6` has 20,873 residues; the predicted gene model was not independently validated.

## EC matches lacking a Signal feature

There were **39** such records. This is an annotation gap/question, not proof of intracellular localization.
Of these, 4 have any UniProt Transmembrane feature.
The complete group is in [`ec_without_signal.csv`](ec_without_signal.csv).

| UniProt accession | Protein name |
|---|---|
| G0R6Z3 | mannan endo-1,6-alpha-mannosidase |
| G0R9W7 | Mannan endo-1,6-alpha-mannosidase |
| G0RA95 | Mannosyl-oligosaccharide glucosidase |
| G0RBU1 | Glycogen debranching enzyme |
| G0RC60 | chitinase |

## Known reference proteins (exploratory check)

The following accessions/gene labels were chosen before inspecting pipeline recovery.
No statement of sensitivity is made from this small reference set.
| Reference | Accession | EC annotation | Meets both strict criteria? |
|---|---|---|---|
| CBH I (Cel7A) | G0RVK1 | 3.2.1.- | no |
| CBH II (Cel6A, UniRef-linked) | G0RB58 | 3.2.1.- | no |
| Endoglucanase II (Cel5A) | G0RB67 | 3.2.1.4 | yes |
| Xylanase II (Xyn2) | G0RUP7 | 3.2.1.8 | yes |
| eIF5A comparison | G0RKY3 | — | no |

The known Cel7A example may be excluded by the EC rule if its accession
contains only an incomplete EC number; this is not a failed secretion experiment.

## Descriptive statistics and effect size

Full-length pI: U=1583906.5, p=5.82e-74, rank-biserial=-0.4266.
Length: U=2537799.5, p=0.000514, rank-biserial=-0.0815.
Tests are exploratory and compare distributions, not causation;
with very large unequal groups, small effects may have very small p-values.

## Boundaries

EC 3.2.1.* refers to an EC-number glycosidase category, not a CAZy GH-family prediction.
The signal feature is not confirmed experimental secretion.
The post-signal sequence is not necessarily the experimentally processed mature protein.
No activity or industrial productivity was measured.
