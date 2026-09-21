# Subtype-enriched genomic alterations in TCGA breast cancer: a fully scripted reproduction and triangulation across statistical and machine-learning arms

**Author:** Reuben P. Bowlby
**Status:** working draft for bioRxiv · 2026-09-21
**Code + data:** `analysis/flagship-reproduction/` (this repository)

## Abstract

**Background.** The 2012 TCGA breast-cancer study reported that
driver mutation frequencies differ by PAM50 intrinsic subtype — TP53
enriched in basal-like, PIK3CA in luminal A, GATA3 essentially absent
from basal — and that HER2-enriched tumors are defined by ERBB2
amplification rather than mutation. Whether these findings reproduce
under a fully scripted, parameter-light pipeline on the later
PanCancer Atlas release is an open practical question for teaching
reproducible cancer genomics.

**Methods.** We queried the cBioPortal public API for the
`brca_tcga_pan_can_atlas_2018` study, restricted to the 981 patients
with PAM50 calls, and tested non-silent mutation status for TP53,
PIK3CA, and GATA3 plus GISTIC discrete copy-number calls for six
genes against subtype using two-sided Fisher exact tests with
Benjamini–Hochberg correction (pure-standard-library Python). We then
asked whether a supervised classifier could recover the same
structure: a 46-gene binary mutation panel, dummy/logistic/random-
forest models, stratified 5-fold out-of-fold evaluation.

**Results.** All primary associations reproduced with matching
direction and near-matching magnitude: TP53-mutant basal 89.5%
(OR 27.7, p=9.4×10⁻⁶¹), PIK3CA-mutant luminal A 47.5% (OR 3.25),
GATA3 absent from basal (0/171). The mutation-only classifier
achieved macro-F1 0.40 (dummy 0.13), recovering Basal (F1 0.66) and
LumA (0.72) but failing on Her2 (0.18) — a *correct* failure, since
Her2 is defined by ERBB2 amplification, which the CNA arm confirmed
(70.5% of Her2, OR 33.6). A flat learning curve (n=196→981) shows the
classifier is feature-limited, not data-limited.

**Conclusions.** The subtype-enrichment structure of TCGA-BRCA
reproduces cleanly on the PanCancer Atlas release with a ~500-line
dependency-light pipeline, and the agreement between three
independent arms (Fisher tests, supervised ML, CNA tests) constitutes
a stronger verification than any single re-analysis. The pipeline
itself is the deliverable: every number in this paper regenerates
from two scripts and a public API.

## 1. Introduction

Replication studies are rare in computational genomics partly because
the cost of re-deriving published cohorts is high and partly because
"close enough" reproductions are hard to certify. The TCGA 2012
breast paper is an ideal target: the claims are discrete
(gene × subtype associations), the cohort is public, and a later
reprocessed release (PanCancer Atlas 2018) provides a natural
out-of-sample check.

The question posed here is deliberately narrow: **do the published
subtype-enrichment associations reproduce under a scripted pipeline
on the Atlas cohort, and does a supervised learner independently
recover the same biological structure?** The second half matters —
convergent failure modes (the classifier failing exactly where the
biology says mutation shouldn't carry signal) are stronger evidence
than convergent successes alone.

## 2. Methods

### 2.1 Cohort and data retrieval

cBioPortal public REST API, study `brca_tcga_pan_can_atlas_2018`.
1084 samples retrieved; PAM50 subtype is a patient-level clinical
attribute (`SUBTYPE`, values `BRCA_*`), present for 981 patients
(LumA 499, LumB 197, Basal 171, Her2 78, Normal 36). Samples without
subtype calls were excluded rather than imputed. Retrieval script:
`scripts/fetch_tcga_brca.py`.

### 2.2 Alteration calls

Mutations: all mutation records for TP53, PIK3CA, GATA3 plus a
43-gene curated driver panel (total 46), excluding Silent, Intron,
UTR, Flank, IGR, and Splice_Region types (Splice_Site retained).
Copy number: GISTIC discrete calls for ERBB2, CCND1, MYC, PTEN,
FGFR1, FGFR3; AMP = value ≥ +2, HOMDEL ≤ −2. Patients with multiple
samples were collapsed to any-call.

### 2.3 Statistics

Per gene × subtype: 2×2 contingency (altered/wt × subtype/rest),
two-sided Fisher exact computed from the hypergeometric
distribution in pure Python (no scipy dependency for the flagship
arm); Benjamini–Hochberg q-values across all tests.

### 2.4 Classifier

Feature matrix: 981 patients × 46 binary non-silent mutation
indicators. Models: DummyClassifier (most-frequent), multinomial
logistic regression (L2, class-balanced), random forest (500 trees,
balanced subsample). Evaluation: StratifiedKFold(5, shuffle,
seed=42) with `cross_val_predict` — every reported prediction is
out-of-fold. Primary metric: macro-F1 (accuracy reported but
de-emphasized due to class imbalance). Learning curve: logistic
macro-F1 at train fractions 0.2–1.0.

## 3. Results

### 3.1 Mutation associations reproduce

| Gene | Subtype | Freq | OR | p | TCGA 2012 |
|---|---|---|---|---|---|
| TP53 | Basal | 89.5% | 27.7 | 9.4e-61 | ~84% ✓ |
| TP53 | Her2 | 70.5% | 5.1 | 3.6e-11 | ~72–75% ✓ |
| TP53 | LumA | 10.6% | 0.08 | 2.3e-63 | ~12–29% ✓ |
| PIK3CA | LumA | 47.5% | 3.25 | 1.8e-17 | ~45% ✓ |
| PIK3CA | Basal | 7.0% | 0.11 | 2.8e-20 | ~7–9% ✓ |
| GATA3 | Basal | 0/171 | 0.0 | 3.2e-11 | near-absent ✓ |
| GATA3 | LumB | 20.3% | 2.2 | 2.4e-04 | luminal ✓ |

### 3.2 The classifier fails exactly where biology predicts

RF macro-F1 0.40 (dummy 0.135). Per-class F1: Basal 0.66, LumA 0.72,
LumB 0.26, Her2 0.18, Normal 0.20. Confusion concentrated in
LumA↔LumB and Her2→LumA. Feature importances (full fit): TP53 0.19,
PIK3CA 0.09, GATA3 0.06 — the same genes the Fisher tests ranked.
Learning curve flat (~0.39) across 196→981 training samples.

### 3.3 CNA arm closes the loop

ERBB2-AMP: Her2 70.5% (OR 33.6); MYC-AMP: Basal 35.7% (OR 4.3);
PTEN-HOMDEL: Basal 16.4% (OR 7.0); CCND1-AMP: LumB 25.9% (OR 2.4).
All match published Atlas patterns and explain the ML failure:
Her2 requires copy-number features the mutation panel lacks.

## 4. Discussion

Three observations generalize beyond this dataset:

1. **Convergent failure is evidence.** The classifier's inability to
   learn Her2 from mutations is not a defect — it is the same fact
   the CNA arm confirms positively (ERBB2-AMP OR 33.6). A reviewer
   can triangulate the biology from three independent directions.
2. **Flat learning curves diagnose feature limitations.** The
   ~0.39 plateau from n=196 to n=981 says more patients cannot fix
   a feature space that lacks expression/CNA — the correct next
   feature, identified *before* fitting, was ERBB2 amplification.
3. **Reproduction cost is now trivial.** The entire analysis —
   retrieval, statistics, ML — is three scripts against a public API
   with no credentials and no local database. The barrier to
   replication studies is no longer infrastructure.

## 5. Limitations

- Atlas cohort (981 subtyped) ≠ 2012 discovery set (~507); frequency
  deltas of a few points are expected and observed.
- Non-silent filter approximates but does not replicate the paper's
  exact mutation-type handling.
- Normal-like (n=36) underpowered; likely a purity artifact biologically.
- No purity/stage/ER-status adjustment — matching the descriptive
  scope of the original claim.
- Subtype calls used as-is from cBioPortal; no re-calling from
  expression data (documented as future work).

## References

1. The Cancer Genome Atlas Network. Comprehensive molecular portraits
   of human breast tumours. *Nature* 490, 61–70 (2012).
2. Hoadley KA et al. Cell-of-Origin Patterns Dominate the Molecular
   Classification of 10,000 Tumors from 33 Types of Cancer. *Cell*
   173, 291–304 (2018). (PanCancer Atlas)
3. Cerami E et al. The cBio Cancer Genomics Portal. *Cancer Discov*
   2, 401–404 (2012).
4. Parker JS et al. Supervised Risk Predictor of Breast Cancer Based
   on Intrinsic Subtypes. *J Clin Oncol* 27, 1160–1167 (2009). (PAM50)
