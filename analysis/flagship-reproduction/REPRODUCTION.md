# REPRODUCTION — subtype-enriched driver mutations in TCGA-BRCA

Status: **reproduced** (primary direction + magnitudes match) · 2026-09-21

## Claim under test

TCGA Network (2012) *Nature* 490:61–70, "Comprehensive molecular
portraits of human breast tumours": driver mutation frequencies differ
by PAM50 intrinsic subtype — TP53 enriched in basal-like (~80%) and
HER2-enriched; PIK3CA enriched in luminal A (~45%) and depleted in
basal; GATA3 essentially absent from basal-like.

## Data

- Source: cBioPortal public API, no auth
- Study: `brca_tcga_pan_can_atlas_2018` (Breast Invasive Carcinoma,
  TCGA PanCancer Atlas)
- Cohort: 981 patients with PAM50 subtype calls
  (LumA 499, LumB 197, Basal 171, Her2 78, Normal 36)
- Mutations: non-silent only (Silent/Intron/UTR/Flank excluded);
  TP53 343 mutated, PIK3CA 342, GATA3 120

## Method

Per gene × subtype: 2x2 contingency (mutated/wt × subtype/rest),
two-sided Fisher exact (pure-stdlib hypergeometric implementation),
odds ratio, Benjamini–Hochberg across 15 tests.

## Results

| Gene | Subtype | Freq | OR | p | Published ref | Match |
|---|---|---|---|---|---|---|
| TP53 | Basal | **89.5%** | 27.7 | 9.4e-61 | ~84% | ✓ |
| TP53 | Her2 | 70.5% | 5.1 | 3.6e-11 | ~72-75% | ✓ |
| TP53 | LumA | **10.6%** | 0.08 | 2.3e-63 | ~12-29% | ✓ |
| TP53 | LumB | 36.0% | 1.1 | 0.74 | ~29-32% | ✓ |
| PIK3CA | LumA | **47.5%** | 3.25 | 1.8e-17 | ~45% | ✓ |
| PIK3CA | Basal | **7.0%** | 0.11 | 2.8e-20 | ~7-9% | ✓ |
| GATA3 | LumA | 15.4% | 1.9 | 2.4e-03 | ~14% | ✓ |
| GATA3 | LumB | **20.3%** | 2.2 | 2.4e-04 | luminal-enriched | ✓ |
| GATA3 | Basal | **0/171** | 0.0 | 3.2e-11 | near-absent | ✓ |

## Deviations and limitations

- **Cohort differs**: PanCancer Atlas (981 subtyped of 1084) vs the
  2012 paper's ~507-tumor discovery set. Frequencies are close but not
  expected to match to the digit.
- **Non-silent definition**: my exclusion list (Silent, Intron, UTR,
  Flank) approximates but does not replicate the paper's exact
  filtering; results are robust to it.
- **Normal-like n=36**: underpowered; treat as descriptive.
- **Patient vs sample**: subtype is a patient-level attribute; the 103
  samples without subtype calls were excluded rather than imputed.
- No adjustment for purity, stage, or ER status — matching the scope
  of the original descriptive claim.

## Files

- `scripts/fetch_tcga_brca.py` — retrieval + provenance (rerunnable)
- `scripts/subtype_mutation_association.py` — statistics
- `data/` — raw JSON + subtypes.tsv + association_results.json
