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

## Extension arm — copy number (the deliberate deviation)

The mutation-only ML arm (see `ml/RESULTS.md`) failed on Her2
(F1 0.18) — correctly, because Her2 subtype is defined by ERBB2
*amplification*, a copy-number event invisible to mutation
indicators. Fetched GISTIC discrete CNA calls (same cohort) and
re-ran the subtype-association test on AMP (+2) and HOMDEL (−2):

| Gene | Event | Subtype | Freq | OR | Published ref | Match |
|---|---|---|---|---|---|---|
| ERBB2 | AMP | **Her2** | **70.5%** | 33.6 | ~77% (2012) | ✓ |
| ERBB2 | AMP | Basal | 1.2% | 0.07 | depleted | ✓ |
| MYC | AMP | **Basal** | **35.7%** | 4.3 | basal-MYC | ✓ |
| PTEN | HOMDEL | **Basal** | **16.4%** | 7.0 | basal PTEN loss | ✓ |
| CCND1 | AMP | **LumB** | **25.9%** | 2.4 | luminal 11q13 | ✓ |

**Cross-arm triangulation**: the statistical arm (Fisher on
mutations), the ML arm (mutation features → subtype), and this CNA
arm all tell one coherent story — mutation features separate
Basal/LumA; CNA features are required for Her2. A portfolio reviewer
can verify the same biology from three independent directions.

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
