# ML results — PAM50 subtype from mutation panel

Data: 981 patients × 46 binary non-silent mutation indicators
(curated BRCA driver panel, cBioPortal PanCancer Atlas 2018).
Evaluation: stratified 5-fold `cross_val_predict` — every reported
prediction is out-of-fold.

## Headline

| Model | Accuracy | Macro-F1 |
|---|---|---|
| Dummy (most-frequent) | 0.509 | 0.135 |
| LogReg L2, balanced | 0.481 | 0.394 |
| **RF 500, balanced** | **0.539** | **0.403** |

## Per-class (RF, out-of-fold)

| Subtype | P | R | F1 | n |
|---|---|---|---|---|
| Basal | 0.59 | **0.77** | **0.66** | 171 |
| LumA | 0.80 | 0.65 | **0.72** | 499 |
| LumB | 0.34 | 0.20 | 0.26 | 197 |
| Her2 | 0.18 | 0.18 | 0.18 | 78 |
| Normal | 0.12 | 0.53 | 0.20 | 36 |

## Error analysis — the failures are biologically correct

1. **Basal and LumA are learnable; LumB, Her2, Normal are not — from
   mutations alone.** This is expected, not a model defect:
   - **Her2** subtype is defined by *ERBB2 amplification* (a copy-number
     event), not mutation. Mutation indicators cannot see it → F1 0.18.
     Correct fix: add CNA features (GISTIC calls), not more trees.
   - **LumB vs LumA** differ mainly in *expression/proliferation*, not
     in which driver genes mutate → they confuse into LumA.
   - **Normal-like** (n=36) is underpowered and biologically contested
     as a subtype (often a purity artifact).
2. **Learning curve is flat** (macro-F1 ~0.35–0.40 from n=196→981):
   the model is *feature-limited, not data-limited*. More patients
   won't help; expression or CNA features would.
3. **Feature importances sanity-check against the flagship stats arm:**
   TP53 (0.19), PIK3CA (0.09), GATA3 (0.06), CDH1, PTEN — the same
   genes Fisher's exact flagged as subtype-enriched. Two independent
   methods, one coherent signal.
4. **Dummy baseline reported honestly**: 0.51 accuracy by always
   predicting LumA. Accuracy is a misleading metric here (imbalanced
   classes) — macro-F1 is the primary number, and RF only reaches
   0.40. Reported as-is rather than tuned into a flattering number.

## Consistency with flagship arm

Mutation-only classification recovering Basal (TP53) and LumA
(PIK3CA) but nothing else *is* the reproduction finding restated in
ML terms: those are the only subtypes with a strong mutation signal.
This agreement between the statistical and ML arms is the result.

## Next iteration

- ~~Add GISTIC amp/del calls per gene~~ — **done in the flagship
  extension arm** (`analysis/flagship-reproduction/REPRODUCTION.md`):
  ERBB2-AMP recovers 70.5% of Her2 (OR 33.6), confirming the
  diagnosed failure mode. Next model iteration should concatenate
  mutation + CNA features.
- Merge Normal into analysis as descriptive-only or drop
- Gradient boosting / calibrated probabilities for subtype
  probability outputs
