# ML for Biology (Domain D6)

One disciplined machine-learning project on biological data. The bar
is not "a model that works" — it is correct experimental hygiene,
honest baselines, and biological interpretation.

## Status: first run complete

**PAM50 subtype classifier from a 46-gene mutation panel** — 981
TCGA-BRCA patients, stratified 5-fold out-of-fold evaluation.
RF macro-F1 0.40 vs dummy 0.13: real signal on Basal (F1 0.66) and
LumA (0.72), correctly unlearnable on Her2 (defined by ERBB2
*amplicon*, not mutation) and LumB (expression-defined). Flat
learning curve → feature-limited. See [`RESULTS.md`](RESULTS.md) and
`results/` (metrics.json + confusion matrix + importances +
learning curve).

## Candidate problems (pick one)

- Variant effect / pathogenicity classification (ClinVar labels +
  gnomAD features)
- Expression-based subtype or outcome classifier (TCGA/GEO)
- Protein function or localization prediction (UniProt + embeddings)

## Non-negotiables

- **Baselines first**: logistic regression / random forest before any
  neural net — the neural result must beat them to matter
- **Proper splits**: split by gene/family or patient, never random-row
  leakage; document the leakage controls
- **Error analysis**: where does the model fail, and is the failure
  biological or technical?
- **Interpretation**: feature attribution or equivalent — connect model
  behavior back to biology
- `RESULTS.md`: metrics with CIs, baselines table, failure modes,
  limitations
