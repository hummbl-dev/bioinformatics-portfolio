# Note — TCGA Network (2012), *Nature* 490:61–70

**Claim under test in my reproduction:** driver mutation frequency
differs by PAM50 subtype.

## What they did

~507 breast tumors assayed across six platforms (exome seq, SNP6 CNA,
methylation, mRNA, miRNA, RPPA). PAM50 expression subtypes used as the
organizing axis; mutation rates and subtype-specific driver
associations reported.

## Key findings (as used in my work)

- TP53 mutated ~80%+ of basal-like, minority of luminal A.
- PIK3CA most frequent in luminal A (~45%).
- GATA3 luminal-enriched, essentially absent in basal.
- HER2-enriched subtype ≠ all clinically HER2+; defined largely by
  ERBB2 amplicon signature.
- Basal-like shares serous-ovarian-like genomic features.

## Strengths

- Multi-platform integration was the real contribution — mutations
  alone would have missed the CNA-defined biology (my ML arm
  rediscovered this the hard way, which is the point).
- Subtype axis made heterogeneous cohort interpretable.

## Weaknesses / caveats

- Discovery cohort ~507 — smaller than later releases; frequencies
  carry binomial noise (my Atlas n=981 gives ~±3% CIs at p=0.5).
- Subtype calls from expression assays with batch effects; some
  discordance with later PAM50 call sets.
- "Normal-like" later questioned as low-purity artifact — consistent
  with my underpowered n=36 result.

## What I verified independently

All three mutation associations + the ERBB2-AMP Her2 definition —
direction and approximate magnitude — on the Atlas 2018 release.
See `analysis/flagship-reproduction/REPRODUCTION.md`.

## Open questions this raises for me

- Do the claims hold per-histology and per-stage strata?
- How much of the subtype–mutation association is mediated by ER
  status rather than subtype per se? (decomposition analysis — future)
