# Flagship: Reproduction Study (Domains D2 + D4 + D5 + D8)

The capstone artifact: a documented re-analysis of a published omics
finding on public data, ending in a bioRxiv preprint. This is the
"rotation zero" — evidence of independent research capability, which
is what PhD committees actually weigh.

## Status: first result in hand

**Reproduced** the subtype-enrichment of driver mutations in TCGA-BRCA
(PanCancer Atlas 2018, n=981): TP53 basal-enriched (89.5%, OR 27.7),
PIK3CA LumA-enriched (47.5%, OR 3.25), GATA3 basal-depleted (0/171).
All directions and magnitudes match TCGA 2012 *Nature* 490:61-70.
See [`REPRODUCTION.md`](REPRODUCTION.md).

## Design (continuing)

1. ~~Pick the dataset~~ — **done**: C1, TCGA-BRCA via cBioPortal
2. ~~State the claim~~ — **done**: subtype-enriched drivers (above)
3. ~~Reproduce the primary result~~ — **done**: see REPRODUCTION.md
4. **Extend** with one deliberate deviation: e.g., add CNA-driven
   calls (GISTIC amp/del), or test robustness to the silent-mutation
   filter definition
5. **Write it up**: journal-format manuscript in `writeup/`, then
   bioRxiv submission

## Honesty contract

A failed reproduction reported rigorously is worth more than a forced
match. Deviations, data differences, and method ambiguity in the
source paper are documented in REPRODUCTION.md — that is the skill
being demonstrated.

