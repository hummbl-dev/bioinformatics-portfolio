# Flagship: Reproduction Study (Domains D2 + D4 + D5 + D8)

The capstone artifact: a documented re-analysis of a published omics
study on public data, ending in a bioRxiv preprint. This is the
"rotation zero" — evidence of independent research capability, which
is what PhD committees actually weigh.

## Design

1. **Pick the dataset** from `data/MANIFEST.md` candidates (default:
   C1 TCGA-BRCA via cBioPortal — canonical, scriptable, richly
   published).
2. **State the claim being reproduced** — e.g., a reported mutation–
   subtype association, a differential-expression signature, or a
   survival association from the original paper.
3. **Reproduce the primary result** with the published methods, then
   **extend** with one deliberate deviation (updated reference release,
   alternative normalization, or an added confounder) and document the
   effect.
4. **Write it up**: intro, methods, results, limitations, reproduction
   notes — journal format, 6–10 pages, then bioRxiv.

## Deliverables

- `notebooks/` or `src/` — the executable analysis
- `writeup/` — manuscript source (Markdown/Quarto → PDF)
- `REPRODUCTION.md` — what matched, what didn't, why
- BioRxiv submission record

## Honesty contract

A failed reproduction reported rigorously is worth more than a forced
match. Document deviations, data differences, and method ambiguity in
the source paper — that is the skill being demonstrated.
