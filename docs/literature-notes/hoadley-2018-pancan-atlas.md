# Note — Hoadley et al. (2018), *Cell* 173:291–304

**Why this paper:** it produced the actual dataset I analyzed —
`brca_tcga_pan_can_atlas_2018` on cBioPortal is this release.

## What they did

Reprocessed all ~10,000 TCGA tumors (33 types) through uniform
pipelines; showed molecular classification groups by cell-of-origin
more than by tissue site (e.g., squamous cancers across organs
cluster together).

## Relevance to my work

- Uniform reprocessing = my "out-of-sample" testbed: the 2012 BRCA
  claims re-tested on the Atlas cohort differ in denominator and
  platform harmonization — a real replication check, not a re-run.
- iCluster/COCA subtypes vs PAM50: my study used PAM50 (`SUBTYPE`
  attribute); COCA labels would be a second subtype axis worth
  cross-tabulating (future).

## Strengths

- Harmonization removed per-study batch artifacts — the reason
  cross-study pooling is legitimate here.
- Public, versioned, API-accessible — the entire reproducibility
  story depends on this release existing.

## Caveats

- "Subtype" in Atlas outputs has multiple sources (PAM50 for BRCA,
  other classifiers elsewhere) — my Aim-1 multi-cohort plan must
  handle heterogeneous subtype axes honestly.
- GISTIC calls are consensus-level; threshold choices (±2 as
  AMP/HOMDEL) are conventions, not biology — documented.

## Takeaway

The Atlas is the right replication substrate: same biology, different
processing pipeline → surviving associations are the robust ones.
