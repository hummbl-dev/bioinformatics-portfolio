# Data Manifest

No raw data is committed to this repository. Every dataset used by an
artifact is registered here with provenance, license, and retrieval
instructions.

## Selection criteria (flagship)

A flagship dataset must be: (1) public, no access application or IRB;
(2) published with a citable analysis to reproduce; (3) sized for a
workstation, not a cluster; (4) rich enough to show QC → statistics →
biological interpretation end to end.

## Candidates

| ID | Dataset | Source | Size | Why | License |
|---|---|---|---|---|---|
| C1 | TCGA-BRCA Pan-Cancer (mutations, CNA, clinical) | cBioPortal / GDC | ~1k samples | Canonical cancer-genomics cohort; cBioPortal R API makes retrieval scriptable; rich published analyses to reproduce | TCGA publication policy |
| C2 | GEO RNA-seq re-analysis (e.g., GSE96058 SCAN-B breast cohort or comparable) | GEO / ArrayExpress | hundreds–3k samples | Classic DESeq2-style differential expression + enrichment reproduction | Series-dependent |
| C3 | 1000 Genomes Phase 3 (VCF, chr subset) | IGSR / EBI FTP | 2.5k samples, subset to chromosomes | Variant QC + population structure; ideal for D1/D2 pipeline validation | Public |

## Rules

- Retrieval must be scripted (wget/ENA/cBioPortalData/GEOquery) — never
  "download by hand"; the script is the provenance.
- Record checksums where the source publishes them.
- Note the retrieval date and source version in the analysis writeup.
