# Note — Cerami et al. (2012), *Cancer Discov* 2:401–404

**Why this paper:** the infrastructure enabling the entire flagship
study — a public REST API over curated cancer-genomics datasets.

## What it is

Web portal + REST API exposing TCGA/other studies as studies →
molecular profiles → sample lists → alteration records, plus clinical
attributes. No credentials needed for public data.

## What I used in practice

- `/studies/{id}/molecular-profiles` — discover mutation/CNA profiles
- `/sample-lists` — cohort definitions (`_all`, per-platform)
- `mutations/fetch`, `discrete-copy-number/fetch` — bulk alteration
  retrieval by Entrez IDs + sample list
- `clinical-data/fetch` — subtype attribute (patient-level for PAM50)

## Practical lessons (worth more than the paper's claims)

- **Attribute levels matter**: `SUBTYPE` exists at SAMPLE level but
  is populated at PATIENT level — my first fetch returned 0 rows;
  the fix was a patient-level fallback. This class of silent-empty
  response is the top reproducibility hazard in API-driven analysis.
- **Label vocabularies aren't documented inline**: subtype values are
  `BRCA_*`-prefixed; assume nothing, census values before coding.
- **Sample↔patient mapping**: TCGA barcode first-12-chars rule is
  reliable but must be asserted, not assumed.
- **Flaky 5xx**: the public API rate-limits/wobbles — retry/backoff
  is required infrastructure, not polish.

## Takeaway

API-first access converts cancer-genomics replication from a data-
engineering project into an afternoon script — but the failure modes
move from "can't get data" to "silently wrong joins." Logging every
denominator (as my fetch script does) is the defense.
