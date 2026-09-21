# Bioinformatics Proof-of-Work Portfolio

Public, reproducible evidence of bioinformatics and computational
biology capability — built before graduate applications. Every number
in this README regenerates from scripts in this repo and public APIs.

## Headline result

**Reproduced the TCGA 2012 subtype–driver associations on the
PanCancer Atlas 2018 cohort (n=981), then verified them two more
ways.** Three independent arms, one biological story:

| Arm | Method | Key finding |
|---|---|---|
| Statistics | Fisher exact + BH, pure stdlib | TP53-basal 89.5% (OR 27.7, p≈1e-60); PIK3CA-LumA 47.5% (OR 3.25); GATA3 absent in basal (0/171) |
| Machine learning | RF, stratified out-of-fold | macro-F1 0.40 — learns Basal/LumA, *correctly fails* on Her2 (F1 0.18): mutations can't see an amplicon |
| Copy number | GISTIC discrete calls | ERBB2-AMP recovers Her2: 70.5% (OR 33.6); MYC-AMP basal 35.7%; PTEN-del basal 16.4% |

→ [`analysis/flagship-reproduction/REPRODUCTION.md`](analysis/flagship-reproduction/REPRODUCTION.md)
· [`ml/RESULTS.md`](ml/RESULTS.md)
· [manuscript draft](analysis/flagship-reproduction/writeup/manuscript.md)

## Quickstart (reproduce the flagship in ~2 min)

```bash
git clone <this-repo> && cd bioinformatics-portfolio
python analysis/flagship-reproduction/scripts/fetch_tcga_brca.py
python analysis/flagship-reproduction/scripts/subtype_mutation_association.py
python analysis/flagship-reproduction/scripts/cna_subtype_association.py
```

Stdlib-only for the flagship arm (Python 3.11+; no packages needed).
ML arm needs `scikit-learn`, `pandas`, `matplotlib`
(`env/environment.yml` pins everything). Algorithm tests:

```bash
python -m pytest algorithms/tests   # 23 tests
```

## What's inside

| Directory | Contents | Status |
|---|---|---|
| `algorithms/` | NW/SW alignment, BLAST-style seed-and-extend, k-mer index, de Bruijn unitigs, profile HMM (Viterbi + forward) — all from scratch, no Bio deps | ✅ 23 tests green, benchmarks in `algorithms/RESULTS.md` |
| `analysis/flagship-reproduction/` | TCGA-BRCA reproduction: fetch → stats → CNA extension → manuscript | ✅ verified vs published values |
| `ml/` | PAM50 classifier from 46-gene mutation panel; baselines, OOF evaluation, error analysis | ✅ `ml/RESULTS.md` + figures |
| `pipelines/ngs-variant/` | Nextflow DSL2 (FastQC→fastp→BWA→bcftools→snpEff→MultiQC) + synthetic test data w/ 40 planted SNVs + truth.vcf | 🟡 written; executes under `env/` (bwa/samtools not on dev host) |
| `docs/proposal/` | F31-format 3-aim proposal extending the preliminary data | ✅ draft |
| `docs/literature-notes/` | Journal-club notes on the papers this work engages | ✅ 4 notes |
| `docs/capability-map.md` | Public mapping of artifacts → top-20 program competencies | ✅ |
| `data/` | Dataset manifests + provenance (raw data fetched by scripts, not committed) | ✅ |
| `env/` | Conda `environment.yml` + Dockerfile + CI | ✅ |

## Standards

- **No claim without a regeneration path**: every number traces to a
  script + public API + committed output.
- **Honest denominators**: 981 subtyped of 1084 samples; patient-level
  attributes; exclusions logged, not silently dropped.
- **Failures are data**: the ML arm's Her2 failure is reported as the
  finding it is, not tuned away.

## What this does *not* claim

Independent computational work only — no wet-lab experience, no
clinical access, no institutional training is implied or asserted.
Subtype calls are used as published (no re-calling). Limitations per
analysis live next to the results.

## License

MIT — see [`LICENSE`](LICENSE).
