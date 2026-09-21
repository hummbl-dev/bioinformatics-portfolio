# Build Roadmap

Ordered to produce citable evidence before Dec–Jan application deadlines.
Each artifact = public repo content + writeup + (where noted) external signal.

## Phase 1 — Foundations (weeks 1–6)

- [x] `algorithms/`: NW/SW alignment, seed-and-extend w/ X-drop, k-mer
      index + de Bruijn unitigs, profile HMM (Viterbi + forward) —
      stdlib-only, 23 pytest tests green, benchmarks in RESULTS.md
- [~] `pipelines/ngs-variant/`: Nextflow DSL2 written (FastQC→fastp→
      BWA→bcftools→snpEff→MultiQC) + synthetic test data (20kb ref,
      40 planted SNVs, truth.vcf) — pending execution in conda/Docker env
- [x] `env/`: environment.yml + Dockerfile + CI workflow

## Phase 2 — Flagship (weeks 7–14)

- [x] `analysis/flagship-reproduction/`: TCGA-BRCA via cBioPortal —
      **reproduced** TCGA 2012 subtype-enrichment (stats arm) +
      **extension arm** (GISTIC CNA: ERBB2-Her2 OR=33.6)
- [x] Manuscript draft: `writeup/manuscript.md` — journal format,
      3-arm triangulation story, ready for preprint polish
- [ ] Preprint: post to bioRxiv (operator action — needs account)
- [ ] 3-aim research proposal (F31/general-exam format, ~6 pages)

## Phase 3 — Breadth + external validation (weeks 15–20)

- [x] `ml/`: PAM50 classifier, 46-gene mutation panel, honest
      baselines + biological error analysis (ml/RESULTS.md)
- [ ] 2–3 merged PRs into nf-core modules / Bioconductor / Biopython /
      scikit-bio (external validation no self-published repo can fake)
- [ ] `docs/literature-notes/`: 4–6 journal-club-style paper writeups

## Ongoing

- [ ] Every artifact: README, LICENSE headers, CI green, plain-English abstract
- [ ] Statement-of-purpose raw material accumulates in `docs/` as real writeups
