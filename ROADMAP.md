# Build Roadmap

Ordered to produce citable evidence before Dec–Jan application deadlines.
Each artifact = public repo content + writeup + (where noted) external signal.

## Phase 1 — Foundations (weeks 1–6)

- [ ] `algorithms/`: Smith-Waterman + BLAST-style seed-and-extend, k-mer index,
      profile HMM build/search — Python implementations w/ pytest + benchmark
      vs Biopython/pysam
- [ ] `pipelines/ngs-variant/`: FASTQ QC (FastQC/MultiQC) → trimming → BWA-MEM
      → GATK-style calling → annotation on 1000 Genomes subset
- [ ] `env/`: environment.yml + Dockerfile + CI workflow running tests

## Phase 2 — Flagship (weeks 7–14)

- [ ] `analysis/flagship-reproduction/`: pick dataset from `data/MANIFEST.md`
      candidates → full re-analysis → methods/results/limitations writeup
- [ ] Preprint: post writeup to bioRxiv (free, timestamped, citable)
- [ ] 3-aim research proposal (F31/general-exam format, ~6 pages) — matches UW
      Genome Sciences general exam and UNC written qual structure

## Phase 3 — Breadth + external validation (weeks 15–20)

- [ ] `ml/`: one ML-on-omics project (variant effect or expression classifier)
      with honest baselines + error analysis
- [ ] 2–3 merged PRs into nf-core modules / Bioconductor / Biopython /
      scikit-bio (external validation no self-published repo can fake)
- [ ] `docs/literature-notes/`: 4–6 journal-club-style paper writeups

## Ongoing

- [ ] Every artifact: README, LICENSE headers, CI green, plain-English abstract
- [ ] Statement-of-purpose raw material accumulates in `docs/` as real writeups
