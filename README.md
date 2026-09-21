# Bioinformatics Proof-of-Work Portfolio

Public evidence of bioinformatics and computational biology capability, built
before graduate applications. Every artifact is reproducible: pinned
environments, documented data provenance, tests, and honest limitations.

## Why this exists

Graduate admissions committees evaluate *demonstrated* work. This repository
implements, in public, the competency core shared by the top U.S. bioinformatics
programs (algorithms for sequence data, NGS pipelines, applied biostatistics,
ML for biology, and an independent research-style study). See
[`docs/capability-map.md`](docs/capability-map.md) for the program-by-program
requirement mapping.

## Structure

| Directory | Contents | Domain |
|---|---|---|
| `algorithms/` | Sequence algorithms implemented from scratch (alignment, k-mer indexing, HMM profiles) with tests and benchmarks | Algorithms |
| `pipelines/ngs-variant/` | End-to-end NGS workflow on public data: QC → alignment → variant calling → annotation | Pipelines |
| `analysis/flagship-reproduction/` | Reproduction of a published omics analysis on public data, with a full writeup | Research |
| `ml/` | Machine learning on biological data: baselines, error analysis, biological interpretation | ML |
| `data/` | Dataset manifests — accessions, licenses, retrieval instructions (no data committed) | Provenance |
| `docs/` | Capability map, literature notes, writeups | Documentation |
| `env/` | Conda environment + Dockerfile for bit-level reproducibility | Reproducibility |

## Standards

- Every pipeline: `environment.yml` + container + CI check
- Every analysis: data manifest entry, methods section, limitations section
- Every claim: reproducible from a clean checkout

## License

MIT — see [`LICENSE`](LICENSE).
