# Results — algorithms suite

Verified: 2026-09-21 · Python 3.14.6 · `pytest algorithms/tests` → 23 passed

## Pairwise alignment (pairwise.py)

Needleman-Wunsch and Smith-Waterman with full DP + traceback; supports
match/mismatch or arbitrary substitution matrices.

Benchmark (`benchmarks/bench_alignment.py`, 42-seeded random DNA):

| n | NW (s) | SW (s) |
|---|---|---|
| 100 | 0.002 | 0.002 |
| 200 | 0.009 | 0.010 |
| 400 | 0.049 | 0.080 |
| 800 | 0.182 | 0.211 |

Scaling is O(n²) as expected (~4x per doubling). Biopython comparison
arm runs automatically when Biopython is installed — the environment in
`env/environment.yml` includes it.

## Seed-and-extend (seed_extend.py)

BLAST-style: k-mer seed index + X-drop ungapped extension producing
HSPs with score/identity. Tests cover exact hits, extension past the
seed, dedup/ranking, and below-threshold rejection.

## k-mer index + de Bruijn (kmer_index.py)

Exact-match index + unitig assembly. Test coverage: k-unique sequence
reconstruction (200bp reassembled exactly from 50bp reads at k=15),
repeat collapse (period-4 → single cyclic contig covering all k-mers),
and branch splitting (shared k-mer → 3 unitigs).

Honest limitation: assembly is unitig-level — real assemblers add
coverage weighting, error correction, and scaffolding.

## Profile HMM (hmm.py)

Durbin topology (M/I/D + Begin/End), built from a multiple alignment.
Emissions and transitions are pseudocount-smoothed; legal-but-unobserved
transitions stay reachable. Viterbi scores member sequences above random;
forward ≥ Viterbi invariant is tested.
