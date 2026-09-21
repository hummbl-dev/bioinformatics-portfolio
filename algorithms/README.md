# Algorithms (Domain D1)

From-scratch implementations of the algorithms every bioinformatics
program teaches in its first-year core — the "could you build it, not
just call it" evidence.

## Scope

| Component | Algorithm(s) | Validation |
|---|---|---|
| Pairwise alignment | Needleman-Wunsch, Smith-Waterman | Score agreement vs `Biopython.Align` on synthetic + real pairs |
| Seed-and-extend search | BLAST-style k-mer seeding, extension | Recall vs BLAST on a small query/db set |
| k-mer index | De Bruijn / FM-index style structures | Query correctness + memory/time benchmark |
| Profile HMM | Build, Viterbi, forward | Agreement vs `hmmer` on a protein family subset |

## Standard

- Pure Python, stdlib only; no wrapping existing tools for the core logic
- `pytest` suite: correctness fixtures + property tests (scoring
  symmetry, gap penalties, edge cases)
- `benchmarks/`: runtime vs reference implementations, honest about
  where naive implementations lose
- A short `RESULTS.md` per component: what it does, complexity,
  where it diverges from the reference and why
