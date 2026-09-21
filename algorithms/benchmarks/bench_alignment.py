"""Alignment benchmarks: ours vs Biopython (if installed) + scaling.

Run: python algorithms/benchmarks/bench_alignment.py
"""

import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from algorithms.pairwise import needleman_wunsch, smith_waterman

try:
    from Bio.Align import PairwiseAligner
except ImportError:
    PairwiseAligner = None


def rand_seq(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def time_it(fn, *args):
    t = time.perf_counter()
    out = fn(*args)
    return time.perf_counter() - t, out


def main():
    rng = random.Random(42)
    print(f"{'n':>6} {'NW (s)':>9} {'SW (s)':>9}" +
          ("  BioNW (s)" if PairwiseAligner else "  (biopython not installed)"))
    for n in (100, 200, 400, 800):
        a, b = rand_seq(n, rng), rand_seq(n, rng)
        t_nw, _ = time_it(needleman_wunsch, a, b)
        t_sw, _ = time_it(smith_waterman, a, b)
        line = f"{n:>6} {t_nw:9.3f} {t_sw:9.3f}"
        if PairwiseAligner:
            aln = PairwiseAligner()
            aln.mode = "global"
            t_bio, _ = time_it(aln.align, a, b)
            line += f" {t_bio:10.3f}"
        print(line)

    # score sanity vs Biopython when available
    if PairwiseAligner:
        aln = PairwiseAligner()
        aln.mode = "local"
        aln.match_score, aln.mismatch_score, aln.open_gap_score = 2, -1, -1
        a, b = rand_seq(200, rng), rand_seq(200, rng)
        ours = smith_waterman(a, b)[0]
        theirs = aln.align(a, b)[0].score
        print(f"\nSW sanity: ours={ours} biopython={theirs} "
              f"{'MATCH' if abs(ours - theirs) < 1e-9 else 'DIFFER'}")


if __name__ == "__main__":
    main()
