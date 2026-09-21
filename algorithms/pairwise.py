"""Pairwise sequence alignment: Needleman-Wunsch (global) and
Smith-Waterman (local), with traceback.

Scoring is either match/mismatch or an explicit substitution matrix
passed as a dict {(a, b): score} (order-insensitive lookups fall back
to the transposed key, then to `mismatch`).
"""


def _score(ca, cb, match, mismatch, subst):
    if subst is not None:
        return subst.get((ca, cb), subst.get((cb, ca), mismatch))
    return match if ca == cb else mismatch


def _traceback(s, i, j, a, b, match, mismatch, gap, subst):
    aa, bb = [], []
    while i > 0 and j > 0 and s[i][j] > 0:
        if s[i][j] == s[i - 1][j - 1] + _score(a[i - 1], b[j - 1], match, mismatch, subst):
            aa.append(a[i - 1]); bb.append(b[j - 1]); i -= 1; j -= 1
        elif s[i][j] == s[i - 1][j] + gap:
            aa.append(a[i - 1]); bb.append("-"); i -= 1
        else:
            aa.append("-"); bb.append(b[j - 1]); j -= 1
    return "".join(reversed(aa)), "".join(reversed(bb))


def needleman_wunsch(a, b, match=1, mismatch=-1, gap=-1, subst=None):
    """Global alignment. Returns (score, aligned_a, aligned_b)."""
    n, m = len(a), len(b)
    s = [[0.0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        s[i][0] = i * gap
    for j in range(1, m + 1):
        s[0][j] = j * gap
    for i in range(1, n + 1):
        ai = a[i - 1]
        row, prev = s[i], s[i - 1]
        for j in range(1, m + 1):
            row[j] = max(
                prev[j - 1] + _score(ai, b[j - 1], match, mismatch, subst),
                prev[j] + gap,
                row[j - 1] + gap,
            )
    i, j = n, m
    aa, bb = [], []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and s[i][j] == s[i - 1][j - 1] + _score(a[i - 1], b[j - 1], match, mismatch, subst):
            aa.append(a[i - 1]); bb.append(b[j - 1]); i -= 1; j -= 1
        elif i > 0 and s[i][j] == s[i - 1][j] + gap:
            aa.append(a[i - 1]); bb.append("-"); i -= 1
        else:
            aa.append("-"); bb.append(b[j - 1]); j -= 1
    return s[n][m], "".join(reversed(aa)), "".join(reversed(bb))


def smith_waterman(a, b, match=2, mismatch=-1, gap=-1, subst=None):
    """Local alignment. Returns (score, aligned_a, aligned_b, (i_end, j_end))."""
    n, m = len(a), len(b)
    s = [[0.0] * (m + 1) for _ in range(n + 1)]
    best, bi, bj = 0.0, 0, 0
    for i in range(1, n + 1):
        ai = a[i - 1]
        row, prev = s[i], s[i - 1]
        for j in range(1, m + 1):
            row[j] = max(
                prev[j - 1] + _score(ai, b[j - 1], match, mismatch, subst),
                prev[j] + gap,
                row[j - 1] + gap,
                0.0,
            )
            if row[j] > best:
                best, bi, bj = row[j], i, j
    aa, bb = _traceback(s, bi, bj, a, b, match, mismatch, gap, subst)
    return best, aa, bb, (bi, bj)
