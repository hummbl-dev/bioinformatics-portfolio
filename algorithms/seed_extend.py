"""BLAST-lite: k-mer seeding + X-drop ungapped extension.

Produces high-scoring segment pairs (HSPs) between a query and a
database of subject sequences — the conceptual core of BLAST.
"""

from dataclasses import dataclass, field


@dataclass
class HSP:
    query_id: str
    subject_id: str
    q_start: int
    q_end: int
    s_start: int
    s_end: int
    score: float
    identity: float
    aligned_q: str = field(repr=False, default="")
    aligned_s: str = field(repr=False, default="")


class SeedIndex:
    """k-mer -> positions index over a sequence database."""

    def __init__(self, db, k=3):
        self.db = dict(db)
        self.k = k
        self.table = {}
        for sid, seq in self.db.items():
            for i in range(len(seq) - k + 1):
                self.table.setdefault(seq[i:i + k], []).append((sid, i))

    def query(self, kmer):
        return self.table.get(kmer, [])


def _extend(q, s, qi, si, k, match, mismatch, xdrop):
    """Extend a seed at (qi, si) in both directions with X-drop cutoff."""
    # left
    score = run = 0.0
    l = best_l = 0
    while qi - l - 1 >= 0 and si - l - 1 >= 0:
        run += match if q[qi - l - 1] == s[si - l - 1] else mismatch
        if run > score:
            score, best_l = run, l + 1
        elif run < score - xdrop:
            break
        l += 1
    left = score
    # right
    score = run = 0.0
    r = best_r = 0
    while qi + k + r < len(q) and si + k + r < len(s):
        run += match if q[qi + k + r] == s[si + k + r] else mismatch
        if run > score:
            score, best_r = run, r + 1
        elif run < score - xdrop:
            break
        r += 1
    right = score
    return best_l, k + best_r, left + match * k + right


def find_hsps(query_id, query, index, match=5, mismatch=-4, xdrop=10,
              min_score=None, top=10):
    """Return up to `top` HSPs sorted by score."""
    k = index.k
    min_score = min_score if min_score is not None else match * k
    hits = {}
    for qi in range(len(query) - k + 1):
        for sid, si in index.table.get(query[qi:qi + k], ()):
            l, span, score = _extend(query, index.db[sid], qi, si, k,
                                     match, mismatch, xdrop)
            if score < min_score:
                continue
            key = (sid, qi - l, si - l)
            prev = hits.get(key)
            if prev is not None and prev[0] >= score:
                continue
            aln_q = query[qi - l:qi - l + span]
            aln_s = index.db[sid][si - l:si - l + span]
            ident = sum(c1 == c2 for c1, c2 in zip(aln_q, aln_s)) / len(aln_q)
            hits[key] = (score, HSP(query_id, sid, qi - l, qi - l + span,
                                    si - l, si - l + span, score, ident,
                                    aln_q, aln_s))
    return [h for _, h in sorted(hits.values(), key=lambda t: -t[0])][:top]
