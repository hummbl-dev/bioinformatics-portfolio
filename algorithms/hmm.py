"""Profile HMM: build from a multiple alignment, score sequences with
Viterbi and the forward algorithm (Durbin et al., ch. 5 topology).

States: M_j (match), I_j (insert), D_j (delete) for j = 1..L, plus a
silent Begin. End is reached implicitly via a termination transition
score folded into the final-column max.
"""

import math

NEG = float("-inf")


def _log(x):
    return math.log(x) if x > 0 else NEG


def _allowed_transitions(L):
    """Legal successors per state in the Durbin profile-HMM topology."""
    allowed = {"B": ["M1", "I0", "D1", "E"],
               "I0": ["M1", "I0", "D1"]}
    for j in range(1, L):
        for s in "MID":
            allowed[f"{s}{j}"] = [f"M{j + 1}", f"I{j}", f"D{j + 1}"]
    for s in "MID":
        allowed[f"{s}{L}"] = [f"I{L}", "E"]
    return allowed


def _logsumexp(xs):
    xs = [x for x in xs if x > NEG]
    if not xs:
        return NEG
    m = max(xs)
    return m + math.log(sum(math.exp(x - m) for x in xs))


class ProfileHMM:
    def __init__(self, length, alphabet, emit_m, emit_i, trans):
        self.L = length
        self.alphabet = alphabet
        self.emit_m = emit_m          # [{sym: prob}] length L
        self.emit_i = emit_i          # [{sym: prob}] length L+1 (I_0..I_L)
        self.trans = trans            # {(from_state, to_state): prob}

    @classmethod
    def from_alignment(cls, aln, alphabet="ACGT", pseudocount=0.5,
                       trans_pseudocount=0.1, gap_thresh=0.5):
        """Build a profile HMM from a gap-containing alignment."""
        aln = [s.upper() for s in aln]
        width = len(aln[0])
        # match columns: gap fraction below threshold
        is_match = [sum(s[c] == "-" for s in aln) / len(aln) < gap_thresh
                    for c in range(width)]
        match_cols = [c for c in range(width) if is_match[c]]
        L = len(match_cols)

        # emission counts
        emit_m = [{a: pseudocount for a in alphabet} for _ in range(L)]
        emit_i = [{a: pseudocount for a in alphabet} for _ in range(L + 1)]
        # transition counts between state types
        tcount = {}
        for s in aln:
            prev = "B"
            j = 0
            for c in range(width):
                if is_match[c]:
                    j += 1
                    if s[c] == "-":
                        st = f"D{j}"
                    else:
                        st = f"M{j}"
                        emit_m[j - 1][s[c]] += 1
                else:
                    if s[c] == "-":
                        continue
                    st = f"I{j}"
                    emit_i[j][s[c]] += 1
                tcount[(prev, st)] = tcount.get((prev, st), 0) + 1
                prev = st
            tcount[(prev, "E")] = tcount.get((prev, "E"), 0) + 1

        def norm(d):
            tot = sum(d.values())
            return {k: v / tot for k, v in d.items()}

        # transitions: observed counts + pseudocount over the legal
        # successor set so unobserved-but-legal paths stay reachable
        trans = {}
        for f, tos in _allowed_transitions(L).items():
            counts = {t: tcount.get((f, t), 0.0) + trans_pseudocount
                      for t in tos}
            for t, p in norm(counts).items():
                trans[(f, t)] = p
        return cls(L, alphabet, [norm(e) for e in emit_m],
                   [norm(e) for e in emit_i], trans)

    def _t(self, f, t):
        return _log(self.trans.get((f, t), 0.0))

    def _e(self, table, j, sym):
        return _log(table[j].get(sym, 0.0))

    def viterbi(self, seq):
        """Log-probability of the most likely state path."""
        n, L = len(seq), self.L
        vm = [[NEG] * (L + 2) for _ in range(n + 1)]
        vi = [[NEG] * (L + 2) for _ in range(n + 1)]
        vd = [[NEG] * (L + 2) for _ in range(n + 1)]
        vm[0][0] = vi[0][0] = vd[0][0] = 0.0  # Begin
        for i in range(1, n + 1):
            x = seq[i - 1]
            for j in range(0, L + 1):
                e_i = self._e(self.emit_i, j, x)
                p_m = f"M{j}" if j else "B"
                vi[i][j] = e_i + max(
                    vm[i - 1][j] + self._t(p_m, f"I{j}"),
                    vi[i - 1][j] + self._t(f"I{j}", f"I{j}"),
                    vd[i - 1][j] + self._t(f"D{j}", f"I{j}"),
                )
                if j == 0:
                    continue
                e_m = self._e(self.emit_m, j - 1, x)
                p_m, p_i, p_d = (f"M{j-1}" if j > 1 else "B"), f"I{j-1}", f"D{j-1}"
                vm[i][j] = e_m + max(
                    vm[i - 1][j - 1] + self._t(p_m, f"M{j}"),
                    vi[i - 1][j - 1] + self._t(p_i, f"M{j}"),
                    vd[i - 1][j - 1] + self._t(p_d, f"M{j}"),
                )
            for j in range(1, L + 1):
                p_m, p_i, p_d = (f"M{j-1}" if j > 1 else "B"), f"I{j-1}", f"D{j-1}"
                vd[i][j] = max(
                    vm[i][j - 1] + self._t(p_m, f"D{j}"),
                    vi[i][j - 1] + self._t(p_i, f"D{j}"),
                    vd[i][j - 1] + self._t(p_d, f"D{j}"),
                )
        return max(
            vm[n][L] + self._t(f"M{L}", "E"),
            vi[n][L] + self._t(f"I{L}", "E"),
            vd[n][L] + self._t(f"D{L}", "E"),
        )

    def forward(self, seq):
        """Log-probability of the sequence summed over all paths."""
        n, L = len(seq), self.L
        fm = [[NEG] * (L + 2) for _ in range(n + 1)]
        fi = [[NEG] * (L + 2) for _ in range(n + 1)]
        fd = [[NEG] * (L + 2) for _ in range(n + 1)]
        fm[0][0] = fi[0][0] = fd[0][0] = 0.0
        for i in range(1, n + 1):
            x = seq[i - 1]
            for j in range(0, L + 1):
                if j == 0:
                    srcs = [("B", fm[i - 1][0]), ("I0", fi[i - 1][0]),
                            ("D0", fd[i - 1][0])]
                else:
                    srcs = [(f"M{j}", fm[i - 1][j]), (f"I{j}", fi[i - 1][j]),
                            (f"D{j}", fd[i - 1][j])]
                fi[i][j] = self._e(self.emit_i, j, x) + _logsumexp(
                    [v + self._t(s, f"I{j}") for s, v in srcs])
                if j == 0:
                    continue
                srcs = [(f"M{j-1}" if j > 1 else "B", fm[i - 1][j - 1]),
                        (f"I{j-1}", fi[i - 1][j - 1]),
                        (f"D{j-1}", fd[i - 1][j - 1])]
                fm[i][j] = self._e(self.emit_m, j - 1, x) + _logsumexp(
                    [v + self._t(s, f"M{j}") for s, v in srcs])
            for j in range(1, L + 1):
                srcs = [(f"M{j-1}" if j > 1 else "B", fm[i][j - 1]),
                        (f"I{j-1}", fi[i][j - 1]),
                        (f"D{j-1}", fd[i][j - 1])]
                fd[i][j] = _logsumexp(
                    [v + self._t(s, f"D{j}") for s, v in srcs])
        return _logsumexp([
            fm[n][L] + self._t(f"M{L}", "E"),
            fi[n][L] + self._t(f"I{L}", "E"),
            fd[n][L] + self._t(f"D{L}", "E"),
        ])
