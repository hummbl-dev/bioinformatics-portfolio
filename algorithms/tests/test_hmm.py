import random
from algorithms.hmm import ProfileHMM

MSA = [
    "ACGT-ACGT",
    "ACGT-ACGT",
    "ACGTAACGT",
    "ACGT-ACGA",
    "TCGT-ACGT",
]


def test_build_length_and_emissions():
    h = ProfileHMM.from_alignment(MSA)
    assert h.L == 8  # 9 cols, one insert col (pos 4 has '-' except one 'A')
    # first column is mostly A, one T
    assert h.emit_m[0]["A"] > h.emit_m[0]["T"]


def test_member_scores_higher_than_random():
    h = ProfileHMM.from_alignment(MSA)
    member = "ACGTAACGT"
    rng = random.Random(1)
    random_seq = "".join(rng.choice("ACGT") for _ in range(9))
    assert h.viterbi(member) > h.viterbi(random_seq)


def test_forward_at_least_viterbi():
    # forward sums over all paths, so it must be >= the best single path
    h = ProfileHMM.from_alignment(MSA)
    seq = "ACGTACGT"
    assert h.forward(seq) >= h.viterbi(seq)


def test_insert_column_tolerated():
    h = ProfileHMM.from_alignment(MSA)
    # sequence with an extra symbol where the MSA has the insert column
    assert h.viterbi("ACGTAAACGT") != float("-inf")


def test_scores_are_negative_log_probs():
    h = ProfileHMM.from_alignment(MSA)
    assert h.viterbi("ACGTACGT") < 0
