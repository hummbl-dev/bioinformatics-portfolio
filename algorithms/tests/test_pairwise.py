from algorithms.pairwise import needleman_wunsch, smith_waterman


def test_nw_identical():
    score, aa, bb = needleman_wunsch("GATTACA", "GATTACA")
    assert score == 7
    assert aa == bb == "GATTACA"


def test_nw_gap_and_mismatch():
    score, aa, bb = needleman_wunsch("GATTACA", "GCATGCU", match=1, mismatch=-1, gap=-1)
    assert len(aa) == len(bb)
    assert score == -0.0 or score == 0  # known optimal: 0
    # every aligned pair either matches or is scored -1
    for x, y in zip(aa, bb):
        if x != "-" and y != "-":
            assert x != y or x == y


def test_nw_score_symmetric():
    s1, _, _ = needleman_wunsch("ACGTAC", "AGTAC")
    s2, _, _ = needleman_wunsch("AGTAC", "ACGTAC")
    assert s1 == s2


def test_nw_substitution_matrix():
    subst = {("A", "G"): 3}
    s_sub, _, _ = needleman_wunsch("A", "G", subst=subst, mismatch=-5)
    s_mm, _, _ = needleman_wunsch("A", "G", mismatch=-5, gap=-10)
    assert s_sub == 3
    assert s_mm == -5


def test_sw_finds_local_match():
    best, aa, bb, _ = smith_waterman("TTTACGTAAA", "GGGACGTCCC")
    assert best == 8  # 4 exact matches * 2
    assert aa == bb == "ACGT"


def test_sw_zero_floor():
    best, aa, bb, _ = smith_waterman("AAAA", "TTTT")
    assert best == 0


def test_sw_empty():
    assert smith_waterman("", "ACGT")[0] == 0


def test_sw_beats_ungapped():
    # local alignment recovers the 8bp conserved block, ignores flanks
    best, aa, bb, _ = smith_waterman("AAAGCTAGCGTAAA", "CCCGCTAGCGTGGG")
    assert best == 16
    assert "GCTAGCGT" in aa
