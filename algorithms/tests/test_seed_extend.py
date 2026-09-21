from algorithms.seed_extend import SeedIndex, find_hsps


def test_seed_index_positions():
    idx = SeedIndex({"s1": "ACGTACGT"}, k=3)
    assert idx.query("ACG") == [("s1", 0), ("s1", 4)]
    assert idx.query("TTT") == []


def test_hsp_exact_hit():
    db = {"subj": "GGGGGACGTACGTAGGGGG"}
    idx = SeedIndex(db, k=3)
    hsps = find_hsps("q1", "TTTACGTACGTATTT", idx, match=5, mismatch=-4)
    assert hsps
    top = hsps[0]
    assert top.subject_id == "subj"
    assert top.aligned_q == top.aligned_s  # ungapped, exact middle
    assert "ACGTACGTA" in top.aligned_s


def test_hsp_extension_past_seed():
    db = {"subj": "AAACCCGGGTTT"}
    idx = SeedIndex(db, k=2)
    hsps = find_hsps("q", "CCCGGG", idx, match=5, mismatch=-4)
    assert hsps[0].s_start == 3 and hsps[0].s_end == 9
    assert hsps[0].score == 6 * 5


def test_hsp_dedup_and_ranking():
    db = {"a": "ACGTACGT", "b": "ACGTACGA"}
    idx = SeedIndex(db, k=2)
    hsps = find_hsps("q", "ACGTACGT", idx, match=5, mismatch=-4)
    subs = {h.subject_id for h in hsps}
    assert "a" in subs
    assert hsps[0].subject_id == "a"  # better score wins


def test_no_hsp_below_threshold():
    idx = SeedIndex({"s": "GGGGGGGG"}, k=3)
    assert find_hsps("q", "AAAAAA", idx) == []
