import random
from algorithms.kmer_index import KmerIndex, de_bruijn_contigs


def _k_unique_seq(length, k, seed=7):
    """Random sequence with no repeated k-mer (unambiguous assembly)."""
    rng = random.Random(seed)
    while True:
        seq = "".join(rng.choice("ACGT") for _ in range(length))
        kmers = {seq[i:i + k] for i in range(len(seq) - k + 1)}
        if len(kmers) == len(seq) - k + 1:
            return seq


def test_kmer_index_locate():
    idx = KmerIndex({"a": "ACGTACGTACGT", "b": "TTTTGGGGCCCC"}, k=4)
    hits = idx.locate("GTACGT")
    assert "a" in hits and "b" not in hits


def test_de_bruijn_reconstructs_kunique():
    k = 15
    seq = _k_unique_seq(200, k)
    reads = [seq[i:i + 50] for i in range(0, len(seq) - 49, 10)]
    contigs = de_bruijn_contigs(reads, k=k)
    assert contigs[0] == seq


def test_de_bruijn_cycle_covers_all_kmers():
    # period-4 sequence collapses to one cyclic contig containing
    # every distinct k-mer — repeats fundamentally cannot resolve
    seq = "ACGT" * 40
    reads = [seq[i:i + 20] for i in range(0, len(seq) - 19, 7)]
    contigs = de_bruijn_contigs(reads, k=4)
    kmers = {seq[i:i + 4] for i in range(len(seq) - 3)}
    assert kmers <= {contigs[0][i:i + 4] for i in range(len(contigs[0]) - 3)}


def test_de_bruijn_branch_split():
    # shared k-mer creates a branch point -> two unitigs
    contigs = de_bruijn_contigs(["AACGTA", "ACGTTT"], k=4)
    assert len(contigs) == 3
    # branch node CGT terminates the incoming unitig and anchors both
    # outgoing unitigs — standard non-branching-path decomposition
    assert {"AACGT", "CGTA", "CGTTT"} <= set(contigs)


def test_kmer_index_accepts_list():
    idx = KmerIndex(["AAACCC", "CCCGGG"], k=3)
    assert idx.query("CCC") == [(0, 3), (1, 0)]
