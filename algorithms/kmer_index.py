"""k-mer index and de Bruijn graph assembly.

KmerIndex: exact-match lookup of k-mers across a sequence set.
de_bruijn_contigs: unitig extraction from a de Bruijn graph built
from reads — the core of short-read assembly.
"""

from collections import defaultdict


class KmerIndex:
    """k-mer -> sorted list of (seq_id, position) occurrences."""

    def __init__(self, seqs, k=15):
        self.k = k
        self.table = defaultdict(list)
        for sid, seq in seqs.items() if hasattr(seqs, "items") else enumerate(seqs):
            for i in range(len(seq) - k + 1):
                self.table[seq[i:i + k]].append((sid, i))

    def query(self, kmer):
        return self.table.get(kmer, [])

    def locate(self, seq, min_hits=1):
        """Map a query sequence to index hits: {seq_id: [positions]}."""
        hits = defaultdict(list)
        for i in range(len(seq) - self.k + 1):
            for sid, pos in self.table.get(seq[i:i + self.k], ()):
                hits[sid].append(pos)
        return {sid: ps for sid, ps in hits.items() if len(ps) >= min_hits}


def de_bruijn_contigs(reads, k=31, min_count=1):
    """Assemble reads into unitigs via a de Bruijn graph.

    Nodes are (k-1)-mers; each observed k-mer is an edge prefix->suffix.
    Unitigs are maximal non-branching paths. Returns sorted contigs
    (longest first). Requires error-free reads for exact reconstruction.
    """
    edges = defaultdict(set)
    indeg, outdeg = defaultdict(int), defaultdict(int)
    for read in reads:
        for i in range(len(read) - k + 1):
            km = read[i:i + k]
            u, v = km[:-1], km[1:]
            if v in edges[u]:
                continue  # deduplicate repeated k-mers
            edges[u].add(v)
            outdeg[u] += 1
            indeg[v] += 1
            indeg.setdefault(u, indeg.get(u, 0))
            outdeg.setdefault(v, outdeg.get(v, 0))
            edges.setdefault(v, edges.get(v, set()))

    contigs = []
    used = defaultdict(set)
    # start at branching/end nodes (indeg != 1 or outdeg != 1), then
    # cover any remaining cycles
    starts = [n for n in edges if indeg[n] != 1 or outdeg[n] != 1]
    starts += [n for n in edges if n not in starts]

    for node in starts:
        for nxt in sorted(edges[node] - used[node]):
            path = [node]
            cur, cur_in = node, nxt
            while True:
                path.append(cur_in)
                used[cur].add(cur_in)
                if indeg[cur_in] != 1 or outdeg[cur_in] != 1:
                    break
                nxt_edge = next(iter(edges[cur_in] - used[cur_in]), None)
                if nxt_edge is None:
                    break
                cur, cur_in = cur_in, nxt_edge
            contigs.append(path[0] + "".join(n[-1] for n in path[1:]))
    return sorted(contigs, key=len, reverse=True)
