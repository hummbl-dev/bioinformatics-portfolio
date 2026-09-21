"""Reproduce: subtype-enriched driver mutations in TCGA-BRCA.

Reference claim (TCGA Network 2012, Nature 490:61-70): TP53 mutations
are enriched in basal-like tumours (~80%+) vs luminal A (~12-29%);
PIK3CA enriched in luminal A (~45%) vs basal (~7-9%).

Method: per-gene 2x2 contingency (subtype vs rest) -> Fisher exact
(two-sided, hypergeometric), odds ratio, Benjamini-Hochberg across all
tests. Pure stdlib.

Run: python analysis/flagship-reproduction/scripts/subtype_mutation_association.py
"""

import json
import math
from collections import defaultdict
from pathlib import Path

DATA = Path(__file__).resolve().parents[1] / "data"
GENES = ["tp53", "pik3ca", "gata3"]
SUBTYPES = ["BRCA_Basal", "BRCA_Her2", "BRCA_LumA", "BRCA_LumB",
            "BRCA_Normal"]


def fisher_two_sided(a, b, c, d):
    """Exact two-sided Fisher p-value: sum of hypergeometric
    probabilities <= P(observed table), fixing margins."""
    n = a + b + c + d
    r1, c1 = a + b, a + c
    lo = max(0, r1 - (n - c1))
    hi = min(r1, c1)

    def prob(x):
        return (math.comb(c1, x) * math.comb(n - c1, r1 - x)
                / math.comb(n, r1))

    p_obs = prob(a)
    return min(1.0, sum(prob(x) for x in range(lo, hi + 1)
                        if prob(x) <= p_obs * (1 + 1e-9)))


def odds_ratio(a, b, c, d):
    if b * c == 0:
        return float("inf") if a * d > 0 else float("nan")
    return (a * d) / (b * c)


def main():
    # patientId -> pam50 (subtype attr is patient-level)
    lines = (DATA / "subtypes.tsv").read_text().strip().splitlines()
    key = lines[0].split("\t")[0]
    pam = {r.split("\t")[0]: r.split("\t")[1] for r in lines[1:]}
    assert key in ("patientId", "sampleId"), key

    # gene -> set of mutated patients (non-silent mutations only)
    mut_by_gene = {}
    for g in GENES:
        muts = json.loads((DATA / f"mutations_{g}.json").read_text())
        ns = {m["patientId"] for m in muts
              if m.get("mutationType", "") not in
              ("Silent", "Intron", "IGR", "3'UTR", "5'UTR", "5'Flank",
               "3'Flank")}
        mut_by_gene[g.upper()] = ns

    # cohort = patients with a subtype call
    cohort = {p: s for p, s in pam.items() if s in SUBTYPES}
    n_cohort = len(cohort)
    by_subtype = defaultdict(set)
    for p, s in cohort.items():
        by_subtype[s].add(p)

    rows = []
    for gene, muts in mut_by_gene.items():
        overall = len(muts & set(cohort)) / n_cohort
        print(f"\n{gene}: {len(muts & set(cohort))}/{n_cohort} "
              f"mutated ({overall:.1%})")
        for st in SUBTYPES:
            sub = by_subtype[st]
            a = len(muts & sub)                    # mutated & subtype
            b = len(sub) - a                       # wt & subtype
            c = len(muts & set(cohort) - sub)      # mutated & other
            d = len(set(cohort) - sub - muts)      # wt & other
            p = fisher_two_sided(a, b, c, d)
            orr = odds_ratio(a, b, c, d)
            freq = a / len(sub)
            rows.append((gene, st, a, len(sub), freq, orr, p))
            print(f"  {st:7s} {a:4d}/{len(sub):4d} ({freq:5.1%}) "
                  f"OR={orr:6.2f} p={p:.2e}")

    # Benjamini-Hochberg
    order = sorted(range(len(rows)), key=lambda i: rows[i][6])
    m = len(rows)
    out = []
    for rank, i in enumerate(order, 1):
        g, st, a, n_sub, freq, orr, p = rows[i]
        out.append({"gene": g, "subtype": st, "mutated": a,
                    "n": n_sub, "freq": round(freq, 4),
                    "odds_ratio": round(orr, 3), "p": p,
                    "p_adj_bh": round(min(1.0, p * m / rank), 6)})
    (DATA / "association_results.json").write_text(
        json.dumps(out, indent=1))
    print(f"\nwrote {DATA / 'association_results.json'}")


if __name__ == "__main__":
    main()
