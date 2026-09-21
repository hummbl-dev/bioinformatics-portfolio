#!/usr/bin/env python3
"""Extension arm: ERBB2/CCND1/MYC copy-number by PAM50 subtype.

The mutation-only ML arm showed Her2 is invisible to mutation
features (F1 0.18) because Her2 subtype is defined by ERBB2
*amplification*. This script fetches discrete CNA (GISTIC) calls and
tests amplification enrichment per subtype — the documented
deviation/extension of the reproduction design.
"""

import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

BASE = "https://www.cbioportal.org/api"
STUDY = "brca_tcga_pan_can_atlas_2018"
DATA = Path(__file__).resolve().parent.parent / "data"

GENES = {2064: "ERBB2", 595: "CCND1", 4609: "MYC", 5728: "PTEN",
         703: "FGFR1", 2261: "FGFR3"}
SUBTYPES = ["BRCA_Basal", "BRCA_Her2", "BRCA_LumA", "BRCA_LumB",
            "BRCA_Normal"]


def _open(req, timeout):
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read())
        except urllib.error.HTTPError as e:
            if e.code >= 500 and attempt < 4:
                time.sleep(2 ** (attempt + 1))
                continue
            raise


def get(path, params=None):
    url = BASE + path
    if params:
        url += "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    return _open(req, 60)


def post(path, body, params=None):
    url = BASE + path
    if params:
        url += "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(
        url, data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json",
                 "Accept": "application/json"}, method="POST")
    return _open(req, 120)


def fisher_two_sided(a, b, c, d):
    """Same implementation as subtype_mutation_association.py."""
    from math import comb, exp, lgamma

    def lchoose(n, k):
        if k < 0 or k > n:
            return float("-inf")
        return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)

    n = a + b + c + d
    r1 = a + b
    c1 = a + c

    def p_table(x):
        return exp(lchoose(c1, x) + lchoose(n - c1, r1 - x)
                   - lchoose(n, r1))

    lo = max(0, r1 - (n - c1))
    hi = min(r1, c1)
    p_obs = p_table(a)
    return min(1.0, sum(p_table(x) for x in range(lo, hi + 1)
                        if p_table(x) <= p_obs * (1 + 1e-9)))


def main():
    # Find the discrete CNA profile
    profiles = get(f"/studies/{STUDY}/molecular-profiles")
    cna = [p for p in profiles
           if "gistic" in p["molecularProfileId"].lower()
           or p.get("molecularAlterationType") == "COPY_NUMBER_ALTERATION"]
    for p in profiles:
        if p.get("molecularAlterationType") == "COPY_NUMBER_ALTERATION":
            cna = cna or [p]
    assert cna, "no CNA profile found"
    profile = cna[0]["molecularProfileId"]
    print(f"using CNA profile {profile}")

    lists = get(f"/studies/{STUDY}/sample-lists")
    cna_lists = [l for l in lists if "cna" in l["sampleListId"].lower()]
    sample_list = (cna_lists[0] if cna_lists else
                   [l for l in lists if l["sampleListId"].endswith("_all")][0]
                   )["sampleListId"]

    # patient -> subtype
    subs = {}
    for line in (DATA / "subtypes.tsv").read_text().splitlines()[1:]:
        pid, st = line.split("\t")
        subs[pid] = st

    calls = post(
        f"/molecular-profiles/{profile}/discrete-copy-number/fetch",
        {"entrezGeneIds": sorted(GENES), "sampleListId": sample_list},
        {"projection": "DETAILED",
         "discreteCopyNumberEventType": "HOMDEL_AND_AMP"})

    # patient -> gene -> max CNA level seen across samples
    per = {}
    for c in calls:
        pid = c["sampleId"][:12]
        per.setdefault(pid, {}).setdefault(
            c["entrezGeneId"], []).append(c["alteration"])
    amp = {pid: {g: (1 if max(v) >= 2 else 0) for g, v in genes.items()}
           for pid, genes in per.items()}
    dele = {pid: {g: (1 if min(v) <= -2 else 0) for g, v in genes.items()}
            for pid, genes in per.items()}

    results = []
    for entrez, sym in sorted(GENES.items(), key=lambda kv: kv[1]):
        for event, tab in (("AMP", amp), ("HOMDEL", dele)):
            n_mut = sum(1 for p in subs if tab.get(p, {}).get(entrez, 0))
            for st in SUBTYPES:
                cohort = [p for p, s in subs.items() if s == st]
                a = sum(1 for p in cohort if tab.get(p, {}).get(entrez, 0))
                b = len(cohort) - a
                c = n_mut - a
                d = len(subs) - len(cohort) - c
                pv = fisher_two_sided(a, b, c, d)
                or_ = (a * d / (b * c)) if b * c else (float("inf") if a * d else 0.0)
                results.append({"gene": sym, "event": event,
                                "subtype": st, "n": len(cohort),
                                "n_event": a,
                                "freq": round(a / len(cohort), 4),
                                "odds_ratio": round(or_, 3),
                                "fisher_p": pv})

    # BH across all tests
    results.sort(key=lambda r: r["fisher_p"])
    m = len(results)
    for i, r in enumerate(results):
        r["fisher_q"] = round(min(1.0, r["fisher_p"] * m / (i + 1)), 6)
        r["fisher_p"] = round(r["fisher_p"], 6)

    (DATA / "cna_association_results.json").write_text(
        json.dumps({"profile": profile, "n_subtyped": len(subs),
                    "genes": GENES, "results": results}, indent=1))

    for r in results:
        if r["fisher_q"] < 0.05:
            print(f"{r['gene']:6s} {r['event']:6s} {r['subtype']:12s} "
                  f"{r['n_event']:>3}/{r['n']} ({r['freq']*100:5.1f}%) "
                  f"OR={r['odds_ratio']:7.2f} q={r['fisher_q']:.2e}")
    print(f"wrote {DATA/'cna_association_results.json'}")


if __name__ == "__main__":
    main()
