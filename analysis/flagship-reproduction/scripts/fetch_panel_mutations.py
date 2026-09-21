#!/usr/bin/env python3
"""Fetch non-silent mutations for a curated BRCA driver panel.

Builds on fetch_tcga_brca.py outputs (all_samples.json, subtypes.tsv).
Writes data/panel_mutations.json for ml/ feature-matrix construction.
"""

import json
import time
import urllib.request
from pathlib import Path

BASE = "https://www.cbioportal.org/api"
STUDY = "brca_tcga_pan_can_atlas_2018"
PROFILE = f"{STUDY}_mutations"
SAMPLE_LIST = f"{STUDY}_all"
DATA = Path(__file__).resolve().parent.parent / "data"

# Curated TCGA/COSMIC breast-cancer drivers + clinically informative genes.
# Entrez IDs verified against NCBI gene records.
PANEL = {
    7157: "TP53", 5290: "PIK3CA", 2625: "GATA3", 207: "AKT1",
    999: "CDH1", 4217: "MAP3K1", 5872: "MAP2K4", 5728: "PTEN",
    58508: "KMT2C", 6926: "TBX3", 865: "CBFB", 861: "RUNX1",
    4763: "NF1", 23451: "SF3B1", 5295: "PIK3R1", 2064: "ERBB2",
    2065: "ERBB3", 3169: "FOXA1", 2099: "ESR1", 8289: "ARID1A",
    472: "ATM", 672: "BRCA1", 675: "BRCA2", 841: "CASP8",
    10664: "CTCF", 5925: "RB1", 6794: "STK11", 3845: "KRAS",
    4893: "NRAS", 673: "BRAF", 3265: "HRAS", 79728: "PALB2",
    11200: "CHEK2", 4089: "SMAD4", 595: "CCND1", 4609: "MYC",
    4851: "NOTCH1", 324: "APC", 7048: "TGFBR2", 3631: "INPP4B",
    10075: "HUWE1", 367: "AR", 7494: "XBP1", 8405: "SPOP",
    578: "PTPRD", 2321: "FLT3",
}

NON_SILENT_EXCLUDE = {
    "Silent", "Intron", "5'UTR", "3'UTR", "5'Flank", "3'Flank", "IGR",
    "Splice_Region",  # keep Splice_Site (usually consequential)
}


def get(path, params=None):
    url = BASE + path
    if params:
        url += "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read())


def post(path, body, params=None):
    url = BASE + path
    if params:
        url += "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(
        url, data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json",
                 "Accept": "application/json"}, method="POST")
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read())


def main():
    # gene_id -> {sample_id: n_nonsilent}
    out = {}
    entrez = sorted(PANEL)
    for i in range(0, len(entrez), 20):  # batches of 20 genes
        batch = entrez[i:i + 20]
        muts = post(
            f"/molecular-profiles/{PROFILE}/mutations/fetch",
            {"entrezGeneIds": batch, "sampleListId": SAMPLE_LIST},
            {"projection": "DETAILED"},
        )
        for m in muts:
            vt = m.get("mutationType", "")
            if vt in NON_SILENT_EXCLUDE:
                continue
            gid = m["entrezGeneId"]
            out.setdefault(gid, {}).setdefault(m["sampleId"], 0)
            out[gid][m["sampleId"]] += 1
        print(f"batch {i//20 + 1}: {len(batch)} genes, "
              f"cumulative {sum(len(v) for v in out.values())} hits")
        time.sleep(0.3)

    (DATA / "panel_mutations.json").write_text(
        json.dumps({"panel": PANEL,
                    "mutations": {str(k): v for k, v in out.items()}},
                   indent=1))
    print(f"wrote {DATA/'panel_mutations.json'}: "
          f"{len(out)} genes with hits")


if __name__ == "__main__":
    import urllib.parse
    main()
