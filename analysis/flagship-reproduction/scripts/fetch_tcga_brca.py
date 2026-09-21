"""Fetch TCGA-BRCA (PanCancer Atlas 2018) mutation + subtype data from
the cBioPortal public API. No auth required.

Outputs to analysis/flagship-reproduction/data/:
  mutations_<gene>.json  - raw mutation records per gene
  subtypes.tsv           - sampleId -> PAM50 subtype
  samples.json           - sample list

Run: python analysis/flagship-reproduction/scripts/fetch_tcga_brca.py
"""

import json
import sys
import time
import urllib.request
from pathlib import Path

BASE = "https://www.cbioportal.org/api"
STUDY = "brca_tcga_pan_can_atlas_2018"
PROFILE = f"{STUDY}_mutations"
SAMPLE_LIST = f"{STUDY}_all"
GENES = {"TP53": 7157, "PIK3CA": 5290, "GATA3": 2625}

OUT = Path(__file__).resolve().parents[1] / "data"


def get(path, timeout=120):
    req = urllib.request.Request(
        BASE + path,
        headers={"Accept": "application/json"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.load(r)
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    samples = get(f"/sample-lists/{SAMPLE_LIST}")["sampleIds"]
    (OUT / "samples.json").write_text(json.dumps(samples, indent=1))
    print(f"samples: {len(samples)}")

    # PAM50 subtype attribute id varies; find it
    attrs = get(f"/studies/{STUDY}/clinical-attributes"
                f"?clinicalDataType=SAMPLE&projection=SUMMARY&pageSize=500")
    pam = [a for a in attrs
           if "PAM50" in a["clinicalAttributeId"].upper()
           or "SUBTYPE" in a["clinicalAttributeId"].upper()]
    print("candidate subtype attrs:",
          [a["clinicalAttributeId"] for a in pam])
    attr = pam[0]["clinicalAttributeId"]

    clin = []
    for lvl, idkey in (("SAMPLE", "sampleId"), ("PATIENT", "patientId")):
        clin = get(f"/studies/{STUDY}/clinical-data"
                   f"?clinicalDataType={lvl}&attributeId={attr}"
                   f"&projection=SUMMARY&pageSize=2000")
        if clin:
            with open(OUT / "subtypes.tsv", "w") as f:
                f.write(f"{idkey}\tpam50\n")
                for row in clin:
                    f.write(f"{row[idkey]}\t{row['value']}\n")
            print(f"subtyped {lvl.lower()}s: {len(clin)} via {attr}")
            break
    else:
        print("WARNING: no subtype rows at either level")

    for gene, eg in GENES.items():
        muts = get(f"/molecular-profiles/{PROFILE}/mutations"
                   f"?sampleListId={SAMPLE_LIST}&entrezGeneId={eg}"
                   f"&projection=DETAILED")
        (OUT / f"mutations_{gene.lower()}.json").write_text(
            json.dumps(muts, indent=1))
        uniq = {m["sampleId"] for m in muts}
        print(f"{gene}: {len(muts)} mutation records across "
              f"{len(uniq)} samples")


if __name__ == "__main__":
    sys.exit(main())
