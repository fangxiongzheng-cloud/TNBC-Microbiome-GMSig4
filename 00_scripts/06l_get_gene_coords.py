import json
import os
import requests

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

genes_all = [
    "CASP1", "IFNG", "DAPK3", "TCF7L2", "AKT1", "ENO1", "TGFBR1", "FOXP3",
    "TGFB1", "LDHA", "IDO1", "VEGFA", "F11R", "CLDN7", "CXCR2", "CTSD",
    "HERC5", "ST6GAL1", "CXCR6", "MAGI3", "ICAM1", "CXCL10", "CCR6", "LYZ",
    "PEMT", "SFN", "NLRP6", "NLRP3", "MUC1", "IL12B", "MAPK3", "HDAC11", "BCL2",
    "FFAR2", "CNKSR3", "RAC1", "MAPK8"
]

out_path = os.path.join(project_root, "02_results", "Candidate_Genes_GRCh38_Coords.json")

for use_proxy in (True, False):
    try:
        s = requests.Session()
        s.trust_env = use_proxy
        r = s.post(
            "https://rest.ensembl.org/lookup/symbol/homo_sapiens",
            json={"symbols": genes_all},
            headers={"Content-Type": "application/json"},
            timeout=25,
        )
        r.raise_for_status()
        d = r.json()
        clean = {
            k: {
                "ensembl_id": v.get("id"),
                "chr": str(v.get("seq_region_name")),
                "start": int(v.get("start")),
                "end": int(v.get("end")),
                "strand": int(v.get("strand", 1)),
            }
            for k, v in d.items()
            if isinstance(v, dict) and v.get("id")
        }
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(clean, f, indent=2)
        print(f"Saved {len(clean)}/{len(genes_all)} gene coords (proxy={use_proxy}) to {out_path}")
        break
    except Exception as e:
        print(f"Attempt proxy={use_proxy} failed: {e}")
