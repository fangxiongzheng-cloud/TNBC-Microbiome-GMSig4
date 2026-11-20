import ast
import csv
import json
import os
import re
import requests

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

s = requests.Session()
s.trust_env = False

gut_dir = os.path.join(project_root, "01_data", "gutMGene")
gm_dir = os.path.join(project_root, "01_data", "GMrepo")
os.makedirs(gut_dir, exist_ok=True)
os.makedirs(gm_dir, exist_ok=True)

gut_headers = {
    "User-Agent": "Mozilla/5.0",
    "Referer": "http://bio-computing.hrbmu.edu.cn/gutmgene/",
    "Content-Type": "application/json",
}

def write_dicts_csv(path, rows):
    if not rows:
        return []
    cols = list(rows[0].keys())
    for r in rows:
        for k in r.keys():
            if k not in cols:
                cols.append(k)
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    return cols

# 1. Fetch all 3 Homo sapiens tables from gutMGene v2.0 /search/searchdetaildata
datasets = {
    "gutMGene2_human_microbe_metabolite.csv": "microbe_metabolite",
    "gutMGene2_human_microbe_gene.csv": "microbe_gene",
    "gutMGene2_human_metabolite_gene.csv": "metabolite_gene",
}

for fname, ds in datasets.items():
    payload = {
        "data": {
            "dataset": ds,
            "species": "human",
            "gutMicrobe": "",
            "gene": "",
            "metaboliteNames": "",
            "metaboliteChEBIOntology": "",
            "subMetanameslist": "",
            "datatype": "",
        }
    }
    r = s.post(
        "http://bio-computing.hrbmu.edu.cn/gutMGene2.0_api/search/searchdetaildata",
        json=payload,
        headers=gut_headers,
        timeout=60,
    )
    print(f"gutMGene {ds} -> HTTP {r.status_code}, len={len(r.text)}")
    if r.status_code == 200:
        js = r.json()
        rows = js.get("resultSet", [])
        out_csv = os.path.join(gut_dir, fname)
        cols = write_dicts_csv(out_csv, rows)
        print(f"  Saved {fname}: rows={len(rows)}, cols={cols}")
    else:
        print("  Error preview:", r.text[:300])

# 2. Parse GMrepo concatenated dict responses for D001943 and D006262 projects
def parse_concat_dicts(text):
    fixed = "[" + re.sub(r"\}\s*\{", "}, {", text.strip()) + "]"
    try:
        return ast.literal_eval(fixed)
    except Exception:
        return json.loads(fixed)

bc_proj_txt = open(os.path.join(gm_dir, "getAssociatedProjectsByMeshID_D001943.json"), encoding="utf-8").read()
bc_projs = parse_concat_dicts(bc_proj_txt)
write_dicts_csv(os.path.join(gm_dir, "D001943_projects.csv"), bc_projs)
print("\n=== GMrepo Breast Neoplasms (D001943) Projects ===")
bc_map = {p["project_id"]: p for p in bc_projs}
for p in bc_projs:
    print(f"  {p['project_id']}: runs={p.get('runs_per_project')}, loaded={p.get('loaded_samples')}, pubmed={p.get('Related_Publications')}")

r_h = s.post(
    "https://gmrepo.humangut.info/api/getAssociatedProjectsByMeshID/",
    data=json.dumps({"mesh_id": "D006262"}),
    headers={"User-Agent": "Mozilla/5.0"},
    timeout=30,
)
h_projs = parse_concat_dicts(r_h.text)
paired = [p for p in h_projs if p.get("project_id") in bc_map]
write_dicts_csv(os.path.join(gm_dir, "D001943_paired_health_projects.csv"), paired)
print("\n=== Breast Cancer Projects ALSO having Health (D006262) Controls in GMrepo ===")
for p in paired:
    bc = bc_map[p["project_id"]]
    print(f"  {p['project_id']}: BC_loaded={bc.get('loaded_samples')}, Health_loaded={p.get('loaded_samples')}, Health_runs={p.get('runs_per_project')}")

# Check GMrepo v2 frontend JS for bulk download paths
try:
    idx_js = s.get("https://gmrepo.humangut.info/static/index.fbde0e49.js", timeout=15).text
    dl_hits = sorted(set(re.findall(r'["\']([^"\']*(?:download|tsv|gz|RunsByProject|ProjectID)[^"\']*)["\']', idx_js, re.I)))
    print("\nGMrepo v2 download JS hints:", [x for x in dl_hits if len(x) < 120][:30])
except Exception as e:
    print("GMrepo JS check:", e)
