import json, os, glob
from collections import Counter

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

hgmt_dir = os.path.join(project_root, "01_data", "HGMT")
print("=== HGMT 8 Breast Cancer Cohorts Inventory ===")
projs = json.load(open(os.path.join(hgmt_dir, "projects_D001943.json"), encoding="utf-8"))
for p in projs:
    pid = p["project_id"]
    det = json.load(open(os.path.join(hgmt_dir, f"{pid}_details.json"), encoding="utf-8"))
    runs = json.load(open(os.path.join(hgmt_dir, f"{pid}_runs.json"), encoding="utf-8"))
    phenos = Counter(r.get("phenotype_name") for r in runs)
    has_ko = os.path.exists(os.path.join(hgmt_dir, f"{pid}_feature_ko.json"))
    print(f"  {pid} | {p.get('assay_type')} | summary={p.get('summary')} | runs={len(runs)} | phenos={dict(phenos)} | has_KO={has_ko} | PMID={det.get('public_accession')}")
    if runs:
        # Check extra metadata fields
        sample_desc = [r.get("sample_description") for r in runs if r.get("sample_description")]
        if sample_desc:
            print(f"    sample_desc examples: {sample_desc[:4]}")

# Check GCST90454344 files
import requests, re
s = requests.Session(); s.trust_env = True
r = s.get("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90454001-GCST90455000/GCST90454344/", timeout=20)
files = [x for x in re.findall(r'href="([^"]+)"', r.text) if not x.startswith('?') and not x.startswith('/')]
print("\nGCST90454344 files:", files)
for fn in files:
    if fn.endswith('.tsv.gz') or fn.endswith('.tsv') or 'yaml' in fn:
        u = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90454001-GCST90455000/GCST90454344/" + fn
        if 'yaml' in fn:
            print(fn, ":\n", s.get(u, timeout=15).text[:800])
        else:
            h = s.head(u, timeout=15)
            print(fn, "size_MB=", round(int(h.headers.get("Content-Length", 0))/1048576, 1))
