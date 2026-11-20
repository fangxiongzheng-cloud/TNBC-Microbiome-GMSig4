import csv
import json
import os
import requests
from importlib.machinery import SourceFileLoader

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

mod = SourceFileLoader("tbi_mod", os.path.join(os.path.dirname(__file__), "06k_tbi_range_extract.py")).load_module()

data_dir = os.path.join(project_root, "01_data")
res_dir  = os.path.join(project_root, "02_results")
out_csv  = os.path.join(data_dir, "eQTL", "GTEx_V8_Breast_target_genes_500kb.csv")

with open(os.path.join(res_dir, "Candidate_Genes_GRCh38_Coords.json"), "r", encoding="utf-8") as f:
    coords = json.load(f)

missing_genes = ["TCF7L2", "CCR6", "FOXP3", "CXCR6", "TGFBR1", "NLRP6", "ICAM1", "IFNG", "TGFB1", "LYZ", "DAPK3", "AKT1"]
ensg_set = {coords[g]["ensembl_id"] for g in missing_genes if g in coords}
queries = []
for g in missing_genes:
    if g in coords:
        info = coords[g]
        tss = info["start"] if info["strand"] == 1 else info["end"]
        queries.append((info["chr"], max(1, tss - 500000), tss + 500000))

url = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/Breast_Mammary_Tissue.tsv.gz"
tbi_path = os.path.join(data_dir, "eQTL", "GTEx_V8_Breast_Mammary_Tissue.tsv.gz.tbi")

s = requests.Session()
s.trust_env = True
h = s.head(url, timeout=30, allow_redirects=True, headers={"User-Agent": "Mozilla/5.0"})
total_size = int(h.headers["Content-Length"])

idx = mod.TabixIndex(tbi_path)
ranges = idx.byte_ranges(queries, merge_gap=131072, max_span=1048576)
final_ranges = [(a, min(total_size - 1, b)) for a, b in ranges]
print(f"Retrying {len(missing_genes)} missing GTEx Breast genes across {len(final_ranges)} 1MB ranges (threads=3)...", flush=True)

# Read existing rows & header
with open(out_csv, "r", encoding="utf-8") as f:
    reader = csv.reader(f)
    header = next(reader)
    existing_rows = list(reader)

seen = {"\t".join(r[:5]) for r in existing_rows}
g_idx = header.index("gene_id")

import concurrent.futures, time

def fetch_robust(rng):
    a, b = rng
    ls = requests.Session()
    ls.trust_env = True
    for att in range(10):
        try:
            r = ls.get(url, headers={"User-Agent": "Mozilla/5.0", "Range": f"bytes={a}-{b}"}, timeout=90)
            if r.status_code in (200, 206) and len(r.content) == (b - a + 1):
                out = []
                for ln in mod.parse_bgzf_lines(r.content):
                    if not ln or ln.startswith(header[0] + "\t"):
                        continue
                    parts = ln.split("\t")
                    if len(parts) > g_idx and parts[g_idx].split(".")[0] in ensg_set:
                        out.append(parts)
                return out
        except Exception:
            time.sleep(1.5 * (att + 1))
    print(f"  WARN still failed range {a}-{b}", flush=True)
    return []

added = 0
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
    for k, res_parts in enumerate(ex.map(fetch_robust, final_ranges), 1):
        for parts in res_parts:
            key = "\t".join(parts[:5])
            if key not in seen:
                seen.add(key)
                existing_rows.append(parts)
                added += 1
        if k % 5 == 0 or k == len(final_ranges):
            print(f"  [{k}/{len(final_ranges)}] +{added} rows (total={len(existing_rows)})", flush=True)

with open(out_csv, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(header)
    w.writerows(existing_rows)

print(f"DONE retry! Total GTEx Breast rows: {len(existing_rows)} (+{added})", flush=True)
