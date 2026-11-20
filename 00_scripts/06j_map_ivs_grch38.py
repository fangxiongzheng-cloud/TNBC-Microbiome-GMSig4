import csv
import os
import requests

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

res_dir = os.path.join(project_root, "02_results")

def read_col(path, colname):
    vals = set()
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        r = csv.DictReader(f)
        for row in r:
            v = (row.get(colname) or "").strip()
            if v.startswith("rs"):
                vals.add(v)
    return vals

rsids = sorted(
    read_col(os.path.join(res_dir, "MiBioGen_clumped_IVs_1e5.csv"), "rsID") |
    read_col(os.path.join(res_dir, "eQTLGen_candidate_instruments.csv"), "SNP")
)
print(f"Total unique IV rsIDs to map to GRCh38: {len(rsids)}", flush=True)

map_cache = os.path.join(res_dir, "IV_rsids_GRCh38_map.csv")
existing = {}
if os.path.exists(map_cache) and os.path.getsize(map_cache) > 50:
    with open(map_cache, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            existing[row["rsid"]] = (row["chr38"], int(row["pos38"]))

todo = [r for r in rsids if r not in existing]
print(f"Mapping {len(todo)} rsIDs via MyVariant hg38 / Ensembl REST...", flush=True)

s = requests.Session()
s.trust_env = True
valid_chrs = {str(x) for x in range(1, 23)} | {"X", "Y"}

for i in range(0, len(todo), 250):
    batch = todo[i:i+250]
    try:
        resp = s.post(
            "https://myvariant.info/v1/query",
            params={"scopes": "dbsnp.rsid", "fields": "chrom,hg38.start,dbsnp.chrom,dbsnp.hg38", "assembly": "hg38"},
            json={"ids": batch},
            timeout=30
        )
        if resp.status_code == 200:
            for item in resp.json():
                q = item.get("query")
                chrom = item.get("chrom")
                pos = None
                hg38 = item.get("hg38")
                if isinstance(hg38, list) and hg38:
                    hg38 = hg38[0]
                if isinstance(hg38, dict):
                    pos = hg38.get("start")
                if (not chrom or not pos) and isinstance(item.get("dbsnp"), dict):
                    chrom = chrom or item["dbsnp"].get("chrom")
                    d38 = item["dbsnp"].get("hg38")
                    if isinstance(d38, list) and d38:
                        d38 = d38[0]
                    if isinstance(d38, dict):
                        pos = pos or d38.get("start")
                if chrom and pos:
                    c_clean = str(chrom).replace("chr", "")
                    if c_clean in valid_chrs:
                        existing[q] = (c_clean, int(pos))
    except Exception as e:
        print("MyVariant error:", e, flush=True)
    if (i // 250 + 1) % 3 == 0 or i + 250 >= len(todo):
        print(f"  MyVariant {min(i+250, len(todo))}/{len(todo)} -> mapped {len(existing)}", flush=True)

rem = [r for r in rsids if r not in existing]
print(f"Ensembl REST fallback for {len(rem)} remaining rsIDs...", flush=True)
for i in range(0, len(rem), 150):
    batch = rem[i:i+150]
    try:
        resp = s.post(
            "https://rest.ensembl.org/variation/homo_sapiens",
            headers={"Content-Type": "application/json", "Accept": "application/json"},
            json={"ids": batch},
            timeout=45
        )
        if resp.status_code == 200:
            for rsid, info in resp.json().items():
                if not isinstance(info, dict):
                    continue
                for m in info.get("mappings", []):
                    seq = str(m.get("seq_region_name", ""))
                    if m.get("coord_system") == "chromosome" and seq in valid_chrs:
                        existing[rsid] = (seq, int(m["start"]))
                        break
    except Exception as e:
        print("Ensembl error:", e, flush=True)

with open(map_cache, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["rsid", "chr38", "pos38"])
    for rsid in sorted(existing):
        c, p = existing[rsid]
        w.writerow([rsid, c, p])

print(f"DONE GRCh38 mapping: {len(existing)} / {len(rsids)}", flush=True)
