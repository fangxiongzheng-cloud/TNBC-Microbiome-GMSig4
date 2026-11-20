import concurrent.futures
import csv
import os
import requests
from threading import local

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

out_csv = os.path.join(project_root, "01_data", "GWAS_MR", "MiBioGen_tophits_all.csv")
os.makedirs(os.path.dirname(out_csv), exist_ok=True)

tl = local()
def get_sess():
    if not hasattr(tl, "s"):
        s = requests.Session()
        s.trust_env = False
        tl.s = s
    return tl.s

r0 = get_sess().get("https://mibiogen.gcc.rug.nl/api/v2/mibiogen_tophits?start=0&num=10000", timeout=30).json()
total = r0["total"]
print(f"Total MiBioGen tophits: {total}", flush=True)

starts = list(range(10000, total, 10000))

def fetch_page(st):
    s = get_sess()
    for _ in range(3):
        try:
            r = s.get(f"https://mibiogen.gcc.rug.nl/api/v2/mibiogen_tophits?start={st}&num=10000", timeout=60)
            if r.status_code == 200:
                return st, r.json().get("items", [])
        except Exception:
            pass
    return st, []

pages = {0: r0["items"]}
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
    for st, items in ex.map(fetch_page, starts):
        pages[st] = items
        print(f"  Fetched page start={st}: {len(items)} rows", flush=True)

def extract(it):
    ref = it.get("ref")
    alt = it.get("alt")
    bac = it.get("bac")
    return {
        "bac": bac.get("taxon", "") if isinstance(bac, dict) else (bac or ""),
        "chr": it.get("chr"),
        "pos": it.get("pos"),
        "rsID": it.get("rsID"),
        "ref": ref.get("id", "") if isinstance(ref, dict) else (ref or ""),
        "alt": alt.get("id", "") if isinstance(alt, dict) else (alt or ""),
        "beta": it.get("beta"),
        "SE": it.get("SE"),
        "metaZ_wz": it.get("metaZ_wz"),
        "P": it.get("P"),
        "N": it.get("N"),
        "Ncohorts": it.get("Ncohorts"),
    }

all_rows = []
for st in sorted(pages):
    for it in pages[st]:
        all_rows.append(extract(it))

cols = ["bac", "chr", "pos", "rsID", "ref", "alt", "beta", "SE", "metaZ_wz", "P", "N", "Ncohorts"]
with open(out_csv, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    w.writerows(all_rows)
print(f"SAVED {out_csv}: {len(all_rows)} rows", flush=True)
