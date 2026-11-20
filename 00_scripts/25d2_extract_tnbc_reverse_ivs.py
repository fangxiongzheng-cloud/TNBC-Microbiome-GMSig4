import os
import csv
import time
import requests
import concurrent.futures
import pandas as pd
from importlib.machinery import SourceFileLoader
from pyliftover import LiftOver

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

mod = SourceFileLoader("tbi_mod", os.path.join(project_root, "00_scripts", "06k_tbi_range_extract.py")).load_module()

GWAS_DIR = os.path.join(project_root, "01_data", "GWAS_MR")
RES_DIR  = os.path.join(project_root, "02_results")

bc_ivs = pd.read_csv(os.path.join(RES_DIR, "ReverseMR_BC_Instruments_OpenGWAS.csv"))
lo = LiftOver(os.path.join(GWAS_DIR, "hg19ToHg38.over.chain.gz"))

# Prioritize ER- & Overall BC loci lifted to GRCh38
queries = []
target_snps = {}
for _, r in bc_ivs.drop_duplicates("SNP").iterrows():
    ch = str(r["chr.exposure"]).replace("chr", "")
    pos = int(r["pos.exposure"])
    res = lo.convert_coordinate(f"chr{ch}", pos - 1)
    if res:
        c38 = res[0][0].replace("chr", "")
        p38 = int(res[0][1]) + 1
        queries.append((c38, max(1, p38 - 500), p38 + 500))
        target_snps[r["SNP"]] = (c38, p38)

tnbc_url = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010100/harmonised/GCST010100.h.tsv.gz"
tnbc_tbi = os.path.join(GWAS_DIR, "GCST010100.h.tsv.gz.tbi")

coloc_tnbc = pd.read_csv(os.path.join(GWAS_DIR, "GCST010100_coloc_regions_500kb.csv"))
header = list(coloc_tnbc.columns)
chr_i = header.index("chromosome")
pos_i = header.index("base_pair_location")
rs_i  = header.index("rsid")
p_i   = header.index("p_value")

s = requests.Session()
s.trust_env = True
h = s.head(tnbc_url, timeout=30, allow_redirects=True, headers={"User-Agent": "Mozilla/5.0"})
total_size = int(h.headers["Content-Length"])

idx = mod.TabixIndex(tnbc_tbi)
ranges = idx.byte_ranges(queries, merge_gap=65536, max_span=524288)
final_ranges = [(a, min(total_size - 1, b)) for a, b in ranges]
print(f"Querying {len(target_snps)} lifted BC loci across {len(final_ranges)} compressed ranges in GCST010100...", flush=True)

rs_set = set(target_snps.keys())
pos_map = {(c, p) for c, p in target_snps.values()}

def fetch_rng(rng):
    a, b = rng
    ls = requests.Session()
    ls.trust_env = True
    for att in range(6):
        try:
            r = ls.get(tnbc_url, headers={"User-Agent": "Mozilla/5.0", "Range": f"bytes={a}-{b}"}, timeout=60)
            if r.status_code in (200, 206):
                out = []
                for ln in mod.parse_bgzf_lines(r.content):
                    if not ln or ln.startswith(header[0] + "\t"):
                        continue
                    parts = ln.split("\t")
                    if len(parts) <= max(chr_i, pos_i, rs_i, p_i):
                        continue
                    c = parts[chr_i].replace("chr", "")
                    try:
                        pv = float(parts[p_i])
                        ps = int(parts[pos_i])
                    except Exception:
                        continue
                    if parts[rs_i] in rs_set or (c, ps) in pos_map or pv < 5e-8:
                        out.append(parts)
                return out
        except Exception:
            time.sleep(1.2 * (att + 1))
    return []

rows_ex = []
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
    for k, res_p in enumerate(ex.map(fetch_rng, final_ranges), 1):
        rows_ex.extend(res_p)
        if k % 25 == 0 or k == len(final_ranges):
            print(f"  [{k}/{len(final_ranges)}] collected {len(rows_ex)} rows", flush=True)

ex_df = pd.DataFrame(rows_ex, columns=header)
for col in ["base_pair_location", "effect_allele_frequency", "beta", "standard_error", "p_value"]:
    if col in ex_df.columns:
        ex_df[col] = pd.to_numeric(ex_df[col], errors="coerce")

tn_coloc_gws = coloc_tnbc[coloc_tnbc["p_value"] < 5e-8].sort_values("p_value").drop_duplicates("chromosome")
comb = pd.concat([ex_df[ex_df["p_value"] < 1e-5], tn_coloc_gws], ignore_index=True)
comb = comb.sort_values("p_value").drop_duplicates("rsid")
comb.to_csv(os.path.join(RES_DIR, "ReverseMR_TNBC_GCST010100_CandidateIVs.csv"), index=False)
print(f"SAVED {len(comb)} TNBC GCST010100 candidate IVs (GWS p<5e-8: {(comb['p_value']<5e-8).sum()})", flush=True)
