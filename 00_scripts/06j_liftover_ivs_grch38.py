import pandas as pd
from pyliftover import LiftOver
import os

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

res_dir = os.path.join(project_root, "02_results")
chain_fp = os.path.join(project_root, "01_data", "GWAS_MR", "hg19ToHg38.over.chain.gz")
lo = LiftOver(chain_fp)

mbg = pd.read_csv(f"{res_dir}/MiBioGen_clumped_IVs_1e5.csv", dtype={"chr": str})
eq  = pd.read_csv(f"{res_dir}/eQTLGen_candidate_instruments.csv", dtype={"chr.exposure": str})

pairs = {}
for _, r in mbg.iterrows():
    rs = str(r["rsID"]).strip()
    if rs.startswith("rs") and pd.notna(r["chr"]) and pd.notna(r["pos"]):
        pairs[rs] = (str(r["chr"]).replace("chr", ""), int(r["pos"]))

for _, r in eq.iterrows():
    rs = str(r["SNP"]).strip()
    if rs.startswith("rs") and pd.notna(r["chr.exposure"]) and pd.notna(r["pos.exposure"]):
        pairs[rs] = (str(r["chr.exposure"]).replace("chr", ""), int(r["pos.exposure"]))

rows = []
for rs, (ch, pos19) in sorted(pairs.items()):
    # UCSC coordinates are 0-based
    res = lo.convert_coordinate("chr" + ch, pos19 - 1)
    if res:
        c38 = str(res[0][0]).replace("chr", "")
        p38 = int(res[0][1]) + 1
        rows.append({"rsid": rs, "chr19": ch, "pos19": pos19, "chr38": c38, "pos38": p38})

df = pd.DataFrame(rows)
df.to_csv(f"{res_dir}/IV_rsids_GRCh38_map.csv", index=False)
print(f"Lifted {len(df)} / {len(pairs)} unique IV SNPs to GRCh38 locally in <1s!")
