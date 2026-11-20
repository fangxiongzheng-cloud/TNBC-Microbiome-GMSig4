import csv
import json
import os
import sys

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

sys.path.insert(0, os.path.dirname(__file__))
from importlib.machinery import SourceFileLoader
mod = SourceFileLoader("tbi_mod", os.path.join(os.path.dirname(__file__), "06k_tbi_range_extract.py")).load_module()

data_dir = os.path.join(project_root, "01_data")
res_dir  = os.path.join(project_root, "02_results")

with open(os.path.join(res_dir, "Candidate_Genes_GRCh38_Coords.json"), "r", encoding="utf-8") as f:
    coords = json.load(f)

target_genes = [
    "PEMT", "TCF7L2", "CCR6", "MAGI3", "FOXP3", "CXCR6", "ST6GAL1", "LDHA",
    "TGFBR1", "NLRP6", "NLRP3", "MUC1", "CASP1", "ICAM1", "IFNG", "IDO1",
    "CXCL10", "ENO1", "TGFB1", "LYZ", "DAPK3", "AKT1", "VEGFA", "CLDN7", "F11R",
    "FFAR2", "CNKSR3", "SFN"
]

import requests
missing = [g for g in target_genes if g not in coords]
if missing:
    try:
        r = requests.post("https://mygene.info/v3/query", data={"q": ",".join(missing), "scopes": "symbol", "species": "human", "fields": "ensembl.gene,genomic_pos"}, timeout=10)
        if r.status_code == 200:
            for item in r.json():
                q = item.get("query")
                gp = item.get("genomic_pos")
                if isinstance(gp, list): gp = gp[0]
                if isinstance(gp, dict) and gp.get("chr") and gp.get("ensemblgene"):
                    coords[q] = {"ensembl_id": gp["ensemblgene"], "chr": str(gp["chr"]), "start": int(gp["start"]), "end": int(gp["end"]), "strand": int(gp.get("strand", 1))}
            print(f"Resolved missing genes via mygene.info: {[g for g in missing if g in coords]}", flush=True)
    except Exception as e:
        print(f"mygene fallback note: {e}", flush=True)

ensg_set = {coords[g]["ensembl_id"] for g in target_genes if g in coords}
regions = []
for g in target_genes:
    if g not in coords:
        continue
    info = coords[g]
    tss = info["start"] if info["strand"] == 1 else info["end"]
    beg = max(1, tss - 500000)
    end = tss + 500000
    regions.append((g, info["ensembl_id"], info["chr"], beg, end))

queries = [(r[2], r[3], r[4]) for r in regions]
print(f"Prepared {len(regions)} gene locus (+/-500kb TSS) queries across GTEx Breast + GCST010100 TNBC.", flush=True)

def filter_gtex(parts, hdr):
    if len(parts) < 12:
        return False
    g_idx = hdr.index("gene_id") if "gene_id" in hdr else 6
    return parts[g_idx].split(".")[0] in ensg_set

# Build fast interval lookup per chrom for GCST010100 region extraction
reg_by_chr = {}
for g, ensg, ch, beg, end in regions:
    reg_by_chr.setdefault(ch.replace("chr", ""), []).append((beg, end))

def filter_tnbc_region(parts, hdr):
    if len(parts) < 10:
        return False
    ch = parts[0].replace("chr", "")
    if ch not in reg_by_chr:
        return False
    try:
        pos = int(parts[1])
        for b, e in reg_by_chr[ch]:
            if b <= pos <= e:
                return True
    except Exception:
        pass
    return False

# 1. Extract GTEx V8 Breast_Mammary_Tissue cis-eQTLs (+/- 500kb) for target genes
mod.remote_tabix_extract(
    "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/Breast_Mammary_Tissue.tsv.gz",
    os.path.join(data_dir, "eQTL", "GTEx_V8_Breast_Mammary_Tissue.tsv.gz.tbi"),
    queries, filter_gtex,
    os.path.join(data_dir, "eQTL", "GTEx_V8_Breast_target_genes_500kb.csv"),
    threads=16
)

# 2. Extract GCST010100 TNBC GWAS dense locus regions (+/- 500kb) for true Breast-eQTL coloc.abf
# Exclude chrX if GCST010100 doesn't index X/23; TabixIndex skips missing refs cleanly
mod.remote_tabix_extract(
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010100/harmonised/GCST010100.h.tsv.gz",
    os.path.join(data_dir, "GWAS_MR", "GCST010100.h.tsv.gz.tbi"),
    queries, filter_tnbc_region,
    os.path.join(data_dir, "GWAS_MR", "GCST010100_coloc_regions_500kb.csv"),
    threads=16
)
