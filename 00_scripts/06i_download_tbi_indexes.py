import os
import requests

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

s = requests.Session()
s.trust_env = True
base = os.path.join(project_root, "01_data")

targets = [
    ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010100/harmonised/GCST010100.h.tsv.gz.tbi",
     os.path.join(base, "GWAS_MR", "GCST010100.h.tsv.gz.tbi")),
    ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90454001-GCST90455000/GCST90454344/harmonised/GCST90454344.h.tsv.gz.tbi",
     os.path.join(base, "GWAS_MR", "GCST90454344.h.tsv.gz.tbi")),
    ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010098/harmonised/GCST010098.h.tsv.gz.tbi",
     os.path.join(base, "GWAS_MR", "GCST010098.h.tsv.gz.tbi")),
    ("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/Breast_Mammary_Tissue.tsv.gz.tbi",
     os.path.join(base, "eQTL", "GTEx_V8_Breast_Mammary_Tissue.tsv.gz.tbi")),
]

for url, dest in targets:
    if os.path.exists(dest) and os.path.getsize(dest) > 100000:
        print(f"SKIP {os.path.basename(dest)}: {os.path.getsize(dest):,}")
        continue
    r = s.get(url, timeout=60)
    r.raise_for_status()
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    with open(dest, "wb") as f:
        f.write(r.content)
    print(f"SAVED {os.path.basename(dest)}: {len(r.content):,} bytes")
