import os
import re
import subprocess
import requests

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

base = os.path.join(project_root, "01_data")

def curl_dl(url, dest, use_proxy=True):
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        print(f"SKIP (exists): {dest} ({os.path.getsize(dest):,} bytes)", flush=True)
        return True
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    part = dest + ".part"
    cmd = ["curl", "-L", "-f", "-sS", "--connect-timeout", "15", "--retry", "3", "-o", part, url]
    if not use_proxy:
        cmd[1:1] = ["--noproxy", "*"]
    print(f"Downloading (proxy={use_proxy}): {os.path.basename(dest)}", flush=True)
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode == 0 and os.path.exists(part) and os.path.getsize(part) > 500:
        os.replace(part, dest)
        print(f"  DONE: {os.path.basename(dest)} ({os.path.getsize(dest):,} bytes)", flush=True)
        return True
    print(f"  FAIL rc={res.returncode}", flush=True)
    if use_proxy:
        return curl_dl(url, dest, use_proxy=False)
    return False

# 1. TCGA-BRCA HiSeqV2 (64 MB)
curl_dl("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FHiSeqV2.gz",
        os.path.join(base, "TCGA_GEO", "TCGA_BRCA_HiSeqV2.tsv.gz"))

# 2. GSE194040 (I-SPY2-990)
curl_dl("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194040/matrix/GSE194040-GPL30493_series_matrix.txt.gz",
        os.path.join(base, "TCGA_GEO", "GSE194040-GPL30493_series_matrix.txt.gz"))
curl_dl("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194040/suppl/GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt.gz",
        os.path.join(base, "TCGA_GEO", "GSE194040_geneLevel_n988.txt.gz"))

# 3. GSE58812, GSE103091, GSE173839, GSE96058 series matrices
s = requests.Session(); s.trust_env = True
for gse, prefix in [("GSE58812", "GSE58nnn"), ("GSE103091", "GSE103nnn"), ("GSE173839", "GSE173nnn"), ("GSE96058", "GSE96nnn")]:
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/matrix/"
    r = s.get(url, timeout=20)
    links = [x for x in re.findall(r'href="([^"]+)"', r.text) if x.endswith(".gz") and not x.startswith("?")]
    print(f"{gse}/matrix: {links}", flush=True)
    for fn in links:
        curl_dl(url + fn, os.path.join(base, "TCGA_GEO", fn))

# 4. Check GCST010100 harmonised folder for .h.tsv.gz (much smaller than raw 1.02GB .tsv!)
r_h = s.get("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010100/harmonised/", timeout=20)
h_files = [x for x in re.findall(r'href="([^"]+)"', r_h.text) if x.endswith(".gz") and not x.startswith("?")]
print("GCST010100 harmonised:", h_files, flush=True)
