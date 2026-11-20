import os
import re
import requests

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

s = requests.Session()
s.trust_env = True  # Use the environment's proxy settings; fallback handled per host

def download_file(url, dest_path, allow_direct_fallback=True):
    if os.path.exists(dest_path) and os.path.getsize(dest_path) > 1000:
        print(f"SKIP (exists): {dest_path} ({os.path.getsize(dest_path)} bytes)")
        return True
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    for use_proxy in ([True, False] if allow_direct_fallback else [True]):
        sess = requests.Session()
        sess.trust_env = use_proxy
        try:
            with sess.get(url, stream=True, timeout=60, headers={"User-Agent": "Mozilla/5.0"}) as r:
                if r.status_code != 200:
                    print(f"HTTP {r.status_code} (proxy={use_proxy}) for {url}")
                    continue
                tmp = dest_path + ".part"
                total = 0
                with open(tmp, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1024 * 256):
                        if chunk:
                            f.write(chunk)
                            total += len(chunk)
                os.replace(tmp, dest_path)
                print(f"SAVED (proxy={use_proxy}): {dest_path} ({total:,} bytes)")
                return True
        except Exception as e:
            print(f"ERR (proxy={use_proxy}) for {url}: {e}")
    return False

base = os.path.join(project_root, "01_data")

# 1. MiBioGen allHits
download_file(
    "https://molgenis26.gcc.rug.nl/downloads/MiBioGen/MBG.allHits.p1e4.txt",
    os.path.join(base, "GWAS_MR", "MiBioGen_MBG.allHits.p1e4.txt"),
)

# 2. TCGA-BRCA (Xena)
download_file(
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix",
    os.path.join(base, "TCGA_GEO", "TCGA_BRCA_clinicalMatrix.tsv"),
)
download_file(
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival%2FBRCA_survival.txt",
    os.path.join(base, "TCGA_GEO", "TCGA_BRCA_survival.tsv"),
)
download_file(
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FHiSeqV2.gz",
    os.path.join(base, "TCGA_GEO", "TCGA_BRCA_HiSeqV2.tsv.gz"),
)

# 3. GSE194040 (I-SPY2-990 gene-level + metadata)
download_file(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194040/matrix/GSE194040-GPL30493_series_matrix.txt.gz",
    os.path.join(base, "TCGA_GEO", "GSE194040-GPL30493_series_matrix.txt.gz"),
)
download_file(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194040/suppl/GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt.gz",
    os.path.join(base, "TCGA_GEO", "GSE194040_geneLevel_n988.txt.gz"),
)

# 4. Inspect GSE58812, GSE103091, GSE173839 directory listings & download
for gse, prefix in [("GSE58812", "GSE58nnn"), ("GSE103091", "GSE103nnn"), ("GSE173839", "GSE173nnn")]:
    for sub in ["matrix", "suppl"]:
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/{sub}/"
        r = requests.get(url, timeout=20)
        links = [x for x in re.findall(r'href="([^"]+)"', r.text) if x.endswith(".gz") or x.endswith(".txt") or x.endswith(".tar")]
        print(f"{gse}/{sub}: {links}")
        for fn in links:
            if sub == "matrix" or ("series_matrix" in fn or gse in fn):
                if fn.endswith(".tar"):
                    continue
                download_file(url + fn, os.path.join(base, "TCGA_GEO", fn))
