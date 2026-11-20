import concurrent.futures
import os
import re
import requests

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

def range_download(url, dest, threads=8, use_proxy=True):
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        print(f"SKIP (exists): {os.path.basename(dest)} ({os.path.getsize(dest):,} bytes)", flush=True)
        return True
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    s = requests.Session()
    s.trust_env = use_proxy
    h = s.head(url, timeout=20, allow_redirects=True, headers={"User-Agent": "Mozilla/5.0"})
    if h.status_code != 200:
        print(f"HEAD failed {h.status_code} for {url}", flush=True)
        return False
    size = int(h.headers.get("Content-Length", 0))
    accept_ranges = h.headers.get("Accept-Ranges", "").lower() == "bytes"
    part = dest + ".part"

    if not accept_ranges or size < 2 * 1024 * 1024:
        print(f"Single-stream downloading {os.path.basename(dest)} ({size:,} bytes)...", flush=True)
        with s.get(url, stream=True, timeout=60, headers={"User-Agent": "Mozilla/5.0"}) as r:
            r.raise_for_status()
            with open(part, "wb") as f:
                for ch in r.iter_content(256 * 1024):
                    if ch:
                        f.write(ch)
        os.replace(part, dest)
        print(f"  DONE: {os.path.basename(dest)} ({os.path.getsize(dest):,} bytes)", flush=True)
        return True

    print(f"Multi-range ({threads}x) downloading {os.path.basename(dest)} ({size/1048576:.1f} MB)...", flush=True)
    with open(part, "wb") as f:
        f.truncate(size)

    chunk_size = (size + threads - 1) // threads
    ranges = [(i * chunk_size, min(size - 1, (i + 1) * chunk_size - 1)) for i in range(threads)]

    def dl_chunk(idx_rng):
        idx, (start, end) = idx_rng
        ls = requests.Session()
        ls.trust_env = use_proxy
        for attempt in range(4):
            try:
                hdrs = {"User-Agent": "Mozilla/5.0", "Range": f"bytes={start}-{end}"}
                with ls.get(url, headers=hdrs, stream=True, timeout=60) as r:
                    if r.status_code not in (200, 206):
                        continue
                    data = r.content
                    if len(data) != (end - start + 1):
                        continue
                    return idx, start, data
            except Exception:
                pass
        raise RuntimeError(f"Chunk {idx} failed for {url}")

    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as ex:
        for idx, start, data in ex.map(dl_chunk, enumerate(ranges)):
            with open(part, "r+b") as f:
                f.seek(start)
                f.write(data)
            print(f"  chunk {idx+1}/{threads} ok ({len(data)/1048576:.1f} MB)", flush=True)

    os.replace(part, dest)
    print(f"  DONE: {os.path.basename(dest)} ({os.path.getsize(dest):,} bytes)", flush=True)
    return True

base = os.path.join(project_root, "01_data")

# Check GSE58812, GSE103091, GSE173839 series matrix listings cleanly
s = requests.Session(); s.trust_env = True
for gse, prefix in [("GSE58812", "GSE58nnn"), ("GSE103091", "GSE103nnn"), ("GSE173839", "GSE173nnn")]:
    u = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/matrix/"
    r = s.get(u, timeout=15)
    print(f"{gse}/matrix len={len(r.text)}, hrefs={re.findall(r'href=\"([^\"]+)\"', r.text)}")

# Download core transcriptomic & ICI cohorts first
targets = [
    ("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FHiSeqV2.gz",
     os.path.join(base, "TCGA_GEO", "TCGA_BRCA_HiSeqV2.tsv.gz"), 8, True),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194040/matrix/GSE194040-GPL30493_series_matrix.txt.gz",
     os.path.join(base, "TCGA_GEO", "GSE194040-GPL30493_series_matrix.txt.gz"), 1, True),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194040/suppl/GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt.gz",
     os.path.join(base, "TCGA_GEO", "GSE194040_geneLevel_n988.txt.gz"), 8, True),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173839/suppl/GSE173839_ISPY2_AgilentGeneExp_durvaPlusCtr_FFPE_probeLevel_n105.txt.gz",
     os.path.join(base, "TCGA_GEO", "GSE173839_durva_probeLevel_n105.txt.gz"), 4, True),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173839/suppl/GSE173839_ISPY2_DurvalumabOlaparibArm_biomarkers.csv.gz",
     os.path.join(base, "TCGA_GEO", "GSE173839_biomarkers.csv.gz"), 1, True),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58812/matrix/GSE58812_series_matrix.txt.gz",
     os.path.join(base, "TCGA_GEO", "GSE58812_series_matrix.txt.gz"), 4, True),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103091/matrix/GSE103091_series_matrix.txt.gz",
     os.path.join(base, "TCGA_GEO", "GSE103091_series_matrix.txt.gz"), 4, True),
]

for url, dest, th, prx in targets:
    try:
        range_download(url, dest, threads=th, use_proxy=prx)
    except Exception as e:
        print("Download error:", dest, e, flush=True)
