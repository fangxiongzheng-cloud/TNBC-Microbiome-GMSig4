import concurrent.futures
import os
import requests

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

def range_download(url, dest, threads=8, use_proxy=True):
    if os.path.exists(dest) and os.path.getsize(dest) > 1000000:
        print(f"SKIP (exists): {dest} ({os.path.getsize(dest):,} bytes)", flush=True)
        return True
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    s = requests.Session(); s.trust_env = use_proxy
    h = s.head(url, timeout=30, allow_redirects=True, headers={"User-Agent": "Mozilla/5.0"})
    if h.status_code != 200:
        print(f"HEAD failed {h.status_code} for {url}", flush=True)
        return False
    size = int(h.headers.get("Content-Length", 0))
    part = dest + ".part"
    print(f"Multi-range ({threads}x) downloading {os.path.basename(dest)} ({size/1048576:.1f} MB)...", flush=True)
    with open(part, "wb") as f:
        f.truncate(size)
    chunk_size = (size + threads - 1) // threads
    ranges = [(i * chunk_size, min(size - 1, (i + 1) * chunk_size - 1)) for i in range(threads)]

    def dl_chunk(idx_rng):
        idx, (start, end) = idx_rng
        ls = requests.Session(); ls.trust_env = use_proxy
        for attempt in range(5):
            try:
                with ls.get(url, headers={"User-Agent": "Mozilla/5.0", "Range": f"bytes={start}-{end}"}, timeout=120) as r:
                    if r.status_code not in (200, 206):
                        continue
                    data = r.content
                    if len(data) == (end - start + 1):
                        return idx, start, data
            except Exception as e:
                pass
        raise RuntimeError(f"Chunk {idx} failed for {url}")

    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as ex:
        for idx, start, data in ex.map(dl_chunk, enumerate(ranges)):
            with open(part, "r+b") as f:
                f.seek(start)
                f.write(data)
            print(f"  [{os.path.basename(dest)}] chunk {idx+1}/{threads} ok", flush=True)

    os.replace(part, dest)
    print(f"DONE: {dest} ({os.path.getsize(dest):,} bytes)", flush=True)
    return True

# 1. Check MAGMA g1000_eur.zip vs MRC 1kg.v3.tgz, and EBI eQTL Catalogue / eqtlgen
s = requests.Session(); s.trust_env = True
for name, url in [
    ("MAGMA_g1000_eur", "https://vu.data.surfsara.nl/index.php/s/VZN6t7vX9kR7sP8/download"),
    ("EBI_eQTLGen_blood", "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/eQTLGen/"),
    ("GCST010100_tbi", "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010100/harmonised/GCST010100.h.tsv.gz.tbi"),
]:
    try:
        r = s.head(url, timeout=15, allow_redirects=True)
        print(f"CHECK {name}: HTTP {r.status_code}, size={r.headers.get('Content-Length')}", flush=True)
    except Exception as e:
        print(f"CHECK {name}: ERR {e}", flush=True)

# 2. Download GCST010100.h.tsv.gz (TNBC GWAS, 365.4 MB) & .tbi
gwas_dir = os.path.join(project_root, "01_data", "GWAS_MR")
range_download(
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010100/harmonised/GCST010100.h.tsv.gz",
    os.path.join(gwas_dir, "GCST010100_TNBC.h.tsv.gz"),
    threads=8, use_proxy=True
)
