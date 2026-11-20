import concurrent.futures
import gzip
import json
import os
import requests

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

def verify_gzip(path):
    try:
        with gzip.open(path, "rb") as f:
            while f.read(1024 * 1024):
                pass
        return True
    except Exception as e:
        print(f"GZIP verify failed for {os.path.basename(path)}: {e}", flush=True)
        return False

def download_chunked(url, dest, chunk_mb=4, threads=8, use_proxy=True, check_gz=True):
    if os.path.exists(dest) and os.path.getsize(dest) > 1024:
        if not check_gz or verify_gzip(dest):
            print(f"SKIP (verified): {os.path.basename(dest)} ({os.path.getsize(dest):,} bytes)", flush=True)
            return True
        else:
            os.remove(dest)

    os.makedirs(os.path.dirname(dest), exist_ok=True)
    s = requests.Session()
    s.trust_env = use_proxy
    h = s.head(url, timeout=30, allow_redirects=True, headers={"User-Agent": "Mozilla/5.0"})
    if h.status_code != 200:
        print(f"HEAD failed {h.status_code} for {url}", flush=True)
        return False
    size = int(h.headers.get("Content-Length", 0))
    part = dest + ".part"
    meta = dest + ".meta.json"

    if size < 4 * 1024 * 1024 or h.headers.get("Accept-Ranges", "").lower() != "bytes":
        with s.get(url, stream=True, timeout=180, headers={"User-Agent": "Mozilla/5.0"}) as r:
            r.raise_for_status()
            with open(part, "wb") as f:
                for ch in r.iter_content(512 * 1024):
                    if ch:
                        f.write(ch)
        os.replace(part, dest)
        print(f"DONE (single): {os.path.basename(dest)} ({os.path.getsize(dest):,} bytes)", flush=True)
        return True

    chunk_size = chunk_mb * 1024 * 1024
    ranges = []
    pos = 0
    idx = 0
    while pos < size:
        end = min(size - 1, pos + chunk_size - 1)
        ranges.append((idx, pos, end))
        pos = end + 1
        idx += 1

    done_chunks = set()
    if os.path.exists(meta) and os.path.exists(part) and os.path.getsize(part) == size:
        try:
            with open(meta, "r") as f:
                m = json.load(f)
            if m.get("size") == size and m.get("chunk_size") == chunk_size:
                done_chunks = set(m.get("done", []))
        except Exception:
            pass

    if not done_chunks:
        with open(part, "wb") as f:
            f.truncate(size)

    todo = [r for r in ranges if r[0] not in done_chunks]
    print(f"Chunked ({threads}x, {chunk_mb}MB) {os.path.basename(dest)}: {size/1048576:.1f} MB | {len(done_chunks)}/{len(ranges)} cached, {len(todo)} todo", flush=True)

    def dl_chunk(item):
        cidx, start, end = item
        ls = requests.Session()
        ls.trust_env = use_proxy
        for attempt in range(8):
            try:
                with ls.get(url, headers={"User-Agent": "Mozilla/5.0", "Range": f"bytes={start}-{end}"}, timeout=60) as r:
                    if r.status_code not in (200, 206):
                        continue
                    data = r.content
                    if len(data) == (end - start + 1):
                        return cidx, start, data
            except Exception:
                pass
        raise RuntimeError(f"Chunk {cidx} ({start}-{end}) failed after retries for {url}")

    import threading
    lock = threading.Lock()
    completed = len(done_chunks)
    total = len(ranges)

    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as ex:
        futs = [ex.submit(dl_chunk, item) for item in todo]
        for fut in concurrent.futures.as_completed(futs):
            cidx, start, data = fut.result()
            with lock:
                with open(part, "r+b") as f:
                    f.seek(start)
                    f.write(data)
                done_chunks.add(cidx)
                completed += 1
                if completed % 10 == 0 or completed == total:
                    with open(meta, "w") as mf:
                        json.dump({"size": size, "chunk_size": chunk_size, "done": sorted(done_chunks)}, mf)
                    print(f"  [{os.path.basename(dest)}] {completed}/{total} chunks ({100*completed/total:.1f}%)", flush=True)

    if len(done_chunks) != total:
        raise RuntimeError(f"Incomplete chunks for {dest}: {len(done_chunks)}/{total}")

    if check_gz and not verify_gzip(part):
        raise RuntimeError(f"Corrupt gzip after download: {part}")

    os.replace(part, dest)
    if os.path.exists(meta):
        os.remove(meta)
    print(f"DONE & VERIFIED: {os.path.basename(dest)} ({os.path.getsize(dest):,} bytes)", flush=True)
    return True

base = os.path.join(project_root, "01_data")

# 1. GSE96058 expression (564MB)
download_chunked(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/suppl/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz",
    os.path.join(base, "TCGA_GEO", "GSE96058_gene_expression_n3273.csv.gz"),
    chunk_mb=4, threads=8, check_gz=True
)

# 2. GCST010100 Zhang 2020 TNBC GWAS harmonized (365MB) + .tbi (1.5MB)
download_chunked(
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010100/harmonised/32424353-GCST010100-EFO_0005537-build37.f.tsv.gz",
    os.path.join(base, "GWAS_MR", "GCST010100_TNBC.h.tsv.gz"),
    chunk_mb=4, threads=8, check_gz=True
)
download_chunked(
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010100/harmonised/32424353-GCST010100-EFO_0005537-build37.f.tsv.gz.tbi",
    os.path.join(base, "GWAS_MR", "GCST010100_TNBC.h.tsv.gz.tbi"),
    chunk_mb=4, threads=1, check_gz=False
)

# 3. GSE176078 Wu et al. 2021 BRCA scRNA-seq atlas (533MB)
download_chunked(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz",
    os.path.join(base, "scRNA", "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"),
    chunk_mb=4, threads=8, check_gz=True
)
