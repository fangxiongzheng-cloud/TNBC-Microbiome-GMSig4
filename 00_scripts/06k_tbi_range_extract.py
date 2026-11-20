import concurrent.futures
import csv
import gzip
import os
import struct
import zlib
import requests

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

def reg2bins(beg, end):
    # Standard UCSC / Tabix 0-based half-open [beg, end)
    bins = [0]
    end -= 1
    for shift, offset in [(26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)]:
        for k in range(offset + (beg >> shift), offset + (end >> shift) + 1):
            bins.append(k)
    return bins

class TabixIndex:
    def __init__(self, tbi_path):
        with gzip.open(tbi_path, "rb") as f:
            data = f.read()
        assert data[:4] == b"TBI\x01", f"Bad magic in {tbi_path}"
        off = 4
        self.n_ref, self.fmt, self.col_seq, self.col_beg, self.col_end, self.meta, self.skip, l_nm = struct.unpack_from("<8i", data, off)
        off += 32
        raw_names = data[off:off+l_nm].split(b"\x00")
        self.names = [x.decode("utf-8") for x in raw_names[:self.n_ref]]
        off += l_nm
        self.refs = {}
        for rname in self.names:
            n_bin = struct.unpack_from("<i", data, off)[0]
            off += 4
            bins = {}
            for _ in range(n_bin):
                bid, n_chunk = struct.unpack_from("<Ii", data, off)
                off += 8
                chunks = []
                for _ in range(n_chunk):
                    cbeg, cend = struct.unpack_from("<QQ", data, off)
                    off += 16
                    chunks.append((cbeg, cend))
                bins[bid] = chunks
            n_intv = struct.unpack_from("<i", data, off)[0]
            off += 4
            ioff = list(struct.unpack_from(f"<{n_intv}Q", data, off)) if n_intv > 0 else []
            off += 8 * n_intv
            self.refs[rname] = (bins, ioff)

    def byte_ranges(self, queries, merge_gap=524288, max_span=3145728):
        # queries: list of (chr_name, beg_1based, end_1based)
        spans = []
        for chrom, sbeg, send in queries:
            cname = chrom if chrom in self.refs else ("chr" + chrom if "chr" + chrom in self.refs else chrom.replace("chr", ""))
            if cname not in self.refs:
                continue
            bins, ioff = self.refs[cname]
            b0 = max(0, sbeg - 1)
            e0 = max(b0 + 1, send)
            min_voff = 0
            idx = b0 >> 14
            if ioff:
                min_voff = ioff[min(idx, len(ioff) - 1)]
            for bid in reg2bins(b0, e0):
                if bid in bins:
                    for cbeg, cend in bins[bid]:
                        if cend > min_voff:
                            c0 = cbeg >> 16
                            c1 = (cend >> 16) + 65536
                            spans.append((c0, c1))
        if not spans:
            return []
        spans.sort()
        merged = [spans[0]]
        for s0, s1 in spans[1:]:
            if s0 <= merged[-1][1] + merge_gap and (max(merged[-1][1], s1) - merged[-1][0]) <= max_span:
                merged[-1] = (merged[-1][0], max(merged[-1][1], s1))
            else:
                merged.append((s0, s1))
        return merged

def parse_bgzf_lines(buf):
    pos = 0
    n = len(buf)
    lines = []
    # Find first BGZF header
    while pos + 18 <= n:
        if buf[pos:pos+4] != b"\x1f\x8b\x08\x04":
            nxt = buf.find(b"\x1f\x8b\x08\x04", pos + 1)
            if nxt == -1:
                break
            pos = nxt
        xlen = struct.unpack_from("<H", buf, pos + 10)[0]
        if pos + 12 + xlen > n:
            break
        bsize = None
        sp = pos + 12
        end_sub = sp + xlen
        while sp + 4 <= end_sub:
            si1, si2, slen = struct.unpack_from("<BBH", buf, sp)
            if si1 == 66 and si2 == 67 and slen == 2:
                bsize = struct.unpack_from("<H", buf, sp + 4)[0]
                break
            sp += 4 + slen
        if bsize is None:
            pos += 1
            continue
        blk_len = bsize + 1
        if pos + blk_len > n:
            break
        cdata = buf[pos + 12 + xlen : pos + blk_len - 8]
        if len(cdata) >= 2:
            try:
                raw = zlib.decompress(cdata, -zlib.MAX_WBITS)
                lines.append(raw)
            except Exception:
                pass
        pos += blk_len
    return b"".join(lines).decode("utf-8", errors="replace").splitlines()

def remote_tabix_extract(url, tbi_path, queries, row_filter_fn, out_csv, threads=16, use_proxy=True):
    if os.path.exists(out_csv) and os.path.getsize(out_csv) > 10000:
        print(f"SKIP (already extracted): {os.path.basename(out_csv)}", flush=True)
        return
    s = requests.Session()
    s.trust_env = use_proxy
    h = s.head(url, timeout=30, allow_redirects=True, headers={"User-Agent": "Mozilla/5.0"})
    h.raise_for_status()
    total_size = int(h.headers["Content-Length"])

    # Fetch header from byte 0..32767
    hr = s.get(url, headers={"User-Agent": "Mozilla/5.0", "Range": "bytes=0-32767"}, timeout=30)
    hr.raise_for_status()
    hdr_lines = parse_bgzf_lines(hr.content)
    header = hdr_lines[0].split("\t")

    idx = TabixIndex(tbi_path)
    ranges = idx.byte_ranges(queries, merge_gap=524288, max_span=3145728)
    final_ranges = [(a, min(total_size - 1, b)) for a, b in ranges]
    total_mb = sum(b - a + 1 for a, b in final_ranges) / 1048576.0
    print(f"Fetching {len(final_ranges)} merged BGZF byte ranges ({total_mb:.2f} MB of {total_size/1048576:.1f} MB) from {os.path.basename(url)}...", flush=True)

    def fetch_one(rng):
        a, b = rng
        ls = requests.Session()
        ls.trust_env = use_proxy
        for att in range(6):
            try:
                r = ls.get(url, headers={"User-Agent": "Mozilla/5.0", "Range": f"bytes={a}-{b}"}, timeout=45)
                if r.status_code in (200, 206) and len(r.content) == (b - a + 1):
                    out = []
                    for ln in parse_bgzf_lines(r.content):
                        if not ln or ln.startswith(header[0] + "\t"):
                            continue
                        parts = ln.split("\t")
                        if row_filter_fn(parts, header):
                            out.append(parts)
                    return out
            except Exception:
                pass
        print(f"  WARN failed byte range {a}-{b}", flush=True)
        return []

    seen = set()
    kept = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as ex:
        futs = [ex.submit(fetch_one, r) for r in final_ranges]
        for k, fut in enumerate(concurrent.futures.as_completed(futs), 1):
            for parts in fut.result():
                key = "\t".join(parts[:5])
                if key not in seen:
                    seen.add(key)
                    kept.append(parts)
            if k % 50 == 0 or k == len(final_ranges):
                print(f"  [{os.path.basename(out_csv)}] {k}/{len(final_ranges)} ranges -> {len(kept)} rows", flush=True)

    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(kept)
    print(f"SAVED {out_csv}: {len(kept)} unique matched rows", flush=True)

if __name__ == "__main__":
    data_dir = os.path.join(project_root, "01_data")
    res_dir  = os.path.join(project_root, "02_results")
    m = {}
    with open(os.path.join(res_dir, "IV_rsids_GRCh38_map.csv"), "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            m[row["rsid"]] = (row["chr38"], int(row["pos38"]))
    pos_set = {(c, p) for c, p in m.values()}
    rsid_set = set(m.keys())
    queries = [(c, max(1, p - 1), p + 1) for c, p in m.values()]

    def filter_gwas(parts, hdr):
        if len(parts) < 8:
            return False
        try:
            if (parts[0].replace("chr", ""), int(parts[1])) in pos_set:
                return True
            if "rsid" in hdr:
                ridx = hdr.index("rsid")
                if ridx < len(parts) and parts[ridx] in rsid_set:
                    return True
        except Exception:
            pass
        return False

    # 1. Extract Zhang 2020 TNBC+BRCA1 (GCST010100)
    remote_tabix_extract(
        "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010100/harmonised/GCST010100.h.tsv.gz",
        os.path.join(data_dir, "GWAS_MR", "GCST010100.h.tsv.gz.tbi"),
        queries, filter_gwas,
        os.path.join(data_dir, "GWAS_MR", "GCST010100_TNBC_IVs_extracted.csv"),
        threads=16
    )
    # 2. Extract Zhang 2020 Pure TNBC (GCST90454344)
    remote_tabix_extract(
        "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90454001-GCST90455000/GCST90454344/harmonised/GCST90454344.h.tsv.gz",
        os.path.join(data_dir, "GWAS_MR", "GCST90454344.h.tsv.gz.tbi"),
        queries, filter_gwas,
        os.path.join(data_dir, "GWAS_MR", "GCST90454344_PureTNBC_IVs_extracted.csv"),
        threads=16
    )
    # 3. Extract Zhang 2020 Overall BCAC (GCST010098)
    remote_tabix_extract(
        "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010098/harmonised/GCST010098.h.tsv.gz",
        os.path.join(data_dir, "GWAS_MR", "GCST010098.h.tsv.gz.tbi"),
        queries, filter_gwas,
        os.path.join(data_dir, "GWAS_MR", "GCST010098_OverallBC_IVs_extracted.csv"),
        threads=16
    )
