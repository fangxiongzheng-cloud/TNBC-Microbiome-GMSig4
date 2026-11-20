import requests, zlib

s = requests.Session()
s.trust_env = True

files = [
    ("GCST010100_meta", "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010100/harmonised/GCST010100.h.tsv.gz-meta.yaml", False),
    ("GCST90454344_meta", "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90454001-GCST90455000/GCST90454344/harmonised/GCST90454344.h.tsv.gz-meta.yaml", False),
    ("GCST010098_meta", "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010098/harmonised/GCST010098.h.tsv.gz-meta.yaml", False),
    ("GCST010100_head", "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010100/harmonised/GCST010100.h.tsv.gz", True),
    ("GCST90454344_head", "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90454001-GCST90455000/GCST90454344/harmonised/GCST90454344.h.tsv.gz", True),
    ("GTEx_Breast_head", "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/Breast_Mammary_Tissue.tsv.gz", True),
]

for name, url, is_gz in files:
    print(f"\n=== {name} ===")
    try:
        hdrs = {"User-Agent": "Mozilla/5.0"}
        if is_gz:
            hdrs["Range"] = "bytes=0-32767"
        r = s.get(url, headers=hdrs, timeout=20)
        print("Status:", r.status_code, "CR:", r.headers.get("Content-Range"), "CL:", r.headers.get("Content-Length"))
        if is_gz:
            d = zlib.decompressobj(16 + zlib.MAX_WBITS)
            txt = d.decompress(r.content).decode("utf-8", errors="replace")
            for ln in txt.splitlines()[:4]:
                print(" ", ln[:300])
        else:
            print(r.text[:1200])
    except Exception as e:
        print("ERR:", e)
