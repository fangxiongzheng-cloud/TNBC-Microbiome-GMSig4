import re
import requests

s = requests.Session()
s.trust_env = True

urls = [
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010100/harmonised/",
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010100/",
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010098/harmonised/",
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90454001-GCST90455000/GCST90454344/harmonised/",
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90454001-GCST90455000/GCST90454344/",
]

for u in urls:
    try:
        r = s.get(u, timeout=20)
        print("===", u, r.status_code)
        for m in re.findall(r'href="([^"]+\.(?:tsv|txt|gz|tbi|yaml)[^"]*)"', r.text):
            print(" ", m)
    except Exception as e:
        print("ERR", u, e)
