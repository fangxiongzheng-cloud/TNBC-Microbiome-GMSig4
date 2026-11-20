import json
import os
import re
import urllib.parse
import requests

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

s = requests.Session()
s.trust_env = False  # Direct connection for domestic HMU/HUST servers

gut_dir = os.path.join(project_root, "01_data", "gutMGene")
gm_dir = os.path.join(project_root, "01_data", "GMrepo")
os.makedirs(gut_dir, exist_ok=True)
os.makedirs(gm_dir, exist_ok=True)

# 1. Discover ALL download file names in gutMGene v2.0 mydata chunks
base_js = "http://bio-computing.hrbmu.edu.cn/gutmgene/static/js/"
r_idx = s.get(base_js + "index-89483b6d.js", timeout=15).text
mydata_chunks = sorted(set(re.findall(r'(?:onepart|twopart|threepart|fourpart)-[A-Za-z0-9]+\.js', r_idx)))
print("Mydata chunks:", mydata_chunks)

download_files = []
for jf in mydata_chunks:
    txt = s.get(base_js + jf, timeout=15).text
    # Look for fileName:"..." and filepath patterns
    for fn in re.findall(r'fileName\s*:\s*["\']([^"\']+)["\']', txt):
        download_files.append((jf, fn))
    for fp in re.findall(r'downloadFolder\?filepath=([^"`\']+)', txt):
        print(f"  {jf} downloadFolder pattern: {fp}")
    for dt in re.findall(r'downloadtable[^;]{0,200}', txt):
        print(f"  {jf} downloadtable: {dt}")

print("Discovered gutMGene files:", download_files)

# Download discovered files via downloadFolder endpoint
for jf, fn in download_files:
    filepath = "allfile/" + fn
    url = "http://bio-computing.hrbmu.edu.cn/gutMGene2.0_api/dow/downloadFolder?filepath=" + urllib.parse.quote(filepath)
    out_path = os.path.join(gut_dir, fn.replace(" ", "_"))
    r = s.get(url, timeout=60)
    print(f"GET {fn} -> {r.status_code}, bytes={len(r.content)}, content-type={r.headers.get('Content-Type')}")
    if r.status_code == 200 and len(r.content) > 100:
        with open(out_path, "wb") as f:
            f.write(r.content)
        # Print first 2 lines
        preview = r.content[:400].decode("utf-8", "ignore").splitlines()[:3]
        print("  preview:", preview)

# 2. Inspect existing projects_D001943.json and test GMrepo API payloads
pfile = os.path.join(gm_dir, "projects_D001943.json")
if os.path.exists(pfile):
    print("\nExisting projects_D001943.json:", open(pfile, encoding="utf-8").read()[:500])

print("\nTesting GMrepo POST parameter styles:")
headers = {"User-Agent": "Mozilla/5.0"}
for style, kw in [
    ("form data=dict", {"data": {"mesh_id": "D001943"}}),
    ("json_str data=dumps + header", {"data": json.dumps({"mesh_id": "D001943"}), "headers": {**headers, "Content-Type": "application/json"}}),
]:
    r = s.post("https://gmrepo.humangut.info/api/getStatisticsByProjectsByMeshID/", timeout=30, **kw)
    print(f"  {style} -> {r.status_code}, len={len(r.text)}, text={r.text[:350]}")

for proj in ["PRJNA658160", "PRJNA383849", "PRJNA872152"]:
    r = s.post(
        "https://gmrepo.humangut.info/api/getMicrobeAbundancesByPhenotypeMeshIDAndProjectID/",
        data={"project_id": proj, "mesh_id": "D001943"},
        headers=headers,
        timeout=30,
    )
    print(f"  Proj {proj} form -> {r.status_code}, len={len(r.text)}, preview={r.text[:250]}")
    if r.status_code == 200 and len(r.text) < 200:
        r2 = s.post(
            "https://gmrepo.humangut.info/api/getMicrobeAbundancesByPhenotypeMeshIDAndProjectID/",
            data=json.dumps({"project_id": proj, "mesh_id": "D001943"}),
            headers={**headers, "Content-Type": "application/json"},
            timeout=30,
        )
        print(f"  Proj {proj} json -> {r2.status_code}, len={len(r2.text)}, preview={r2.text[:250]}")
