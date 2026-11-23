"""Small public metadata only; no raw downloads, scores, or associations.
Usage: python 00_scripts/26e_spatial_proteomic_feasibility.py URL relative_filename [max_MiB]
"""
from pathlib import Path
import requests, sys, json, hashlib, datetime
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'01_data/JTM_upgrade/spatial_metadata'
OUT=ROOT/'02_results/JTM_upgrade/spatial_feasibility'
DATA.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
def fetch(url,name,limit=10*1024**2):
    dest=DATA/name; dest.parent.mkdir(parents=True,exist_ok=True)
    manifest=DATA/'source_manifest.json'
    recs=json.loads(manifest.read_text()) if manifest.exists() else []
    if dest.exists():
        print('EXISTS',name,dest.stat().st_size);return dest
    with requests.get(url,timeout=90,stream=True) as r:
        r.raise_for_status()
        if int(r.headers.get('Content-Length',0))>limit:raise ValueError('Content-Length exceeds cap')
        blob=bytearray()
        for chunk in r.iter_content(65536):
            blob.extend(chunk)
            if len(blob)>limit:raise ValueError('Stream exceeds cap')
        # No attempts to solve or bypass access challenges.
        low=bytes(blob[:20000]).lower()
        if any(x in low for x in [b'verify you are human',b'cf-chl-',b'hcaptcha',b'recaptcha/api',b'recaptcha/challengepage',b'checking your browser']):
            raise RuntimeError('Human/anti-bot challenge: stopped')
        dest.write_bytes(blob)
        recs.append(dict(source_url=url,resolved_url=r.url,file=str(dest.relative_to(ROOT)),bytes=len(blob),sha256=hashlib.sha256(blob).hexdigest(),retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),content_type=r.headers.get('Content-Type')))
        manifest.write_text(json.dumps(recs,indent=2))
        print(name,len(blob),flush=True)
    return dest
if __name__=='__main__': fetch(sys.argv[1],sys.argv[2],int(float(sys.argv[3])*1024**2) if len(sys.argv)>3 else 10*1024**2)
