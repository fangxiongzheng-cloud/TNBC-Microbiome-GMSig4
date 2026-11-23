"""Public processed-data acquisition and eligibility audit only; never computes risk scores.
Run:
  python 00_scripts/26b_download_processed.py
  python 00_scripts/26b_fetch_reference_evidence.py
  python 00_scripts/26b_external_cohort_eligibility.py audit
All transfers capped at 500 MiB; no controlled data and no raw assay files.
"""
from pathlib import Path
import urllib.request
import json
import hashlib
import datetime
import time
import sys

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / '01_data/JTM_upgrade/external'
OUT = ROOT / '02_results/JTM_upgrade/cohort_eligibility'
DATA.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)
MAN = DATA / 'source_manifest.json'

TARGET_GENES = ['IFNG', 'IDO1', 'CXCL10', 'ENO1', 'STAT1', 'PDCD1', 'CD274', 'GZMA', 'PRF1']
GEP_18_GENES = [
    'PSMB10', 'HLA-DQA1', 'HLA-DRB1', 'CMKLR1', 'HLA-E', 'NKG7',
    'CD8A', 'CCL5', 'CXCL9', 'CD27', 'CXCR6', 'IDO1',
    'STAT1', 'TIGIT', 'LAG3', 'CD274', 'PDCD1LG2', 'CD276'
]
ALL_PANEL_GENES = TARGET_GENES + [g for g in GEP_18_GENES if g not in TARGET_GENES]


def _load_manifest():
    return json.loads(MAN.read_text(encoding='utf-8')) if MAN.exists() else []


def _save_manifest(records):
    MAN.write_text(json.dumps(records, indent=2), encoding='utf-8')


def _append_manifest(entry):
    records = _load_manifest()
    records = [r for r in records if r.get('path') != entry.get('path')]
    records.append(entry)
    _save_manifest(records)


def fetch(url, name, limit=500 * 1024 ** 2, terms='Public endpoint; source-specific terms documented in eligibility report.'):
    dest = DATA / name
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    last_err = None
    for attempt in range(4):
        try:
            req = urllib.request.Request(
                url,
                headers={'User-Agent': 'Academic reproducibility data audit; urllib Python'}
            )
            with urllib.request.urlopen(req, timeout=180) as r:
                n = int(r.headers.get('Content-Length') or 0)
                if n > limit:
                    raise ValueError(f'Size cap: {n} > {limit}')
                blob = bytearray()
                while True:
                    b = r.read(1024 * 1024)
                    if not b:
                        break
                    blob.extend(b)
                    if len(blob) > limit:
                        raise ValueError('Size cap reached')
                dest.write_bytes(blob)
                _append_manifest(dict(
                    url=url,
                    resolved_url=r.url,
                    path=str(dest.relative_to(ROOT)),
                    retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    bytes=len(blob),
                    sha256=hashlib.sha256(blob).hexdigest(),
                    content_type=r.headers.get('Content-Type'),
                    terms=terms
                ))
                print(name, len(blob), flush=True)
                return dest
        except Exception as e:
            last_err = e
            if attempt < 3:
                time.sleep(2 ** attempt)
    raise last_err


def fetch_api_json(url, name, payload=None, limit=100 * 1024 ** 2, terms='cBioPortal REST API public endpoint.'):
    dest = DATA / name
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    last_err = None
    for attempt in range(4):
        try:
            headers = {
                'User-Agent': 'Academic reproducibility data audit; urllib Python',
                'Accept': 'application/json',
            }
            data_bytes = None
            if payload is not None:
                data_bytes = json.dumps(payload).encode('utf-8')
                headers['Content-Type'] = 'application/json'
            req = urllib.request.Request(url, data=data_bytes, headers=headers)
            with urllib.request.urlopen(req, timeout=180) as r:
                blob = r.read(limit + 1)
                if len(blob) > limit:
                    raise ValueError('API size cap reached')
                dest.write_bytes(blob)
                _append_manifest(dict(
                    url=url,
                    resolved_url=r.url,
                    method='POST' if payload is not None else 'GET',
                    path=str(dest.relative_to(ROOT)),
                    retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    bytes=len(blob),
                    sha256=hashlib.sha256(blob).hexdigest(),
                    content_type=r.headers.get('Content-Type'),
                    terms=terms
                ))
                print(name, len(blob), flush=True)
                return dest
        except Exception as e:
            last_err = e
            if attempt < 3:
                time.sleep(2 ** attempt)
    raise last_err


def fetch_initial_metadata():
    for gse in ['GSE25066', 'GSE21653', 'GSE103091', 'GSE58812']:
        fetch(
            f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&targ=self&form=text&view=full',
            f'{gse}/{gse}_metadata.txt',
            10 * 1024 ** 2
        )
    fetch(
        'https://api.github.com/repos/cBioPortal/datahub/contents/public/brca_metabric',
        'METABRIC/datahub_directory.json',
        10 * 1024 ** 2
    )
    for fn in ['data_clinical_patient.txt', 'data_clinical_sample.txt', 'meta_study.txt', 'meta_mrna_illumina_microarray.txt']:
        fetch(
            'https://raw.githubusercontent.com/cBioPortal/datahub/master/public/brca_metabric/' + fn,
            'METABRIC/' + fn,
            20 * 1024 ** 2
        )


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'fetch'
    if cmd == 'fetch':
        fetch_initial_metadata()
    else:
        raise SystemExit(f'Unknown command: {cmd}. Audit module will be invoked after downloads.')
