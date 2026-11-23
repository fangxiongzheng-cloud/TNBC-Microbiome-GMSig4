#!/usr/bin/env python3
"""
Stage C2: Download brca_cptac_2020 (Krug et al., Cell 2020) clinical, RNA-seq FPKM,
and mass-spec proteomics files from cBioPortal datahub GitHub LFS and record
exact provenance (URL, byte size, SHA256) in cptac_download_manifest.csv.
"""

from __future__ import annotations

import csv
import hashlib
import os
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

BASE_DIR = Path(r"D:\Bioinformatics\TNBC_Microbiome")
DATA_DIR = BASE_DIR / "01_data" / "JTM_upgrade" / "cptac"
OUT_DIR = BASE_DIR / "02_results" / "JTM_upgrade" / "cptac"

BASE_URL = "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/brca_cptac_2020"
FILES = [
    "data_clinical_patient.txt",
    "data_clinical_sample.txt",
    "data_mrna_seq_fpkm.txt",
    "data_protein_quantification.txt",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download_file(url: str, dest: Path) -> tuple[int, int, str]:
    req = urllib.request.Request(url, headers={"User-Agent": "Antigravity-JTM-CPTAC-Validation/1.0"})
    with urllib.request.urlopen(req, timeout=180) as resp:
        status = getattr(resp, "status", 200)
        h = hashlib.sha256()
        total = 0
        tmp = dest.with_suffix(dest.suffix + ".part")
        with tmp.open("wb") as out:
            while True:
                chunk = resp.read(1024 * 1024)
                if not chunk:
                    break
                out.write(chunk)
                h.update(chunk)
                total += len(chunk)
        os.replace(tmp, dest)
    return status, total, h.hexdigest()


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows = []
    for fname in FILES:
        url = f"{BASE_URL}/{fname}"
        dest = DATA_DIR / fname
        ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        status, nbytes, digest = download_file(url, dest)
        # Re-verify SHA256 from disk
        disk_digest = sha256_file(dest)
        assert disk_digest == digest, f"SHA256 mismatch for {fname}"
        rows.append(
            {
                "filename": fname,
                "url": url,
                "local_path": str(dest),
                "http_status": status,
                "bytes": nbytes,
                "mib": f"{nbytes / (1024 * 1024):.4f}",
                "sha256": disk_digest,
                "downloaded_utc": ts,
            }
        )
        print(f"[OK] {fname}: {nbytes:,} bytes ({nbytes / (1024*1024):.2f} MiB) | SHA256={disk_digest}")

    manifest_path = OUT_DIR / "cptac_download_manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["filename", "url", "local_path", "http_status", "bytes", "mib", "sha256", "downloaded_utc"],
        )
        writer.writeheader()
        writer.writerows(rows)

    total_bytes = sum(r["bytes"] for r in rows)
    print(f"Total downloaded: {total_bytes:,} bytes ({total_bytes / (1024*1024):.2f} MiB)")
    print(f"Manifest saved: {manifest_path}")


if __name__ == "__main__":
    main()
