#!/usr/bin/env python3
"""
Stage C1 — Fetch GSE210616 (Bassiouni et al., Cancer Res 2023) Visium per-GSM files
Strictly downloads ONLY:
  - *_filtered_feature_bc_matrix.h5
  - *_tissue_positions_list.csv.gz
  - *_scalefactors_json.json.gz
Avoids *.cloupe.gz and GSE210616_RAW.tar (37.69 GB).
Also supports --hires-gsm GSM1,GSM2 to fetch tissue_hires_image.png.gz for representative slides.
Logs exact URLs, expected vs observed byte sizes, and SHA256 checksums to provenance manifest.
"""
import argparse
import concurrent.futures
import hashlib
import json
import os
import re
import time
import urllib.request
from pathlib import Path

BASE_DIR = Path(r"D:\Bioinformatics\TNBC_Microbiome")
GSM_BRIEF = BASE_DIR / "01_data" / "JTM_upgrade" / "spatial_metadata" / "GEO" / "GSE210616_gsm_brief.txt"
FILELIST = BASE_DIR / "01_data" / "JTM_upgrade" / "spatial_metadata" / "GEO" / "GSE210616_filelist.txt"
OUT_DIR = BASE_DIR / "01_data" / "JTM_upgrade" / "spatial_gse210616"


def parse_filelist():
    sizes = {}
    with open(FILELIST, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4 and parts[0] == "File":
                sizes[parts[1]] = int(parts[3])
    return sizes


def parse_gsm_brief():
    txt = GSM_BRIEF.read_text(encoding="utf-8")
    blocks = txt.split("^SAMPLE = ")[1:]
    records = []
    for b in blocks:
        gsm = b.splitlines()[0].strip()
        title_m = re.search(r"!Sample_title = (.+)", b)
        title = title_m.group(1).strip() if title_m else ""
        pt_m = re.search(r"patient\s+(\d+),\s*section\s+(\d+)", title, re.I)
        patient_num = int(pt_m.group(1)) if pt_m else -1
        section_num = int(pt_m.group(2)) if pt_m else -1
        cohort = "Reference" if patient_num <= 15 else "Validation"
        suppl_urls = re.findall(r"!Sample_supplementary_file_\d+ = (\S+)", b)
        files = {}
        for u in suppl_urls:
            https_u = u.replace("ftp://ftp.ncbi.nlm.nih.gov/", "https://ftp.ncbi.nlm.nih.gov/")
            fn = https_u.rsplit("/", 1)[-1]
            if fn.endswith("_filtered_feature_bc_matrix.h5"):
                files["h5"] = (fn, https_u)
            elif fn.endswith("_tissue_positions_list.csv.gz"):
                files["positions"] = (fn, https_u)
            elif fn.endswith("_scalefactors_json.json.gz"):
                files["scalefactors"] = (fn, https_u)
            elif fn.endswith("_tissue_hires_image.png.gz"):
                files["hires"] = (fn, https_u)
        sample_code = files["h5"][0].split("_")[1] if "h5" in files else ""
        records.append({
            "gsm": gsm,
            "sample_code": sample_code,
            "patient_num": patient_num,
            "patient_id": f"P{patient_num:02d}",
            "section_num": section_num,
            "cohort": cohort,
            "title": title,
            "files": files,
        })
    return records


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download_one(item):
    gsm, ftype, fn, url, expected_size, dest = item
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and expected_size is not None and dest.stat().st_size == expected_size:
        return {
            "gsm": gsm,
            "file_type": ftype,
            "filename": fn,
            "url": url,
            "expected_bytes": expected_size,
            "actual_bytes": dest.stat().st_size,
            "sha256": sha256_file(dest),
            "status": "cached_verified",
        }
    tmp = dest.with_suffix(dest.suffix + ".part")
    last_err = None
    for attempt in range(1, 4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (TNBC-JTM-Spatial-Audit)"})
            with urllib.request.urlopen(req, timeout=120) as resp, open(tmp, "wb") as out:
                while True:
                    chunk = resp.read(256 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
            actual = tmp.stat().st_size
            if expected_size is not None and actual != expected_size:
                raise IOError(f"Size mismatch for {fn}: got {actual}, expected {expected_size}")
            os.replace(tmp, dest)
            return {
                "gsm": gsm,
                "file_type": ftype,
                "filename": fn,
                "url": url,
                "expected_bytes": expected_size,
                "actual_bytes": actual,
                "sha256": sha256_file(dest),
                "status": "downloaded_verified",
            }
        except Exception as e:
            last_err = str(e)
            if tmp.exists():
                tmp.unlink(missing_ok=True)
            time.sleep(2 * attempt)
    return {
        "gsm": gsm,
        "file_type": ftype,
        "filename": fn,
        "url": url,
        "expected_bytes": expected_size,
        "actual_bytes": -1,
        "sha256": "",
        "status": f"ERROR: {last_err}",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--hires-gsm", type=str, default="", help="Comma-separated GSM IDs to also download tissue_hires_image.png.gz")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sizes = parse_filelist()
    records = parse_gsm_brief()
    print(f"Parsed {len(records)} GSM sections across {len(set(r['patient_id'] for r in records))} unique patients.")

    hires_set = {x.strip() for x in args.hires_gsm.split(",") if x.strip()}
    tasks = []
    for r in records:
        gsm = r["gsm"]
        sec_dir = OUT_DIR / "raw" / gsm
        for ftype in ["h5", "positions", "scalefactors"]:
            fn, url = r["files"][ftype]
            tasks.append((gsm, ftype, fn, url, sizes.get(fn), sec_dir / fn))
        if gsm in hires_set and "hires" in r["files"]:
            fn, url = r["files"]["hires"]
            tasks.append((gsm, "hires", fn, url, sizes.get(fn), sec_dir / fn))

    total_expected = sum(t[4] for t in tasks if t[4] is not None)
    print(f"Queued {len(tasks)} files totaling {total_expected / (1024**2):.2f} MiB (workers={args.workers})...")

    results = []
    t0 = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
        for i, res in enumerate(ex.map(download_one, tasks), 1):
            results.append(res)
            if res["status"].startswith("ERROR") or i % 15 == 0 or i == len(tasks):
                print(f"[{i}/{len(tasks)}] {res['filename']} -> {res['status']} ({res['actual_bytes']} B)")

    errors = [r for r in results if r["status"].startswith("ERROR")]
    elapsed = time.time() - t0
    print(f"Completed {len(results)} files in {elapsed:.1f}s; errors={len(errors)}")

    # Save sample metadata CSV + provenance JSON
    meta_csv = OUT_DIR / "gse210616_gsm_metadata.csv"
    with open(meta_csv, "w", encoding="utf-8") as f:
        f.write("gsm,sample_code,patient_num,patient_id,section_num,cohort,title,h5_file,positions_file,scalefactors_file,hires_file\n")
        for r in records:
            f.write(
                f"{r['gsm']},{r['sample_code']},{r['patient_num']},{r['patient_id']},{r['section_num']},"
                f"{r['cohort']},\"{r['title']}\",{r['files']['h5'][0]},{r['files']['positions'][0]},"
                f"{r['files']['scalefactors'][0]},{r['files']['hires'][0]}\n"
            )

    manifest_path = OUT_DIR / "download_provenance_manifest.json"
    # Merge with existing manifest if hires added later
    existing = {}
    if manifest_path.exists():
        old = json.loads(manifest_path.read_text(encoding="utf-8"))
        for entry in old.get("files", []):
            existing[entry["filename"]] = entry
    for r in results:
        existing[r["filename"]] = r
    manifest = {
        "dataset": "GSE210616",
        "publication": "Bassiouni et al., Cancer Res 2023 (PMID:36283023)",
        "n_sections": len(records),
        "n_patients": len(set(r["patient_id"] for r in records)),
        "total_verified_bytes": sum(v["actual_bytes"] for v in existing.values() if v["actual_bytes"] > 0),
        "n_files": len(existing),
        "n_errors": sum(1 for v in existing.values() if v["status"].startswith("ERROR")),
        "files": sorted(existing.values(), key=lambda x: (x["gsm"], x["file_type"])),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote {meta_csv} and {manifest_path}")
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
