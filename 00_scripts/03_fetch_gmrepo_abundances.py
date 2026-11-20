import concurrent.futures
import csv
import json
import os
import requests
from threading import local

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

gm_dir = os.path.join(project_root, "01_data", "GMrepo")
projects = ["PRJNA658160", "PRJNA383849", "PRJNA872152"]

thread_local = local()

def get_session():
    if not hasattr(thread_local, "sess"):
        s = requests.Session()
        s.trust_env = False
        thread_local.sess = s
    return thread_local.sess

def fetch_run(row):
    rid = row["run_id"]
    s = get_session()
    for attempt in range(3):
        try:
            r = s.post(
                "https://gmrepo.humangut.info/api/getFullTaxonomicProfileByRunID/",
                data=json.dumps({"run_id": rid}),
                headers={"User-Agent": "Mozilla/5.0"},
                timeout=25,
            )
            if r.status_code == 200:
                js = r.json()
                genus = {g["scientific_name"]: float(g["relative_abundance"]) for g in js.get("genus", []) if g.get("scientific_name")}
                return rid, genus, None
        except Exception as e:
            err = str(e)
    return rid, {}, "failed"

all_meta = []
for proj in projects:
    fp = os.path.join(gm_dir, f"all_runs_in_project_{proj}.tsv")
    lines = [ln for ln in open(fp, encoding="utf-8", errors="ignore") if not ln.startswith("#")]
    rows = [r for r in csv.DictReader(lines, delimiter="\t") if str(r.get("QCstatus")) == "1"]
    print(f"Fetching {proj}: {len(rows)} QC-passed runs...")
    all_meta.extend(rows)

    run_profiles = {}
    failed = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=12) as ex:
        for rid, genus, err in ex.map(fetch_run, rows):
            if err:
                failed.append(rid)
            else:
                run_profiles[rid] = genus

    print(f"  {proj} done: ok={len(run_profiles)}, failed={len(failed)}")

    # Collect all genera present in this project
    genera = sorted({g for prof in run_profiles.values() for g in prof})
    out_mat = os.path.join(gm_dir, f"genus_abundance_{proj}.csv")
    with open(out_mat, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["run_id", "project_id", "sample_id", "disease", "term", "host_age", "sex", "BMI", "country"] + genera)
        for r in rows:
            rid = r["run_id"]
            if rid not in run_profiles:
                continue
            prof = run_profiles[rid]
            w.writerow(
                [rid, r["project_id"], r["sample_id"], r["disease"], r["term"], r.get("host_age", ""), r.get("sex", ""), r.get("BMI", ""), r.get("country", "")]
                + [prof.get(g, 0.0) for g in genera]
            )
    print(f"  Saved {out_mat}: samples={len(run_profiles)}, genera={len(genera)}")
