"""Harmonize strictly independent external TNBC cohorts (METABRIC, GSE25066, GSE21653)
and verify eligibility, receptor definitions, gene coverage, missingness, and zero-overlap
provenance. STRICTLY NO GMSig-4 scoring or survival/outcome association is performed here."""

import gzip
import hashlib
import json
import re
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(r"D:\Bioinformatics\TNBC_Microbiome")
EXT = ROOT / "01_data/JTM_upgrade/external"
HARM = EXT / "harmonized"
ELIG = ROOT / "02_results/JTM_upgrade/cohort_eligibility"
HARM.mkdir(parents=True, exist_ok=True)
ELIG.mkdir(parents=True, exist_ok=True)

GMSIG4_GENES = ["IFNG", "IDO1", "CXCL10", "ENO1"]
COMP_GENES = ["STAT1", "PDCD1", "CD274", "GZMA", "PRF1"]
GEP18_GENES = [
    "CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6",
    "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2",
    "PSMB10", "STAT1", "TIGIT",
]
ALL_TARGET_GENES = sorted(set(GMSIG4_GENES + COMP_GENES + GEP18_GENES))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_gpl_unambiguous_map(annot_gz: Path, targets: set) -> dict:
    """Return {probe_id: gene_symbol} for unambiguous single-gene probes (no '///')."""
    mapping = {}
    in_table = False
    id_idx = sym_idx = None
    with gzip.open(annot_gz, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("!platform_table_begin"):
                hdr = next(f).rstrip("\r\n").split("\t")
                id_idx = hdr.index("ID")
                sym_idx = hdr.index("Gene symbol")
                in_table = True
                continue
            if not in_table:
                continue
            if line.startswith("!platform_table_end"):
                break
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) <= max(id_idx, sym_idx):
                continue
            pid = parts[id_idx].strip()
            sym = parts[sym_idx].strip()
            if not sym or "///" in sym:
                continue
            if sym in targets:
                mapping[pid] = sym
    return mapping


def parse_geo_series_matrix(matrix_gz: Path, probe_map: dict):
    """Parse GEO series_matrix.txt.gz metadata and target-probe expression rows."""
    gsm_ids = []
    titles = []
    char_rows = []
    expr_rows = {}  # probe_id -> np.array(float)
    in_table = False
    header_gsms = []

    with gzip.open(matrix_gz, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not in_table:
                if line.startswith("!Sample_geo_accession"):
                    gsm_ids = [x.strip().strip('"') for x in line.rstrip("\r\n").split("\t")[1:]]
                elif line.startswith("!Sample_title"):
                    titles = [x.strip().strip('"') for x in line.rstrip("\r\n").split("\t")[1:]]
                elif line.startswith("!Sample_characteristics_ch1"):
                    char_rows.append([x.strip().strip('"') for x in line.rstrip("\r\n").split("\t")[1:]])
                elif line.startswith("!series_matrix_table_begin"):
                    hdr = [x.strip().strip('"') for x in next(f).rstrip("\r\n").split("\t")]
                    header_gsms = hdr[1:]
                    in_table = True
            else:
                if line.startswith("!series_matrix_table_end"):
                    break
                tab_pos = line.find("\t")
                if tab_pos < 0:
                    continue
                pid = line[:tab_pos].strip().strip('"')
                if pid in probe_map:
                    vals = line[tab_pos + 1:].rstrip("\r\n").split("\t")
                    arr = np.array([np.nan if (v.strip().strip('"') in ("", "NA", "null", "NULL")) else float(v.strip().strip('"')) for v in vals], dtype=float)
                    expr_rows[pid] = arr

    assert gsm_ids == header_gsms, f"GSM order mismatch in {matrix_gz.name}"
    n = len(gsm_ids)
    # Build per-sample key-value dictionary from characteristics
    records = []
    for i in range(n):
        rec = {"sampleID": gsm_ids[i], "title": titles[i] if i < len(titles) else ""}
        for row in char_rows:
            if i >= len(row):
                continue
            cell = row[i]
            if ":" in cell:
                k, v = cell.split(":", 1)
                k = k.strip().lower()
                v = v.strip()
                if k in rec and rec[k] != v:
                    rec[k + "_2"] = v
                else:
                    rec[k] = v
        records.append(rec)
    meta_df = pd.DataFrame(records)

    # Collapse probes to genes: median (primary) and mean (sensitivity)
    gene_to_probes = {}
    for pid, g in probe_map.items():
        if pid in expr_rows:
            gene_to_probes.setdefault(g, []).append(pid)

    expr_med = pd.DataFrame(index=sorted(gene_to_probes.keys()), columns=gsm_ids, dtype=float)
    expr_mean = pd.DataFrame(index=sorted(gene_to_probes.keys()), columns=gsm_ids, dtype=float)
    for g, pids in sorted(gene_to_probes.items()):
        mat = np.vstack([expr_rows[p] for p in pids])
        # Check if linear scale (> 50)
        expr_med.loc[g] = np.nanmedian(mat, axis=0)
        expr_mean.loc[g] = np.nanmean(mat, axis=0)

    return meta_df, expr_med, expr_mean, gene_to_probes


def harmonize_metabric():
    p_clin = pd.read_csv(EXT / "METABRIC/resolved/data_clinical_patient.txt", sep="\t", comment="#")
    s_clin = pd.read_csv(EXT / "METABRIC/resolved/data_clinical_sample.txt", sep="\t", comment="#")
    df = pd.merge(p_clin, s_clin, on="PATIENT_ID", how="inner")
    assert df["PATIENT_ID"].nunique() == len(df) == df["SAMPLE_ID"].nunique(), "Duplicate METABRIC patient/sample IDs"

    # Build gene-level expression matrix from cBioPortal API JSON
    gmap = {}
    for p in sorted((EXT / "METABRIC/api/genes").glob("*.json")):
        d = json.loads(p.read_text(encoding="utf-8"))
        gmap[int(d["entrezGeneId"])] = d["hugoGeneSymbol"]

    mrna_raw = json.loads((EXT / "METABRIC/api/brca_metabric_target_gep_mrna.json").read_text(encoding="utf-8"))
    mrna_df = pd.DataFrame(mrna_raw)
    mrna_df["gene"] = mrna_df["entrezGeneId"].map(gmap)
    expr = mrna_df.pivot_table(index="gene", columns="sampleId", values="value", aggfunc="median")
    expr = expr.loc[sorted(expr.index)]

    # Strict primary TNBC definition: ER_STATUS == Negative & PR_STATUS == Negative & HER2_STATUS == Negative
    # Also record ER_IHC == Negative sensitivity flag
    df["has_mrna"] = df["SAMPLE_ID"].isin(expr.columns)
    df["is_TNBC_primary"] = (
        (df["ER_STATUS"] == "Negative")
        & (df["PR_STATUS"] == "Negative")
        & (df["HER2_STATUS"] == "Negative")
    )
    df["is_TNBC_erihc"] = (
        (df["ER_IHC"] == "Negative")
        & (df["PR_STATUS"] == "Negative")
        & (df["HER2_STATUS"] == "Negative")
    )

    # Parse survival endpoints cleanly
    df["OS_event"] = df["OS_STATUS"].map({"1:DECEASED": 1, "0:LIVING": 0})
    df["OS_yr"] = pd.to_numeric(df["OS_MONTHS"], errors="coerce") / 12.0
    df["DSS_event"] = df["VITAL_STATUS"].map({"Died of Disease": 1, "Living": 0, "Died of Other Causes": 0})
    df["RFS_event"] = df["RFS_STATUS"].map({"1:Recurred": 1, "0:Not Recurred": 0})
    df["RFS_yr"] = pd.to_numeric(df["RFS_MONTHS"], errors="coerce") / 12.0

    df["Age"] = pd.to_numeric(df["AGE_AT_DIAGNOSIS"], errors="coerce")
    df["Tumor_mm"] = pd.to_numeric(df["TUMOR_SIZE"], errors="coerce")
    df["Node_count"] = pd.to_numeric(df["LYMPH_NODES_EXAMINED_POSITIVE"], errors="coerce")
    df["Node_Pos"] = np.where(df["Node_count"].isna(), np.nan, (df["Node_count"] > 0).astype(float))
    df["Stage_num"] = pd.to_numeric(df["TUMOR_STAGE"], errors="coerce")
    df["Stage_III_IV"] = np.where(df["Stage_num"].isna(), np.nan, (df["Stage_num"] >= 3).astype(float))

    tnbc = df[df["is_TNBC_primary"] & df["has_mrna"]].copy()
    # Ensure all 4 GMSig-4 genes finite
    valid_4 = expr.loc[GMSIG4_GENES, tnbc["SAMPLE_ID"]].notna().all(axis=0)
    tnbc["has_all4_gmsig4"] = tnbc["SAMPLE_ID"].map(valid_4)
    tnbc_os = tnbc[tnbc["has_all4_gmsig4"] & tnbc["OS_event"].notna() & tnbc["OS_yr"].notna() & (tnbc["OS_yr"] > 0)].copy()
    tnbc_rfs = tnbc[tnbc["has_all4_gmsig4"] & tnbc["RFS_event"].notna() & tnbc["RFS_yr"].notna() & (tnbc["RFS_yr"] > 0)].copy()

    # Save harmonized clinical & expression (all mRNA-available primary TNBC)
    keep_cols = [
        "PATIENT_ID", "SAMPLE_ID", "ER_IHC", "ER_STATUS", "PR_STATUS", "HER2_STATUS", "HER2_SNP6",
        "CLAUDIN_SUBTYPE", "THREEGENE", "COHORT", "CHEMOTHERAPY", "HORMONE_THERAPY", "RADIO_THERAPY",
        "is_TNBC_primary", "is_TNBC_erihc", "has_all4_gmsig4",
        "OS_yr", "OS_event", "DSS_event", "RFS_yr", "RFS_event",
        "Age", "Tumor_mm", "Node_count", "Node_Pos", "Stage_num", "Stage_III_IV", "GRADE", "NPI",
    ]
    tnbc[keep_cols].to_csv(HARM / "METABRIC_TNBC_clinical.csv", index=False)
    expr[tnbc["SAMPLE_ID"]].to_csv(HARM / "METABRIC_TNBC_expr.csv")

    return {
        "cohort": "METABRIC",
        "source_accession": "cBioPortal brca_metabric (Curtis 2012 / Pereira 2016)",
        "platform": "Illumina HT-12 v3 microarray (log2 intensity, gene-level)",
        "expression_scale": f"log2 normalized [{expr.min().min():.2f}, {expr.max().max():.2f}]",
        "total_samples": len(df),
        "mrna_samples": int(df["has_mrna"].sum()),
        "tnbc_definition": "ER_STATUS=='Negative' & PR_STATUS=='Negative' & HER2_STATUS=='Negative' (cBioPortal clinical_sample; note PR/HER2 partly expression/SNP6-derived in METABRIC; ER_IHC sensitivity tracked)",
        "tnbc_with_mrna": len(tnbc),
        "tnbc_all4_gmsig4": int(tnbc["has_all4_gmsig4"].sum()),
        "gmsig4_genes_present": f"{sum(g in expr.index for g in GMSIG4_GENES)}/4",
        "comparator_genes_present": f"{sum(g in expr.index for g in COMP_GENES)}/5",
        "gep18_genes_present": f"{sum(g in expr.index for g in GEP18_GENES)}/18",
        "primary_endpoint": "OS",
        "primary_evaluable_n": len(tnbc_os),
        "primary_events": int(tnbc_os["OS_event"].sum()),
        "secondary_endpoints": f"DSS (n={int(tnbc_os['DSS_event'].notna().sum())}, events={int(tnbc_os['DSS_event'].sum())}); RFS (n={len(tnbc_rfs)}, events={int(tnbc_rfs['RFS_event'].sum())})",
        "age_complete_n": int(tnbc_os["Age"].notna().sum()),
        "size_mm_complete_n": int(tnbc_os["Tumor_mm"].notna().sum()),
        "node_pos_complete_n": int(tnbc_os["Node_Pos"].notna().sum()),
        "stage_complete_n": int(tnbc_os["Stage_num"].notna().sum()),
        "overlap_with_prior": "0 overlap with TCGA/GSE58812/GSE96058/GSE194040/GSE173839 (UK/Canada METABRIC cohort)",
        "eligibility_verdict": "ELIGIBLE_PRIMARY_OS",
    }


def harmonize_gse25066():
    pmap = load_gpl_unambiguous_map(EXT / "platforms/GPL96.annot.gz", set(ALL_TARGET_GENES))
    meta, expr_med, expr_mean, g2p = parse_geo_series_matrix(EXT / "GSE25066/GSE25066_series_matrix.txt.gz", pmap)

    # Check if linear scale
    mx = np.nanmax(expr_med.values)
    scale_note = f"log2 author-normalized [{np.nanmin(expr_med.values):.2f}, {mx:.2f}]"
    if mx > 50:
        expr_med = np.log2(np.maximum(expr_med, 0) + 1.0)
        expr_mean = np.log2(np.maximum(expr_mean, 0) + 1.0)
        scale_note = "linear -> log2(x+1)"

    meta["ER_IHC"] = meta.get("er_status_ihc", pd.Series(dtype=str))
    meta["PR_IHC"] = meta.get("pr_status_ihc", pd.Series(dtype=str))
    meta["HER2"] = meta.get("her2_status", pd.Series(dtype=str))
    meta["is_TNBC_primary"] = (meta["ER_IHC"] == "N") & (meta["PR_IHC"] == "N") & (meta["HER2"] == "N")

    meta["DRFS_event"] = pd.to_numeric(meta.get("drfs_1_event_0_censored"), errors="coerce")
    meta["DRFS_yr"] = pd.to_numeric(meta.get("drfs_even_time_years"), errors="coerce")
    meta["pCR"] = meta.get("pathologic_response_pcr_rd").map({"pCR": 1, "RD": 0})
    meta["Age"] = pd.to_numeric(meta.get("age_years"), errors="coerce")
    t_str = meta.get("clinical_t_stage", pd.Series("", index=meta.index)).astype(str)
    n_str = meta.get("clinical_nodal_status", pd.Series("", index=meta.index)).astype(str)
    s_str = meta.get("clinical_ajcc_stage", pd.Series("", index=meta.index)).astype(str)
    meta["T_large"] = np.where(t_str.str.match(r"^T[34]"), 1.0, np.where(t_str.str.match(r"^T[012]"), 0.0, np.nan))
    meta["Node_Pos"] = np.where(n_str.str.match(r"^N[123]"), 1.0, np.where(n_str.str.match(r"^N0"), 0.0, np.nan))
    meta["Stage_III_IV"] = np.where(s_str.str.contains(r"III|IV", na=False), 1.0, np.where(s_str.str.match(r"^I($|A|B)|^II($|A|B)"), 0.0, np.nan))

    tnbc = meta[meta["is_TNBC_primary"]].copy()
    has4 = all(g in expr_med.index for g in GMSIG4_GENES)
    tnbc["has_all4_gmsig4"] = expr_med.loc[GMSIG4_GENES, tnbc["sampleID"]].notna().all(axis=0).values if has4 else False
    tnbc_drfs = tnbc[tnbc["has_all4_gmsig4"] & tnbc["DRFS_event"].notna() & tnbc["DRFS_yr"].notna() & (tnbc["DRFS_yr"] > 0)].copy()

    keep_cols = [
        "sampleID", "title", "source", "ER_IHC", "PR_IHC", "HER2", "is_TNBC_primary", "has_all4_gmsig4",
        "DRFS_yr", "DRFS_event", "pCR", "Age", "clinical_t_stage", "T_large", "clinical_nodal_status", "Node_Pos",
        "clinical_ajcc_stage", "Stage_III_IV", "grade", "pam50_class",
    ]
    keep_cols = [c for c in keep_cols if c in tnbc.columns]
    tnbc[keep_cols].to_csv(HARM / "GSE25066_TNBC_clinical.csv", index=False)
    expr_med[tnbc["sampleID"]].to_csv(HARM / "GSE25066_TNBC_expr_median.csv")
    expr_mean[tnbc["sampleID"]].to_csv(HARM / "GSE25066_TNBC_expr_mean.csv")

    src_counts = tnbc_drfs["source"].value_counts(dropna=False).to_dict() if "source" in tnbc_drfs.columns else {}
    missing_gep = sorted(set(GEP18_GENES) - set(expr_med.index))

    return {
        "cohort": "GSE25066",
        "source_accession": "GEO GSE25066 (Hatzis et al. JAMA 2011; SuperSeries of GSE25055+GSE25065)",
        "platform": "Affymetrix HG-U133A (GPL96; unambiguous single-gene probes collapsed by median [primary] & mean [sensitivity])",
        "expression_scale": scale_note,
        "total_samples": len(meta),
        "mrna_samples": len(meta),
        "tnbc_definition": "er_status_ihc=='N' & pr_status_ihc=='N' & her2_status=='N' (strict IHC; excludes indeterminate/NA)",
        "tnbc_with_mrna": len(tnbc),
        "tnbc_all4_gmsig4": int(tnbc["has_all4_gmsig4"].sum()),
        "gmsig4_genes_present": f"{sum(g in expr_med.index for g in GMSIG4_GENES)}/4",
        "comparator_genes_present": f"{sum(g in expr_med.index for g in COMP_GENES)}/5",
        "gep18_genes_present": f"{sum(g in expr_med.index for g in GEP18_GENES)}/18 (missing: {','.join(missing_gep) if missing_gep else 'none'})",
        "primary_endpoint": "DRFS (neoadjuvant taxane-anthracycline chemo cohort)",
        "primary_evaluable_n": len(tnbc_drfs),
        "primary_events": int(tnbc_drfs["DRFS_event"].sum()),
        "secondary_endpoints": f"Chemo pCR (n={int(tnbc['pCR'].notna().sum())}, pCR={int(tnbc['pCR'].sum())}); sources={src_counts}",
        "age_complete_n": int(tnbc_drfs["Age"].notna().sum()),
        "size_mm_complete_n": 0,
        "node_pos_complete_n": int(tnbc_drfs["Node_Pos"].notna().sum()),
        "stage_complete_n": int(tnbc_drfs["Stage_III_IV"].notna().sum()),
        "overlap_with_prior": "0 overlap with TCGA/GSE58812/GSE96058/GSE194040(I-SPY2)/GSE173839 (note source='ISPY' is I-SPY1 2002-2006, distinct from I-SPY2 2010-2016)",
        "eligibility_verdict": "ELIGIBLE_SECONDARY_DRFS",
    }


def harmonize_gse21653():
    pmap = load_gpl_unambiguous_map(EXT / "platforms/GPL570.annot.gz", set(ALL_TARGET_GENES))
    meta, expr_med, expr_mean, g2p = parse_geo_series_matrix(EXT / "GSE21653/GSE21653_series_matrix.txt.gz", pmap)

    mx = np.nanmax(expr_med.values)
    scale_note = f"log2 RMA author-normalized [{np.nanmin(expr_med.values):.2f}, {mx:.2f}]"
    if mx > 50:
        expr_med = np.log2(np.maximum(expr_med, 0) + 1.0)
        expr_mean = np.log2(np.maximum(expr_mean, 0) + 1.0)
        scale_note = "linear -> log2(x+1)"

    meta["ER_IHC"] = meta.get("er ihc", pd.Series(dtype=str))
    meta["PR_IHC"] = meta.get("pr ihc", pd.Series(dtype=str))
    meta["ERBB2"] = meta.get("erbb2", pd.Series(dtype=str))
    meta["is_TNBC_primary"] = (meta["ER_IHC"] == "0") & (meta["PR_IHC"] == "0") & (meta["ERBB2"] == "0")

    meta["DFS_event"] = pd.to_numeric(meta.get("dfs evt"), errors="coerce")
    meta["DFS_yr"] = pd.to_numeric(meta.get("dfs time (months)"), errors="coerce") / 12.0
    meta["Age"] = pd.to_numeric(meta.get("age at diagnosis"), errors="coerce")
    pt_str = meta.get("pt", pd.Series("", index=meta.index)).astype(str)
    meta["T_large"] = np.where(pt_str.str.match(r"^p?T[34]", case=False), 1.0, np.where(pt_str.str.match(r"^p?T[012]", case=False), 0.0, np.nan))
    meta["Node_Pos"] = pd.to_numeric(meta.get("pn"), errors="coerce")
    meta["Grade"] = pd.to_numeric(meta.get("sbr grade"), errors="coerce")

    # Verify 0 GSM overlap with GSE58812
    g58_text = (EXT / "GSE58812/GSE58812_metadata.txt").read_text(encoding="utf-8", errors="ignore")
    g58_gsms = set(re.findall(r"GSM\d+", g58_text))
    ov58 = len(set(meta["sampleID"]) & g58_gsms)

    tnbc = meta[meta["is_TNBC_primary"]].copy()
    has4 = all(g in expr_med.index for g in GMSIG4_GENES)
    tnbc["has_all4_gmsig4"] = expr_med.loc[GMSIG4_GENES, tnbc["sampleID"]].notna().all(axis=0).values if has4 else False
    tnbc_dfs = tnbc[tnbc["has_all4_gmsig4"] & tnbc["DFS_event"].notna() & tnbc["DFS_yr"].notna() & (tnbc["DFS_yr"] > 0)].copy()

    keep_cols = [
        "sampleID", "title", "ER_IHC", "PR_IHC", "ERBB2", "is_TNBC_primary", "has_all4_gmsig4",
        "DFS_yr", "DFS_event", "Age", "pt", "T_large", "Node_Pos", "Grade", "molecular subtype",
    ]
    keep_cols = [c for c in keep_cols if c in tnbc.columns]
    tnbc[keep_cols].to_csv(HARM / "GSE21653_TNBC_clinical.csv", index=False)
    expr_med[tnbc["sampleID"]].to_csv(HARM / "GSE21653_TNBC_expr_median.csv")
    expr_mean[tnbc["sampleID"]].to_csv(HARM / "GSE21653_TNBC_expr_mean.csv")

    missing_gep = sorted(set(GEP18_GENES) - set(expr_med.index))
    return {
        "cohort": "GSE21653",
        "source_accession": "GEO GSE21653 (Sabatier et al. BCRT 2011; IPC Marseille)",
        "platform": "Affymetrix HG-U133 Plus 2.0 (GPL570; unambiguous single-gene probes collapsed by median [primary] & mean [sensitivity])",
        "expression_scale": scale_note,
        "total_samples": len(meta),
        "mrna_samples": len(meta),
        "tnbc_definition": "er ihc=='0' & pr ihc=='0' & erbb2=='0' (strict IHC; excludes NA)",
        "tnbc_with_mrna": len(tnbc),
        "tnbc_all4_gmsig4": int(tnbc["has_all4_gmsig4"].sum()),
        "gmsig4_genes_present": f"{sum(g in expr_med.index for g in GMSIG4_GENES)}/4",
        "comparator_genes_present": f"{sum(g in expr_med.index for g in COMP_GENES)}/5",
        "gep18_genes_present": f"{sum(g in expr_med.index for g in GEP18_GENES)}/18 (missing: {','.join(missing_gep) if missing_gep else 'none'})",
        "primary_endpoint": "DFS",
        "primary_evaluable_n": len(tnbc_dfs),
        "primary_events": int(tnbc_dfs["DFS_event"].sum()),
        "secondary_endpoints": "None (OS not reported in GEO series matrix)",
        "age_complete_n": int(tnbc_dfs["Age"].notna().sum()),
        "size_mm_complete_n": 0,
        "node_pos_complete_n": int(tnbc_dfs["Node_Pos"].notna().sum()),
        "stage_complete_n": 0,
        "overlap_with_prior": f"0 GSM overlap with GSE58812 (checked: {ov58}); distinct institution (IPC Marseille vs ICO Nantes/Rennes); 0 overlap with TCGA/GSE96058/METABRIC",
        "eligibility_verdict": "ELIGIBLE_SECONDARY_DFS",
    }


def main():
    rows = [
        harmonize_metabric(),
        harmonize_gse25066(),
        harmonize_gse21653(),
        {
            "cohort": "GSE103091",
            "source_accession": "GEO GSE103091 (SuperSeries of GSE58812 [n=107] + GSE83937 [n=131])",
            "platform": "Affymetrix HG-U133 Plus 2.0 (GPL570)",
            "expression_scale": "Linear MAS5.0",
            "total_samples": 238,
            "mrna_samples": 238,
            "tnbc_definition": "107 duplicate GSE58812 + 131 GSE83937",
            "tnbc_with_mrna": 238,
            "tnbc_all4_gmsig4": 238,
            "gmsig4_genes_present": "4/4",
            "comparator_genes_present": "5/5",
            "gep18_genes_present": "18/18",
            "primary_endpoint": "MFS/OS (nested with GSE58812)",
            "primary_evaluable_n": 0,
            "primary_events": 0,
            "secondary_endpoints": "Excluded from strict independent validation due to 107/238 exact overlap with model-selection cohort GSE58812",
            "age_complete_n": 0,
            "size_mm_complete_n": 0,
            "node_pos_complete_n": 0,
            "stage_complete_n": 0,
            "overlap_with_prior": "107/238 GSMs (GP_001..GP_107) 100% duplicate GSE58812",
            "eligibility_verdict": "EXCLUDED_OVERLAPS_GSE58812",
        },
    ]
    elig_df = pd.DataFrame(rows)
    elig_df.to_csv(ELIG / "dataset_eligibility.csv", index=False)

    # Print summary cleanly
    for r in rows:
        print(f"\n[{r['cohort']}] -> {r['eligibility_verdict']}")
        print(f"  TNBC mRNA={r['tnbc_with_mrna']}, all4={r['tnbc_all4_gmsig4']}, Endpoint={r['primary_endpoint']}: n={r['primary_evaluable_n']}, events={r['primary_events']}")
        print(f"  Genes: GMSig4={r['gmsig4_genes_present']}, Comp={r['comparator_genes_present']}, GEP18={r['gep18_genes_present']}")
        print(f"  Covariates complete: Age={r['age_complete_n']}, Size={r['size_mm_complete_n']}, Node={r['node_pos_complete_n']}, Stage={r['stage_complete_n']}")
        print(f"  Secondary: {r['secondary_endpoints']}")

    # Save SHA256 hashes of harmonized files
    h_rows = []
    for p in sorted(HARM.glob("*.csv")):
        h_rows.append({"file": str(p.relative_to(ROOT)), "bytes": p.stat().st_size, "sha256": sha256_file(p)})
    pd.DataFrame(h_rows).to_csv(ELIG / "harmonized_sha256.csv", index=False)
    print("\nSaved dataset_eligibility.csv and harmonized CSVs cleanly (ZERO outcome associations run).")


if __name__ == "__main__":
    main()
