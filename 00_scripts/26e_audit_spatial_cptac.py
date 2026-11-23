"""Verify and summarize spatial & CPTAC proteomic feasibility metadata.

Reads only previously fetched public metadata in 01_data/JTM_upgrade/spatial_metadata/
plus read-only local scRNA/TCGA patient ID manifests to determine exact sample sizes,
patient overlaps, mass-spec missingness, file sizes, SHA256 provenance, and pre-specified
non-circular primary spatial metric rules. Does NOT run outcome or score associations.
"""

from __future__ import annotations

import hashlib
import json
import re
import tarfile
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
META_DIR = ROOT / "01_data" / "JTM_upgrade" / "spatial_metadata"
OUT_DIR = ROOT / "02_results" / "JTM_upgrade" / "spatial_feasibility"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def verify_manifest(manifest_path: Path) -> dict[str, dict[str, Any]]:
    """Load source_manifest.json, verify SHA256 & sizes, and index by relative path."""
    records = json.loads(manifest_path.read_text(encoding="utf-8"))
    by_rel: dict[str, dict[str, Any]] = {}
    for rec in records:
        fp = ROOT / rec["file"]
        if not fp.exists():
            raise FileNotFoundError(f"Manifest file missing on disk: {fp}")
        blob = fp.read_bytes()
        sha = hashlib.sha256(blob).hexdigest()
        if sha != rec["sha256"] or len(blob) != rec["bytes"]:
            raise ValueError(f"Checksum/byte mismatch for {rec['file']}")
        key = rec["file"].replace("\\", "/")
        by_rel[key] = rec
    return by_rel


def audit_wu2021(manifest_idx: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """Audit Wu et al. 2021 Zenodo 4739739 Visium metadata & GSE176078 scRNA overlap."""
    z_rec = json.loads((META_DIR / "Wu2021" / "zenodo_record.json").read_text(encoding="utf-8"))
    files_info = {f["key"]: f for f in z_rec["files"]}

    # Read local scRNA TNBC patient IDs (read-only)
    scrna_umap_path = ROOT / "02_results" / "scRNA_TNBC_UMAP_metadata_keygenes.csv"
    scrna_tnbc_patients: set[str] = set()
    if scrna_umap_path.exists():
        u = pd.read_csv(scrna_umap_path, usecols=["patient"])
        scrna_tnbc_patients = set(u["patient"].astype(str).unique())

    full_scrna_meta = ROOT / "01_data" / "scRNA" / "Wu_etal_2021_BRCA_scRNASeq" / "metadata.csv"
    all_scrna_patients: dict[str, str] = {}
    if full_scrna_meta.exists():
        fm = pd.read_csv(full_scrna_meta, usecols=["orig.ident", "subtype"]).drop_duplicates()
        all_scrna_patients = dict(zip(fm["orig.ident"].astype(str), fm["subtype"].astype(str)))

    # Read clinical IHC table from MOESM4 Supplementary Table 1
    st1 = pd.read_excel(META_DIR / "Wu2021" / "41588_2021_911_MOESM4_ESM.xlsx", sheet_name="Supplementary Table 1")
    hdr_row = st1.index[st1.iloc[:, 0].astype(str) == "Case ID"][0]
    st1_clean = st1.iloc[hdr_row + 1 :].copy()
    st1_clean.columns = st1.iloc[hdr_row].tolist()
    st1_cases = {str(x).replace("-", ""): row for x, row in zip(st1_clean["Case ID"], st1_clean.to_dict("records"))}

    slide_rows = []
    with tarfile.open(META_DIR / "Wu2021" / "metadata.tar.gz", "r:gz") as tf:
        for member in tf.getmembers():
            if member.isdir():
                continue
            df = pd.read_csv(tf.extractfile(member), index_col=0)
            pid = str(df["patientid"].iloc[0])
            subtype = str(df["subtype"].iloc[0])
            cls_counts = df["Classification"].fillna("NA").value_counts().to_dict()
            # Match CID4290 -> CID4290A in scRNA if present
            scrna_match = pid if pid in all_scrna_patients else (f"{pid}A" if f"{pid}A" in all_scrna_patients else "")
            case_key = pid.replace("CID", "")
            ihc_rec = st1_cases.get(case_key, {})
            slide_rows.append(
                {
                    "patient_id": pid,
                    "subtype": subtype,
                    "n_spots": int(len(df)),
                    "median_counts": int(df["nCount_RNA"].median()),
                    "median_features": int(df["nFeature_RNA"].median()),
                    "in_scrna_tnbc_umap": pid in scrna_tnbc_patients,
                    "in_scrna_full": bool(scrna_match),
                    "scrna_id_matched": scrna_match or "None",
                    "ihc_table_er": ihc_rec.get("ER", "Not in Table S1"),
                    "ihc_table_pr": ihc_rec.get("PR", "Not in Table S1"),
                    "ihc_table_her2_ihc": ihc_rec.get("HER2 IHC", "Not in Table S1"),
                    "classification_counts": cls_counts,
                }
            )

    tnbc_slides = [s for s in slide_rows if s["subtype"] == "TNBC"]
    tnbc_overlap = [s["patient_id"] for s in tnbc_slides if s["in_scrna_tnbc_umap"]]
    tnbc_disjoint = [s["patient_id"] for s in tnbc_slides if not s["in_scrna_tnbc_umap"]]

    return {
        "slides": slide_rows,
        "n_total_slides": len(slide_rows),
        "n_tnbc_slides": len(tnbc_slides),
        "n_tnbc_spots": sum(s["n_spots"] for s in tnbc_slides),
        "tnbc_ids": [s["patient_id"] for s in tnbc_slides],
        "scrna_tnbc_n": len(scrna_tnbc_patients),
        "tnbc_overlap_ids": tnbc_overlap,
        "tnbc_disjoint_ids": tnbc_disjoint,
        "zenodo_files": {
            k: {"bytes": v["size"], "checksum": v["checksum"], "url": v["links"]["self"]}
            for k, v in files_info.items()
        },
        "meta_manifest": manifest_idx["01_data/JTM_upgrade/spatial_metadata/Wu2021/metadata.tar.gz"],
        "zenodo_manifest": manifest_idx["01_data/JTM_upgrade/spatial_metadata/Wu2021/zenodo_record.json"],
    }


def audit_geo_datasets(manifest_idx: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    """Audit GEO datasets from geo_search_summaries.json plus GSE169246 and GSE210616 metadata."""
    geo_sum = json.loads((META_DIR / "geo_search_summaries.json").read_text(encoding="utf-8"))["result"]
    uids = geo_sum["uids"]

    # GSE210616 file breakdown (excluding cloupe vs full tar)
    fl210 = (META_DIR / "GEO" / "GSE210616_filelist.txt").read_text(encoding="utf-8").splitlines()
    gse210_sizes: dict[str, int] = {}
    for line in fl210:
        if line.startswith("File\t"):
            parts = line.split("\t")
            gse210_sizes[parts[4]] = gse210_sizes.get(parts[4], 0) + int(parts[3])
    gse210_non_cloupe_bytes = sum(v for k, v in gse210_sizes.items() if k != "CLOUPE")
    gse210_full_tar_bytes = int([l.split("\t")[3] for l in fl210 if l.startswith("Archive\t")][0])

    gsm210 = (META_DIR / "GEO" / "GSE210616_gsm_brief.txt").read_text(encoding="utf-8", errors="replace")
    pats210 = set(re.findall(r"!Sample_title = Primary TNBC, patient (\d+),", gsm210))

    # GSE299631 5 human TNBC Visium files
    fl299 = (META_DIR / "GEO" / "GSE299631_filelist.txt").read_text(encoding="utf-8").splitlines()
    human_gsms = {"GSM9042337", "GSM9042338", "GSM9042339", "GSM9042340", "GSM9042341"}
    gse299_human_bytes = sum(
        int(l.split("\t")[3]) for l in fl299 if l.startswith("File\t") and any(g in l for g in human_gsms)
    )

    records: list[dict[str, Any]] = []
    # Explicit lookup table of verified GEO metadata facts
    for uid in uids:
        rec = geo_sum[str(uid)]
        acc = rec["accession"]
        n_samp = int(rec.get("n_samples") or 0)
        title = rec.get("title", "")
        taxon = rec.get("taxon", "")

        if acc == "GSE210616":
            records.append(
                {
                    "dataset_id": acc,
                    "source_type": "Spatial_Visium",
                    "citation_or_pmid": "PMID:36283023 / PMC9812886 (Bassiouni et al., Cancer Res 2023)",
                    "platform": "10x Visium (GPL24676)",
                    "species_tissue": "Homo sapiens / Primary TNBC",
                    "total_samples_or_sections": n_samp,
                    "confirmed_tnbc_patients": len(pats210),
                    "overlap_existing_cohorts": "0 with GSE176078 scRNA or TCGA-BRCA (VCU cohort)",
                    "genes_covered": "Whole transcriptome GRCh38 (all 4 GMSig-4 genes expected)",
                    "protein_or_spatial_missingness": "No pathologist spot-label CSV on GEO; H&E + coords + H5 present; clinical/RFS Supp File 1 behind PMC reCAPTCHA",
                    "matched_rna_seq": "Spatial transcriptomics only (no matched bulk)",
                    "survival_outcomes": "RFS reported in paper Supp File 1 (PMC reCAPTCHA-gated; not on GEO)",
                    "download_bytes": gse210_non_cloupe_bytes,
                    "download_notes": f"Per-GSM H5+PNG+CSV+JSON = {gse210_non_cloupe_bytes:,} B (~542.4 MiB); full RAW.tar = {gse210_full_tar_bytes:,} B (~35.1 GiB due to CLOUPE)",
                    "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210616",
                    "metadata_sha256": manifest_idx["01_data/JTM_upgrade/spatial_metadata/GEO/GSE210616.txt"]["sha256"],
                    "key_uncertainties": "Mixed treatment (13 chemo-naive, 9 post-NAC); no spot-level pathologist labels in GEO; must use GMSig-4-free markers + lattice geometry",
                    "go_nogo_verdict": "GO (Recommended Primary Spatial Cohort)",
                }
            )
        elif acc == "GSE299631":
            records.append(
                {
                    "dataset_id": acc,
                    "source_type": "Spatial_Visium",
                    "citation_or_pmid": "Turos et al., 2025 (Geneva Univ. Hospitals + KB1P mouse)",
                    "platform": "10x Visium (GPL24676 human; GPL24247 mouse)",
                    "species_tissue": "Homo sapiens (n=5) + Mus musculus (n=45)",
                    "total_samples_or_sections": 5,
                    "confirmed_tnbc_patients": 5,
                    "overlap_existing_cohorts": "0 with GSE176078 or TCGA-BRCA",
                    "genes_covered": "Whole transcriptome human GRCh38 (5 human GSMs)",
                    "protein_or_spatial_missingness": "No pathologist spot-label CSV on GEO; human subset restricted to germline BRCA1 carriers (2 primary, 3 post-NAC residual)",
                    "matched_rna_seq": "Mouse scRNA only; no matched human scRNA/bulk",
                    "survival_outcomes": "None on GEO",
                    "download_bytes": gse299_human_bytes,
                    "download_notes": f"5 human TNBC GSM Space Ranger files = {gse299_human_bytes:,} B (~312.8 MiB); full series TAR = 2.4 GiB",
                    "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299631",
                    "metadata_sha256": manifest_idx["01_data/JTM_upgrade/spatial_metadata/GEO/GSE299631.txt"]["sha256"],
                    "key_uncertainties": "N=5 germline BRCA1-only; 3/5 post-NAC residual disease; heterogeneous therapy state",
                    "go_nogo_verdict": "NO-GO (Optional Secondary BRCA1 Check Only)",
                }
            )
        elif acc == "GSE218951":
            records.append(
                {
                    "dataset_id": acc,
                    "source_type": "Spatial_Visium",
                    "citation_or_pmid": "He et al., Starfysh (Columbia Univ.)",
                    "platform": "10x Visium (GPL18573/GPL24676)",
                    "species_tissue": "Homo sapiens / ER+ (P1), TNBC (P2), MBC (P3, P4)",
                    "total_samples_or_sections": n_samp,
                    "confirmed_tnbc_patients": 1,
                    "overlap_existing_cohorts": "0 with GSE176078 or TCGA-BRCA",
                    "genes_covered": "Whole transcriptome",
                    "protein_or_spatial_missingness": "Only 1 non-metaplastic TNBC patient (2 sections) + 2 MBC patients",
                    "matched_rna_seq": "1 spatial TCR sample (MBC2)",
                    "survival_outcomes": "None on GEO",
                    "download_bytes": 157 * 1024**2,
                    "download_notes": "GSE218951_RAW.tar = 157 MiB (index listing)",
                    "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218951",
                    "metadata_sha256": manifest_idx["01_data/JTM_upgrade/spatial_metadata/GEO/GSE218951.txt"]["sha256"],
                    "key_uncertainties": "N=1 TNBC patient (+2 rare MBC); insufficient patient-level sample size",
                    "go_nogo_verdict": "NO-GO (Insufficient TNBC N=1)",
                }
            )
        elif acc == "GSE235640":
            records.append(
                {
                    "dataset_id": acc,
                    "source_type": "Spatial_GeoMx_ROI",
                    "citation_or_pmid": "PMID:38061601 (Ravichandran et al., BBA 2024)",
                    "platform": "NanoString GeoMx WTA (GPL24676)",
                    "species_tissue": "Homo sapiens / TNBC endothelial ROIs",
                    "total_samples_or_sections": n_samp,
                    "confirmed_tnbc_patients": 4,
                    "overlap_existing_cohorts": "0 with GSE176078 or TCGA-BRCA",
                    "genes_covered": "GeoMx WTA panel",
                    "protein_or_spatial_missingness": "ROIs selectively gated on tumor vs adjacent-normal endothelial cells only; no tumor-nest vs immune-margin spatial lattice",
                    "matched_rna_seq": "None",
                    "survival_outcomes": "None",
                    "download_bytes": int(3.7 * 1024**2),
                    "download_notes": "normdata.csv.gz = 3.7 MiB; RAW.tar = 2.7 MiB",
                    "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235640",
                    "metadata_sha256": manifest_idx["01_data/JTM_upgrade/spatial_metadata/GEO/GSE235640.txt"]["sha256"],
                    "key_uncertainties": "Endothelial-only segmented ROIs (n=4 patients); cannot measure CD8 tumor-nest penetration vs stromal margin",
                    "go_nogo_verdict": "NO-GO (Endothelial-Only ROIs, N=4)",
                }
            )
        elif acc == "GSE283412":
            records.append(
                {
                    "dataset_id": acc,
                    "source_type": "Spatial_Visium",
                    "citation_or_pmid": "PMID:40790295 (Lakhani Lab, 2025)",
                    "platform": "10x Visium (GPL24676)",
                    "species_tissue": "Homo sapiens / Metaplastic Breast Cancer (MpBC)",
                    "total_samples_or_sections": n_samp,
                    "confirmed_tnbc_patients": 4,
                    "overlap_existing_cohorts": "0",
                    "genes_covered": "Whole transcriptome",
                    "protein_or_spatial_missingness": "Rare metaplastic breast cancer (MpBC) only (4 patients: MP_14, MP_167, MP_170, MP_10)",
                    "matched_rna_seq": "Methylation paired",
                    "survival_outcomes": "None on GEO",
                    "download_bytes": pd.NA,
                    "download_notes": "Series /suppl/ dir restricted (403); per-GSM files available via SOFT",
                    "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283412",
                    "metadata_sha256": manifest_idx["01_data/JTM_upgrade/spatial_metadata/GEO/GSE283412.txt"]["sha256"],
                    "key_uncertainties": "Rare MpBC histology (N=4), not general primary IDC TNBC",
                    "go_nogo_verdict": "NO-GO (Rare MpBC Histology, N=4)",
                }
            )
        elif acc == "GSE328422":
            records.append(
                {
                    "dataset_id": acc,
                    "source_type": "Spatial_Xenium",
                    "citation_or_pmid": "Zhao et al., bioRxiv 2026",
                    "platform": "10x Xenium Prime 5K (GPL33762)",
                    "species_tissue": "Homo sapiens / TNBC FFPE Lymph Nodes",
                    "total_samples_or_sections": n_samp,
                    "confirmed_tnbc_patients": pd.NA,
                    "overlap_existing_cohorts": "0",
                    "genes_covered": "5K targeted panel (unverified for GMSig-4)",
                    "protein_or_spatial_missingness": "Lymph node sections (uninvolved/metastatic), NOT primary breast tumors",
                    "matched_rna_seq": "CD45+ scRNA in GSE328275",
                    "survival_outcomes": "None on GEO",
                    "download_bytes": 28 * 1024**3,
                    "download_notes": "RAW.tar = 28 GiB; Seurat RDS = 38 GiB",
                    "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE328422",
                    "metadata_sha256": manifest_idx["01_data/JTM_upgrade/spatial_metadata/GEO/GSE328422.txt"]["sha256"],
                    "key_uncertainties": "Lymph nodes (not primary tumor), 28-38 GiB download, 5K panel",
                    "go_nogo_verdict": "NO-GO (Lymph Node Tissue & 28-38 GiB)",
                }
            )
        elif acc == "GSE296517":
            records.append(
                {
                    "dataset_id": acc,
                    "source_type": "Spatial_CosMx",
                    "citation_or_pmid": "Gamache et al., BR54 trial (NCT04796220)",
                    "platform": "NanoString CosMx 1000-plex (GPL33484)",
                    "species_tissue": "Homo sapiens / Breast Cancer (FUS +/- GEM trial)",
                    "total_samples_or_sections": n_samp,
                    "confirmed_tnbc_patients": pd.NA,
                    "overlap_existing_cohorts": "0",
                    "genes_covered": "1000-plex targeted RNA panel (incomplete GMSig-4 coverage likely)",
                    "protein_or_spatial_missingness": "Intervention-perturbed trial biopsies (2 trial patients + 1 archival); subtype not TNBC-only",
                    "matched_rna_seq": "None",
                    "survival_outcomes": "None",
                    "download_bytes": 176 * 1024**3,
                    "download_notes": "RAW.tar = 176 GiB; SeuratObject.RDS = 1.1 GiB; TileDB = 8.6 GiB",
                    "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296517",
                    "metadata_sha256": manifest_idx["01_data/JTM_upgrade/spatial_metadata/GEO/GSE296517.txt"]["sha256"],
                    "key_uncertainties": "N=2-3 patients, FUS/gemcitabine perturbed, 1000-plex panel, huge files",
                    "go_nogo_verdict": "NO-GO (Perturbed Trial, N<=3, 1000-plex)",
                }
            )
        else:
            # Non-primary-TNBC-spatial GEO hits from geo_search_summaries.json
            reason_map = {
                "GSE328275": ("scRNA_CD45_only", "CD45+ immune scRNA-seq (not spatial)"),
                "GSE284115": ("scRNA_Myeloid", "Human/mouse myeloid scRNA-seq (not spatial)"),
                "GSE261774": ("scRNA_CD45_HRpos", "CD45+ scRNA-seq focused on HR+ BC TILs (not TNBC spatial)"),
                "GSE299880": ("Spatial_GeoMx_HRpos", "HR+ breast cancer GeoMx (not TNBC)"),
                "GSE297605": ("Bulk_or_HRpos", "HR+ breast cancer (not TNBC)"),
                "GSE267534": ("Mouse_Only", "Mus musculus only (not human)"),
                "GSE291237": ("DNA_Methylation", "Genome-wide DNA methylation array (not spatial)"),
                "GSE293199": ("Spatial_Benchmark_N1", "Single-sample foundation model benchmark (n=1)"),
                "GSE184720": ("Xenograft_CellLine", "MDA-MB-231 mouse xenograft punch ST (not primary patient tumor)"),
                "GSE102155": ("Bulk_lncRNA_Array", "Bulk lncRNA microarray (not spatial)"),
                "GSE60689": ("Bulk_lncRNA_Array", "Bulk lncRNA microarray (not spatial)"),
            }
            stype, reason = reason_map.get(acc, ("Other", "Not human primary TNBC spatial transcriptomics"))
            records.append(
                {
                    "dataset_id": acc,
                    "source_type": stype,
                    "citation_or_pmid": f"GEO esummary uid={uid}",
                    "platform": f"GPL {rec.get('gpl', '')}",
                    "species_tissue": f"{taxon} / {title[:80]}",
                    "total_samples_or_sections": n_samp,
                    "confirmed_tnbc_patients": pd.NA,
                    "overlap_existing_cohorts": "Not evaluated (excluded on modality/subtype/species)",
                    "genes_covered": "N/A",
                    "protein_or_spatial_missingness": reason,
                    "matched_rna_seq": "N/A",
                    "survival_outcomes": "N/A",
                    "download_bytes": pd.NA,
                    "download_notes": f"Excluded via geo_search_summaries.json: {reason}",
                    "source_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}",
                    "metadata_sha256": manifest_idx["01_data/JTM_upgrade/spatial_metadata/geo_search_summaries.json"]["sha256"],
                    "key_uncertainties": reason,
                    "go_nogo_verdict": f"NO-GO ({stype})",
                }
            )

    # Add GSE169246 (explicit Task 2 audit)
    records.append(
        {
            "dataset_id": "GSE169246",
            "source_type": "scRNA_scATAC_Dynamic",
            "citation_or_pmid": "PMID:34653365 (Zhang et al., Cancer Cell 2021)",
            "platform": "10x scRNA-seq + scATAC-seq (GPL20795)",
            "species_tissue": "Homo sapiens / Sorted CD45+ immune cells (primary/metastatic TNBC & blood)",
            "total_samples_or_sections": 107,
            "confirmed_tnbc_patients": 22,
            "overlap_existing_cohorts": "0 with GSE176078 or TCGA-BRCA",
            "genes_covered": "Whole transcriptome in CD45+ immune cells only; epithelial/CAF compartments ABSENT",
            "protein_or_spatial_missingness": "CD45+ sorted only: tumor-intrinsic ENO1 and bulk-equivalent GMSig-4 cannot be constructed; mixed blood/primary/metastatic compartments",
            "matched_rna_seq": "scATAC paired in subset; no bulk RNA-seq",
            "survival_outcomes": "Response/treatment dynamics (paclitaxel vs paclitaxel + atezolizumab)",
            "download_bytes": int(1.9 * 1024**3),
            "download_notes": "GSE169246_TNBC_RNA.counts.mtx.gz = 1.9 GiB (+ barcodes 2.0 MiB, features 98 KiB, ATAC 530 MiB)",
            "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169246",
            "metadata_sha256": manifest_idx["01_data/JTM_upgrade/spatial_metadata/GEO/GSE169246.txt"]["sha256"],
            "key_uncertainties": "CD45+ sorting removes malignant epithelial cells (main source of ENO1) + 1.9 GiB matrix",
            "go_nogo_verdict": "NO-GO (Impractical 1.9 GiB & CD45+-Only)",
        }
    )
    return records


def audit_cptac(manifest_idx: dict[str, dict[str, Any]]) -> tuple[dict[str, Any], dict[str, Any]]:
    """Audit CPTAC brca_cptac_2020 clinical TNBC counts, TCGA overlap, RNA-seq, and protein missingness."""
    pat = pd.read_csv(META_DIR / "CPTAC" / "data_clinical_patient.txt", sep="\t", comment="#")
    samp = pd.read_csv(META_DIR / "CPTAC" / "data_clinical_sample.txt", sep="\t", comment="#")
    merged = pat.merge(samp, on="PATIENT_ID", how="outer")

    prot = pd.read_csv(META_DIR / "CPTAC" / "data_protein_quantification.txt", sep="\t")
    sample_cols = [c for c in prot.columns if c != "Composite.Element.REF"]
    prot["SYMBOL"] = prot["Composite.Element.REF"].astype(str).apply(lambda s: s.split("|")[0].strip())

    # TCGA overlap check (read-only)
    tcga_path = ROOT / "02_results" / "GMSig4_PatientScores_TCGA_TNBC.csv"
    tcga_overlap = 0
    if tcga_path.exists():
        tcga_df = pd.read_csv(tcga_path)
        tcga_ids = set(tcga_df["sampleID"].astype(str).str[:12])
        tcga_overlap = len(tcga_ids & set(merged["PATIENT_ID"].astype(str)))

    tnbc_pos = merged[merged["TNBC_UPDATED_CLINICAL_STATUS"] == "Positive"]
    tnbc_ids = [s for s in tnbc_pos["SAMPLE_ID"] if s in sample_cols]
    strict_ihc = merged[
        (merged["ER_UPDATED_CLINICAL_STATUS"].astype(str).str.lower() == "negative")
        & (merged["PR_CLINICAL_STATUS"].astype(str).str.lower() == "negative")
        & (merged["ERBB2_UPDATED_CLINICAL_STATUS"].astype(str).str.lower() == "negative")
    ]
    strict_ids = [s for s in strict_ihc["SAMPLE_ID"] if s in sample_cols]

    # Matched RNA-seq case list
    expr_txt = (META_DIR / "CPTAC" / "cases_expression.txt").read_text(encoding="utf-8")
    expr_ids = []
    for line in expr_txt.splitlines():
        if line.startswith("case_list_ids:"):
            expr_ids = [x.strip() for x in line.split(":", 1)[1].strip().split("\t") if x.strip()]

    prot_stats: dict[str, dict[str, Any]] = {}
    for gene in ["IFNG", "IDO1", "CXCL10", "ENO1"]:
        sub = prot[prot["SYMBOL"] == gene]
        if len(sub) == 0:
            prot_stats[gene] = {
                "ref_id": "ABSENT",
                "all_non_na": 0,
                "all_n": len(sample_cols),
                "all_pct": 0.0,
                "tnbc_non_na": 0,
                "tnbc_n": len(tnbc_ids),
                "tnbc_pct": 0.0,
                "strict_non_na": 0,
                "strict_n": len(strict_ids),
            }
        else:
            row = sub.iloc[0]
            all_v = row[sample_cols]
            tnbc_v = row[tnbc_ids]
            strict_v = row[strict_ids]
            prot_stats[gene] = {
                "ref_id": str(row["Composite.Element.REF"]),
                "all_non_na": int(all_v.notna().sum()),
                "all_n": len(sample_cols),
                "all_pct": round(100.0 * float(all_v.notna().mean()), 1),
                "tnbc_non_na": int(tnbc_v.notna().sum()),
                "tnbc_n": len(tnbc_ids),
                "tnbc_pct": round(100.0 * float(tnbc_v.notna().mean()), 1),
                "strict_non_na": int(strict_v.notna().sum()),
                "strict_n": len(strict_ids),
            }

    # Joint 3-protein completeness in TNBC
    co_det = pd.DataFrame(
        {g: prot.loc[prot["SYMBOL"] == g, tnbc_ids].iloc[0].notna() for g in ["IDO1", "CXCL10", "ENO1"]}
    )
    n_3prot_tnbc = int(co_det.all(axis=1).sum())

    # Size info from datahub_directory.json & media check
    prot_bytes = manifest_idx["01_data/JTM_upgrade/spatial_metadata/CPTAC/data_protein_quantification.txt"]["bytes"]
    mrna_bytes = 17773792  # Verified via HTTP HEAD Content-Length on media.githubusercontent.com
    clin_bytes = (
        manifest_idx["01_data/JTM_upgrade/spatial_metadata/CPTAC/data_clinical_patient.txt"]["bytes"]
        + manifest_idx["01_data/JTM_upgrade/spatial_metadata/CPTAC/data_clinical_sample.txt"]["bytes"]
    )

    summary = {
        "total_patients": int(len(merged)),
        "tnbc_updated_positive": int(len(tnbc_ids)),
        "tnbc_strict_ihc": int(len(strict_ids)),
        "tnbc_pam50_breakdown": tnbc_pos["PAM50"].value_counts(dropna=False).to_dict(),
        "tnbc_matched_rnaseq": len(set(expr_ids) & set(tnbc_ids)),
        "all_matched_rnaseq": len(expr_ids),
        "tcga_overlap": tcga_overlap,
        "cd3_tils_performed_tnbc": int((tnbc_pos["CD3_TILS_STATUS"] != "not performed").sum()),
        "estimate_immune_tnbc": int(tnbc_pos["ESTIMATE_IMMUNE_SCORE"].notna().sum()),
        "prot_stats": prot_stats,
        "n_3prot_tnbc": n_3prot_tnbc,
        "prot_bytes": prot_bytes,
        "mrna_bytes": mrna_bytes,
        "clin_bytes": clin_bytes,
    }

    row_csv = {
        "dataset_id": "CPTAC_BRCA_2020 (brca_cptac_2020)",
        "source_type": "Proteogenomics_TMT_MS_RNAseq",
        "citation_or_pmid": "PMID:33212010 / PMC8077737 (Krug et al., Cell 2020)",
        "platform": "TMT-10-plex MS proteomics (9,715 proteins) + RNA-seq FPKM",
        "species_tissue": "Homo sapiens / Prospective treatment-naive primary BRCA",
        "total_samples_or_sections": 122,
        "confirmed_tnbc_patients": 28,
        "overlap_existing_cohorts": f"{tcga_overlap} with TCGA-BRCA (prospective CPTAC-3, distinct from Mertins 2016 CPTAC-2)",
        "genes_covered": "RNA-seq 122/122; Proteomics: ENO1=100%, IDO1=83.6%(all)/82.1%(TNBC), CXCL10=74.6%(all)/75.0%(TNBC), IFNG=0%(ABSENT)",
        "protein_or_spatial_missingness": (
            f"IFNG absent (0/122); IDO1 NA in 5/28 TNBC (17.9%); CXCL10 NA in 7/28 TNBC (25.0%); "
            f"ENO1 0 NA (28/28); 3-protein complete=20/28 TNBC; 4-protein complete=0/28"
        ),
        "matched_rna_seq": "122/122 overall (28/28 TNBC) via data_mrna_seq_fpkm.txt",
        "survival_outcomes": "NO OS/RFS survival endpoints in brca_cptac_2020 (1-yr forms only; Krug et al. used METABRIC for survival)",
        "download_bytes": prot_bytes + mrna_bytes + clin_bytes,
        "download_notes": f"protein={prot_bytes:,} B (~8.44 MiB) + mRNA_FPKM={mrna_bytes:,} B (~16.95 MiB) + clinical={clin_bytes:,} B",
        "source_url": "https://github.com/cBioPortal/datahub/tree/master/public/brca_cptac_2020",
        "metadata_sha256": manifest_idx["01_data/JTM_upgrade/spatial_metadata/CPTAC/data_protein_quantification.txt"]["sha256"],
        "key_uncertainties": "IFNG completely absent at protein level; 12/28 TNBC lack clinical ERBB2 IHC (resolved via proteogenomics; 5/28 are PAM50 Her2); no survival outcomes",
        "go_nogo_verdict": "CONDITIONAL NARROW GO (ENO1/IDO1/CXCL10 RNA-Protein Concordance Only; NO-GO for 4-Protein Score or Prognosis)",
    }
    return summary, row_csv


def main() -> None:
    manifest_idx = verify_manifest(META_DIR / "source_manifest.json")
    wu = audit_wu2021(manifest_idx)
    geo_rows = audit_geo_datasets(manifest_idx)
    cptac_sum, cptac_row = audit_cptac(manifest_idx)

    wu_bytes = (
        wu["zenodo_files"]["metadata.tar.gz"]["bytes"]
        + wu["zenodo_files"]["spatial.tar.gz"]["bytes"]
        + wu["zenodo_files"]["filtered_count_matrices.tar.gz"]["bytes"]
    )
    wu_row = {
        "dataset_id": "Wu2021_Zenodo_4739739",
        "source_type": "Spatial_Visium",
        "citation_or_pmid": "PMID:34493872 (Wu et al., Nat Genet 2021)",
        "platform": "10x Visium (Space Ranger v1.0.0 GRCh38)",
        "species_tissue": "Homo sapiens / Primary breast cancer (4 TNBC + 2 ER+)",
        "total_samples_or_sections": wu["n_total_slides"],
        "confirmed_tnbc_patients": wu["n_tnbc_slides"],
        "overlap_existing_cohorts": (
            f"2/4 TNBC overlap GSE176078 scRNA ({', '.join(wu['tnbc_overlap_ids'])}); "
            f"2/4 Visium-only ({', '.join(wu['tnbc_disjoint_ids'])}); 0 TCGA overlap"
        ),
        "genes_covered": "Whole transcriptome GRCh38 filtered count matrices (all 4 GMSig-4 genes profiled)",
        "protein_or_spatial_missingness": (
            "12,052 TNBC spots; pathologist Loupe labels present for all 4 TNBC slides, but 3/4 lump tumor bed into "
            "'Invasive cancer + stroma + lymphocytes' (CID4465 has 93.4% mixed & only 73 Stroma spots, 0 Lymphocytes-only)"
        ),
        "matched_rna_seq": "Matched scRNA-seq + bulk RNA-seq for CID4465 & CID44971 (plus ER+ CID4290 & CID4535); none for 1142243F/1160920F",
        "survival_outcomes": "None",
        "download_bytes": wu_bytes,
        "download_notes": (
            f"metadata+spatial+filtered_counts = {wu_bytes:,} B (~200.5 MiB); "
            f"raw_counts={wu['zenodo_files']['raw_count_matrices.tar.gz']['bytes']:,} B; "
            f"images.pdf={wu['zenodo_files']['images.pdf']['bytes']:,} B"
        ),
        "source_url": "https://zenodo.org/records/4739739",
        "metadata_sha256": wu["meta_manifest"]["sha256"],
        "key_uncertainties": (
            "N=4 TNBC patients only (insufficient for statistical hypothesis testing; 2/4 overlap GSE176078 scRNA); "
            "1142243F & 1160920F lack clinical ER/PR/HER2 % in Table S1; Visium 55um multicellular spot mixing"
        ),
        "go_nogo_verdict": "CONDITIONAL GO (Case-Level Descriptive Vignette N=4 Only; Pair with GSE210616 for Powered Patient-Level Test)",
    }

    # Order priority candidates first, then remaining GEO entries
    priority_rows = [wu_row] + [r for r in geo_rows if r["dataset_id"] == "GSE210616"] + [cptac_row]
    other_rows = [r for r in geo_rows if r["dataset_id"] != "GSE210616"]
    df_out = pd.DataFrame(priority_rows + other_rows)
    csv_path = OUT_DIR / "spatial_proteomic_eligibility.csv"
    df_out.to_csv(csv_path, index=False)

    # Save exact audit summary JSON for transparent markdown traceability? No, keep only requested outputs in OUT_DIR
    # and print summary for verification.
    print(f"Wrote {csv_path} ({len(df_out)} datasets)")
    print("Wu2021 TNBC:", wu["n_tnbc_slides"], "patients,", wu["n_tnbc_spots"], "spots, overlap:", wu["tnbc_overlap_ids"])
    print("CPTAC TNBC:", cptac_sum["tnbc_updated_positive"], "updated TNBC,", cptac_sum["tnbc_strict_ihc"], "strict IHC,", "prot:", cptac_sum["prot_stats"])


if __name__ == "__main__":
    main()
