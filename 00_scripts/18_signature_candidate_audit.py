import pandas as pd
import numpy as np
import os

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

res = os.path.join(project_root, "02_results")
ev = pd.read_csv(f"{res}/MultiCohort_TNBC_gutMGene_evidence_matrix.csv")
eq = pd.read_csv(f"{res}/MR_eQTLGen_BCAC_results.csv")
coloc = pd.read_csv(f"{res}/Coloc_eQTLGen_BCAC_results.csv")

# Summarize eQTLGen MR best hit per gene
eq_ivw = eq[eq["method"].isin(["Inverse variance weighted", "Wald ratio"])].copy()
eq_best = eq_ivw.sort_values("pval").groupby("exposure", as_index=False).first()
eq_best = eq_best.rename(columns={
    "exposure": "gene", "id.outcome": "eqtl_best_outcome",
    "OR": "eqtl_OR", "pval": "eqtl_pval", "fdr": "eqtl_fdr", "nsnp": "eqtl_nsnp"
})

coloc_best = coloc.sort_values("PP.H4", ascending=False).groupby("gene", as_index=False).first()[["gene", "outcome_id", "PP.H4"]]

m = ev.merge(eq_best[["gene", "eqtl_best_outcome", "eqtl_OR", "eqtl_pval", "eqtl_fdr", "eqtl_nsnp"]], on="gene", how="left")
m = m.merge(coloc_best, on="gene", how="left")

# Score candidates by multi-evidence support:
# 1) Prognosis: TCGA_OS_p < 0.05 + GSE58812_OS_p < 0.05 same direction
same_os_dir = ((m["TNBC_OS_HR"] > 1) & (m["GSE58812_OS_HR"] > 1)) | ((m["TNBC_OS_HR"] < 1) & (m["GSE58812_OS_HR"] < 1))
m["dual_OS_sig"] = same_os_dir & (m["TNBC_OS_p"] < 0.10) & (m["GSE58812_OS_p"] < 0.10)
m["any_OS_sig"]  = (m["TNBC_OS_p"] < 0.05) | (m["GSE58812_OS_p"] < 0.05)
m["pembro_sig"]  = (m["pembro_TNBC_p"] < 0.05) | (m["pembro_all_p"] < 0.05)
m["durva_sig"]   = (m["durva_TNBC_p"] < 0.05) | (m["durva_all_p"] < 0.05)
m["eqtl_sig"]    = m["eqtl_pval"] < 0.05
m["deg_sig"]     = (m["FDR_TNBC_vs_Normal"] < 0.05) & (m["FDR_TNBC_vs_NonTNBC"] < 0.05)

m["evidence_score"] = (
    2.0 * m["dual_OS_sig"].astype(int)
    + 1.0 * (m["TNBC_OS_p"] < 0.05).fillna(False).astype(int)
    + 1.0 * (m["GSE58812_OS_p"] < 0.05).fillna(False).astype(int)
    + 1.5 * m["pembro_sig"].fillna(False).astype(int)
    + 1.0 * m["durva_sig"].fillna(False).astype(int)
    + 1.5 * m["eqtl_sig"].fillna(False).astype(int)
    + 0.5 * m["anchored_GMrepo"].fillna(False).astype(int)
    + 0.5 * m["deg_sig"].fillna(False).astype(int)
)
m = m.sort_values(["evidence_score", "pembro_TNBC_p", "TCGA_OS_p" if "TCGA_OS_p" in m else "TNBC_OS_p"], ascending=[False, True, True])
m.to_csv(f"{res}/Integrated_gutMGene_MultiEvidence_Ranking.csv", index=False)

cols = ["gene", "anchored_GMrepo", "evidence_score", "TNBC_OS_HR", "TNBC_OS_p", "GSE58812_OS_HR", "GSE58812_OS_p",
        "pembro_TNBC_AUC", "pembro_TNBC_p", "durva_TNBC_AUC", "durva_TNBC_p", "eqtl_OR", "eqtl_pval", "PP.H4", "cor_CytotoxicImmune_rho"]
print("=== Top 35 Integrated Multi-Evidence gutMGene Targets ===")
print(m[cols].head(35).to_string(index=False))

print("\n=== All MR-Significant eQTLGen Genes + Prognosis/ICI overlap ===")
print(m[m["eqtl_sig"].fillna(False)][cols].to_string(index=False))
