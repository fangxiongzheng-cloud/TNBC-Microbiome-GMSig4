import os
import pandas as pd

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

r = os.path.join(project_root, "02_results")
for fn in ["MR_MiBioGen_GCST010098_results.csv", "MR_eQTLGen_GCST010098_results.csv"]:
    df = pd.read_csv(f"{r}/{fn}")
    df = df[df["method"].isin(["Inverse variance weighted", "Wald ratio"]) & (df["pval"] < 0.05)].sort_values("pval")
    print(f"\n=== {fn} (p<0.05 n={len(df)}) ===")
    print(df[["exposure", "method", "nsnp", "OR", "OR_lci", "OR_uci", "pval", "fdr"]].to_string(index=False))

eq_rep = pd.read_csv(f"{r}/MR_eQTLGen_5GWAS_CrossSubtype_Long.csv")
print("\n=== Key eQTLGen targets across 5 GWAS ===")
genes = ["FFAR2", "PEMT", "CCR6", "MUC1", "TCF7L2", "NLRP6", "FOXP3", "MAGI3", "ICAM1", "ST6GAL1", "CXCR6", "CXCR2"]
print(eq_rep[eq_rep["exposure"].isin(genes)].sort_values(["exposure", "pval"])[["exposure", "outcome_id", "nsnp", "OR", "pval", "fdr"]].to_string(index=False))
