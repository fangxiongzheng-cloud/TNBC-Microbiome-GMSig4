import os
import pandas as pd

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

res_dir = os.path.join(project_root, "02_results")
fig_dir = os.path.join(project_root, "03_figures")

print("=== 1. HGMT MICROBIOME COHORTS ===")
m8 = pd.read_csv(f"{res_dir}/HGMT_all8_cohorts_metadata_summary.csv")
print(f"Total cohorts: {len(m8)}, Total bacteria runs: {m8['n_runs_bacteria'].sum()}, BC: {m8['n_BC'].sum()}, Benign: {m8['n_BreastDiseases'].sum()}, Health: {m8['n_Health'].sum()}")
grad = pd.read_csv(f"{res_dir}/HGMT_benign_to_malignant_gradient.csv")
print("PRJNA658160 top gradient genera (|rho|>0.12):")
print(grad[(grad['project_id']=='PRJNA658160') & (grad['spearman_rho'].abs()>0.12)][['feature','spearman_rho','spearman_p','mean_Health','mean_BreastDiseases','mean_BC']].head(8).to_string(index=False))

print("\n=== 2. 6-GWAS MR & COLOCALIZATION & MEDIATION ===")
mib5 = pd.read_csv(f"{res_dir}/MR_MiBioGen_5GWAS_CrossSubtype_Long.csv")
mib98 = pd.read_csv(f"{res_dir}/MR_MiBioGen_GCST010098_results.csv")
print(f"MiBioGen tested traits per GWAS ~{mib5.groupby('outcome_id').size().to_dict()}, GCST010098 p<0.05: {len(mib98[(mib98['method']=='Inverse variance weighted') & (mib98['pval']<0.05)])}")
eq5 = pd.read_csv(f"{res_dir}/MR_eQTLGen_5GWAS_CrossSubtype_Long.csv")
eq98 = pd.read_csv(f"{res_dir}/MR_eQTLGen_GCST010098_results.csv")
print(f"eQTLGen tested genes per GWAS: {eq5.groupby('outcome_id').size().to_dict()}, GCST010098 p<0.05: {len(eq98[eq98['method'].isin(['Inverse variance weighted','Wald ratio']) & (eq98['pval']<0.05)])}")
coloc = pd.read_csv(f"{res_dir}/Coloc_GTExBreast_TNBC_results.csv")
print(f"GTEx Breast coloc genes: {len(coloc)}; top cond_H4 >= 0.35:")
print(coloc[coloc['cond_H4'] >= 0.35][['gene','n_shared_snps','min_gtex_p','min_tnbc_p','PP.H4','cond_H4','lead_snp']].to_string(index=False))

med_path = f"{res_dir}/MR_Mediation_2Step_Causal_Pathways.csv"
if os.path.exists(med_path):
    med = pd.read_csv(med_path)
    print(f"2-Step Mediation rows: {len(med)}")
    print(med.head(15).to_string(index=False))
else:
    print("MR_Mediation_2Step_Causal_Pathways.csv not yet present.")

iv_m = pd.read_csv(f"{res_dir}/MiBioGen_clumped_IVs_1e5.csv")
gx_e = pd.read_csv(os.path.join(project_root, "01_data", "eQTL", "GTEx_V8_Breast_target_genes_500kb.csv"))
m_ov = iv_m.merge(gx_e, left_on="rsID", right_on="rsid")
print(f"Local GTEx Breast x MiBioGen IV overlap rows: {len(m_ov)}")
if len(m_ov) > 0:
    print(m_ov[["bac","rsID","gene_id","beta_x","P","beta_y","pvalue"]].sort_values("pvalue").head(15).to_string(index=False))

print("\n=== 3. GMSIG-4 PROGNOSIS & NOMOGRAM ===")
prog = pd.read_csv(f"{res_dir}/GMSig4_Prognostic_Validation_Summary.csv")
print(prog.to_string(index=False))
nomo = pd.read_csv(f"{res_dir}/GMSig4_Nomogram_Comparison.csv")
print(nomo.to_string(index=False))

print("\n=== 4. IMMUNOTHERAPY pCR & scRNA ===")
ici = pd.read_csv(f"{res_dir}/GMSig4_Immunotherapy_pCR_Validation_Summary.csv")
print(ici[['cohort','n','pCR_n','GMSig4_AUC','GMSig4_wilcox_p','pCR_rate_LowRisk','pCR_rate_HighRisk','Fisher_p','UniBeta_AUC','CD274_AUC','PDCD1_AUC']].to_string(index=False))
umap = pd.read_csv(f"{res_dir}/scRNA_TNBC_UMAP_metadata_keygenes.csv")
print(f"scRNA TNBC total cells: {len(umap)}, patients: {umap['patient'].nunique()}")

print("\n=== 5. VINA DOCKING ===")
dock = pd.read_csv(f"{res_dir}/Docking_Vina_Affinities_Matrix.csv")
print(f"Total docked complexes: {len(dock)}, gutMGene2 annotated pairs: {(dock['gutMGene_evidence'].fillna('') != '').sum()}")
print(dock[dock['gutMGene_evidence'].fillna('') != ''][['gene','metabolite','affinity_kcal_mol','gutMGene_evidence']].to_string(index=False))

print("\n=== 6. FIGURES CHECK ===")
for f in sorted(os.listdir(fig_dir)):
    sz = os.path.getsize(f"{fig_dir}/{f}")
    print(f"  {f}: {sz/1024:.1f} KB")
