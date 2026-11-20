import csv
import os

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

res_dir = os.path.join(project_root, "02_results")
gut_dir = os.path.join(project_root, "01_data", "gutMGene")

# Load GMrepo cross-cohort differential genera
cons = list(csv.DictReader(open(os.path.join(res_dir, "GMrepo_genus_cross_cohort_consistency.csv"), encoding="utf-8")))
diff_by = list(csv.DictReader(open(os.path.join(res_dir, "GMrepo_genus_differential_by_cohort.csv"), encoding="utf-8")))

# Filter candidate genera: either FDR<0.05 in >=1 cohort OR p<0.05 in >=2 cohorts
cand_genera = {}
for row in cons:
    if int(row["n_wilcox_fdr05"]) >= 1 or int(row["n_linda_fdr05"]) >= 1 or int(row["n_wilcox_p05"]) >= 2 or int(row["n_linda_p05"]) >= 2:
        cand_genera[row["genus"]] = row

print(f"Candidate GMrepo differential genera: {len(cand_genera)}")

# Load gutMGene v2.0 human tables
mm = list(csv.DictReader(open(os.path.join(gut_dir, "gutMGene2_human_microbe_metabolite.csv"), encoding="utf-8-sig")))
mg = list(csv.DictReader(open(os.path.join(gut_dir, "gutMGene2_human_microbe_gene.csv"), encoding="utf-8-sig")))
tg = list(csv.DictReader(open(os.path.join(gut_dir, "gutMGene2_human_metabolite_gene.csv"), encoding="utf-8-sig")))

def match_genus(microbe_name, genus):
    m = microbe_name.lower().strip()
    g = genus.lower().strip()
    return m == g or m.startswith(g + " ") or ("[" + g + "]") in m or (" " + g + " ") in m

matched_mm = []
matched_mg = []
linked_metabolites = set()
linked_genes_direct = set()
linked_genes_via_meta = set()

for g in sorted(cand_genera):
    for row in mm:
        if match_genus(row.get("gut_microbiota", ""), g):
            matched_mm.append((g, row))
            if row.get("metabolite") and row["metabolite"] != "nan":
                linked_metabolites.add(row["metabolite"])
    for row in mg:
        if match_genus(row.get("gut_microbiota", ""), g):
            matched_mg.append((g, row))
            if row.get("gene") and row["gene"] != "nan":
                linked_genes_direct.add(row["gene"])

matched_tg = []
for row in tg:
    if row.get("metabolite") in linked_metabolites:
        matched_tg.append(row)
        if row.get("gene") and row["gene"] != "nan":
            linked_genes_via_meta.add(row["gene"])

all_linked_genes = linked_genes_direct | linked_genes_via_meta
print(f"Matched GMrepo genera -> gutMGene microbe_metabolite edges: {len(matched_mm)} (unique metabolites={len(linked_metabolites)})")
print(f"Matched GMrepo genera -> gutMGene microbe_gene edges: {len(matched_mg)} (direct genes={len(linked_genes_direct)})")
print(f"Linked metabolites -> gutMGene metabolite_gene edges: {len(matched_tg)} (metabolite-linked genes={len(linked_genes_via_meta)})")
print(f"Total union target genes anchored to GMrepo BC differential genera: {len(all_linked_genes)}")
print("Matched genera with edges:", sorted({g for g, _ in matched_mm} | {g for g, _ in matched_mg}))
print("Direct genes:", sorted(linked_genes_direct))
print("Metabolites bridging to genes:", sorted({r['metabolite'] for r in matched_tg}))

# Save candidate target genes anchored to GMrepo
out_path = os.path.join(res_dir, "GMrepo_anchored_gutMGene_targets.csv")
with open(out_path, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f)
    w.writerow(["gene", "via_direct_microbe", "via_metabolite", "linked_genera", "linked_metabolites"])
    for gene in sorted(all_linked_genes):
        genera_d = sorted({g for g, r in matched_mg if r.get("gene") == gene})
        metas_i = sorted({r["metabolite"] for r in matched_tg if r.get("gene") == gene})
        genera_i = sorted({g for g, r in matched_mm if r.get("metabolite") in metas_i})
        w.writerow([
            gene,
            gene in linked_genes_direct,
            gene in linked_genes_via_meta,
            "; ".join(sorted(set(genera_d) | set(genera_i))),
            "; ".join(metas_i),
        ])
print("Saved:", out_path)
