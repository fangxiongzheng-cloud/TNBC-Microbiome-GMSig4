import csv
import os
from collections import Counter

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

gut_dir = os.path.join(project_root, "01_data", "gutMGene")
out_dir = os.path.join(project_root, "02_results")
os.makedirs(out_dir, exist_ok=True)

mm = list(csv.DictReader(open(os.path.join(gut_dir, "gutMGene2_human_microbe_metabolite.csv"), encoding="utf-8-sig")))
mg = list(csv.DictReader(open(os.path.join(gut_dir, "gutMGene2_human_microbe_gene.csv"), encoding="utf-8-sig")))
tg = list(csv.DictReader(open(os.path.join(gut_dir, "gutMGene2_human_metabolite_gene.csv"), encoding="utf-8-sig")))

all_genes = set()
gene_sources = {}
for row in mg:
    g = row.get("gene", "").strip()
    if g and g != "nan":
        all_genes.add(g)
        gene_sources.setdefault(g, {"microbes": set(), "metabolites": set(), "modes": set(), "pmids": set()})
        gene_sources[g]["microbes"].add(row.get("gut_microbiota", ""))
        gene_sources[g]["modes"].add(row.get("associative_modes", ""))
        gene_sources[g]["pmids"].add(row.get("pmid", ""))

for row in tg:
    g = row.get("gene", "").strip()
    if g and g != "nan":
        all_genes.add(g)
        gene_sources.setdefault(g, {"microbes": set(), "metabolites": set(), "modes": set(), "pmids": set()})
        gene_sources[g]["metabolites"].add(row.get("metabolite", ""))
        gene_sources[g]["modes"].add(row.get("associative_modes", ""))
        gene_sources[g]["pmids"].add(row.get("pmid", ""))

print(f"Total unique human target genes in gutMGene v2.0: {len(all_genes)}")
both = [g for g, v in gene_sources.items() if v["microbes"] and v["metabolites"]]
print(f"Genes supported by BOTH direct microbe->gene and metabolite->gene: {len(both)}")
print("Top 20 metabolites in metabolite_gene:", Counter(r["metabolite"] for r in tg if r.get("metabolite")).most_common(20))
print("Top 20 microbes in microbe_gene:", Counter(r["gut_microbiota"] for r in mg if r.get("gut_microbiota")).most_common(20))
print("Top 20 microbes in microbe_metabolite:", Counter(r["gut_microbiota"] for r in mm if r.get("gut_microbiota")).most_common(20))

out_csv = os.path.join(out_dir, "gutMGene2_human_target_genes_master.csv")
with open(out_csv, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f)
    w.writerow(["gene", "n_microbes", "n_metabolites", "microbes", "metabolites", "associative_modes", "pmids"])
    for g in sorted(all_genes):
        v = gene_sources[g]
        w.writerow([
            g,
            len(v["microbes"]),
            len(v["metabolites"]),
            "; ".join(sorted(x for x in v["microbes"] if x and x != "nan")),
            "; ".join(sorted(x for x in v["metabolites"] if x and x != "nan")),
            "; ".join(sorted(x for x in v["modes"] if x and x != "nan")),
            "; ".join(sorted(x for x in v["pmids"] if x and x != "nan")),
        ])
print("Saved master target gene summary:", out_csv)
