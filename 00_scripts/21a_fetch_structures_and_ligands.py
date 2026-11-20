import csv
import json
import os
import requests

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

dock_dir = os.path.join(project_root, "01_data", "docking")
res_dir  = os.path.join(project_root, "02_results")
os.makedirs(os.path.join(dock_dir, "proteins"), exist_ok=True)
os.makedirs(os.path.join(dock_dir, "ligands"), exist_ok=True)

# 1. Fetch verified human protein structures (highest-res X-ray PDB w/ ligand or AlphaFold EBI v4)
target_proteins = ["IDO1", "PEMT", "ENO1", "FFAR2", "TGFBR1", "CASP1", "NLRP3", "DAPK3", "TCF7L2", "MAGI3", "MUC1"]

s = requests.Session()
s.trust_env = True

prot_manifest = []
for gene in target_proteins:
    u_url = f"https://rest.uniprot.org/uniprotkb/search?query=gene_exact:{gene}+AND+organism_id:9606+AND+reviewed:true&format=json&size=1"
    r = s.get(u_url, timeout=30)
    r.raise_for_status()
    res = r.json().get("results", [])
    if not res:
        print(f"WARN no UniProt entry for {gene}")
        continue
    entry = res[0]
    acc = entry["primaryAccession"]
    pdb_refs = [x for x in entry.get("uniProtKBCrossReferences", []) if x.get("database") == "PDB"]
    # Pick highest-resolution X-ray crystal structure if available
    best_pdb = None
    best_res = 99.0
    for pr in pdb_refs:
        props = {p["key"]: p["value"] for p in pr.get("properties", [])}
        if props.get("Method") == "X-ray":
            try:
                rv = float(props.get("Resolution", "99").replace(" A", "").strip())
            except Exception:
                rv = 99.0
            if rv < best_res:
                best_res = rv
                best_pdb = pr["id"]

    # Always download AlphaFold structure as fallback/clean full-length reference; download crystal PDB if <= 2.5A
    af_path = os.path.join(dock_dir, "proteins", f"{gene}_{acc}_AF.pdb")
    pdb_path = os.path.join(dock_dir, "proteins", f"{gene}_{best_pdb}.pdb") if best_pdb else None

    if best_pdb and best_res <= 2.5:
        if not os.path.exists(pdb_path) or os.path.getsize(pdb_path) < 1000:
            pr_dl = s.get(f"https://files.rcsb.org/download/{best_pdb}.pdb", timeout=45)
            if pr_dl.status_code == 200 and len(pr_dl.text) > 1000:
                with open(pdb_path, "w", encoding="utf-8") as f:
                    f.write(pr_dl.text)
        source_type = f"RCSB_PDB:{best_pdb} ({best_res:.2f}A)"
        chosen = pdb_path
    else:
        chosen = af_path
        source_type = f"AlphaFold:{acc}"

    # Ensure AlphaFold API query works for membrane/unresolved targets (PEMT, FFAR2, MAGI3, etc.)
    if chosen == af_path or not os.path.exists(chosen):
        af_api = s.get(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", timeout=30)
        if af_api.status_code == 200 and af_api.json():
            pdb_url = af_api.json()[0].get("pdbUrl")
            if pdb_url:
                p_dl = s.get(pdb_url, timeout=45)
                p_dl.raise_for_status()
                with open(af_path, "w", encoding="utf-8") as f:
                    f.write(p_dl.text)
                chosen = af_path
                source_type = f"AlphaFold:{acc}"

    print(f"  Protein {gene:8s} -> {acc:8s} | {source_type:25s} | saved={os.path.exists(chosen)}", flush=True)
    prot_manifest.append({"gene": gene, "uniprot": acc, "best_pdb": best_pdb or "", "resolution": best_res if best_pdb else "", "source": source_type, "pdb_path": os.path.relpath(chosen, project_root).replace(os.sep, "/")})

with open(os.path.join(res_dir, "Docking_Protein_Manifest.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["gene", "uniprot", "best_pdb", "resolution", "source", "pdb_path"])
    w.writeheader()
    w.writerows(prot_manifest)

# 2. Fetch canonical SMILES via PubChem REST for key gutMGene2 microbial metabolites
metabolites = {
    "Indole-3-carboxylic_acid": "69867",
    "Urolithin_A": "5488186",
    "Indole": "798",
    "3-Indolepropionic_acid": "3744",
    "Butyrate": "264",
    "Propionate": "1032",
    "Caffeic_acid": "689043",
    "3-(4-Hydroxyphenyl)propionic_acid": "10394",
    "Trimethylamine_N-oxide": "1145",
    "Lithocholic_acid": "9903",
    "Nicotinic_acid": "938"
}
# Verify & fetch IsomericSMILES directly from PubChem PUG REST
lig_rows = []
cids = ",".join(metabolites.values())
pr = s.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cids}/property/IsomericSMILES,MolecularWeight,MolecularFormula/JSON", timeout=30)
pr.raise_for_status()
props_map = {str(item["CID"]): item for item in pr.json()["PropertyTable"]["Properties"]}

for name, cid in metabolites.items():
    pinfo = props_map.get(cid, {})
    smi = pinfo.get("IsomericSMILES") or pinfo.get("SMILES") or pinfo.get("ConnectivitySMILES", "")
    mw  = pinfo.get("MolecularWeight", "")
    mf  = pinfo.get("MolecularFormula", "")
    lig_rows.append({"metabolite": name, "pubchem_cid": cid, "smiles": smi, "mw": mw, "formula": mf})
    print(f"  Ligand {name:34s} CID={cid:8s} MW={mw:7s} SMILES={smi}", flush=True)

with open(os.path.join(res_dir, "Docking_Ligand_Manifest.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["metabolite", "pubchem_cid", "smiles", "mw", "formula"])
    w.writeheader()
    w.writerows(lig_rows)
print("SAVED protein & ligand docking manifests!", flush=True)
