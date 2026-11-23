"""Download official structural QC references only; original docking files untouched."""
import pathlib, json, urllib.request, datetime, hashlib, csv
ROOT=pathlib.Path(r'D:\Bioinformatics\TNBC_Microbiome')
DATA=ROOT/'01_data/JTM_upgrade/docking'; OUT=ROOT/'02_results/JTM_upgrade/docking'
DATA.mkdir(parents=True,exist_ok=True); OUT.mkdir(parents=True,exist_ok=True)
jobs={
 '2D0T.pdb':'https://files.rcsb.org/download/2D0T.pdb',
 '4PK5.pdb':'https://files.rcsb.org/download/4PK5.pdb',
 '9FOH.cif':'https://files.rcsb.org/download/9FOH.cif',
 'A1IEE.cif':'https://files.rcsb.org/ligands/download/A1IEE.cif',
 '2D0T_entry.json':'https://data.rcsb.org/rest/v1/core/entry/2D0T',
 '4PK5_entry.json':'https://data.rcsb.org/rest/v1/core/entry/4PK5',
 '9FOH_entry.json':'https://data.rcsb.org/rest/v1/core/entry/9FOH',
 'PEMT_uniprot.json':'https://rest.uniprot.org/uniprotkb/Q9UBM1.json',
 'IDO1_uniprot.json':'https://rest.uniprot.org/uniprotkb/P14902.json',
}
rows=[]
for name,url in jobs.items():
 p=DATA/name; row={'filename':name,'url':url,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:
  with urllib.request.urlopen(url,timeout=45) as response: b=response.read()
  p.write_bytes(b); row.update(status='downloaded',bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),error='')
 except Exception as e: row.update(status='blocked_or_failed',bytes=0,sha256='',error=str(e))
 rows.append(row); print(row,flush=True)
with (OUT/'download_provenance.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
