"""Retrieve public processed matrices, clinical metadata, and bounded cBioPortal API mRNA data."""
import importlib.util
import json
from pathlib import Path

spec = importlib.util.spec_from_file_location('acquire', Path(__file__).with_name('26b_external_cohort_eligibility.py'))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

# 1. METABRIC clinical, license, case lists (GitHub / Git LFS)
for fn in ['data_clinical_patient.txt', 'data_clinical_sample.txt']:
    m.fetch(
        'https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/brca_metabric/' + fn,
        'METABRIC/resolved/' + fn,
        terms='ODC Open Database License (ODbL) 1.0 per brca_metabric/LICENSE'
    )

for fn in ['LICENSE', 'Readme.txt', 'meta_clinical_patient.txt', 'meta_clinical_sample.txt',
           'case_lists/cases_all.txt', 'case_lists/cases_mRNA.txt',
           'case_lists/cases_nature_2012.txt', 'case_lists/cases_nat_comm_2016.txt']:
    m.fetch(
        'https://raw.githubusercontent.com/cBioPortal/datahub/master/public/brca_metabric/' + fn,
        'METABRIC/' + fn,
        10 * 1024 ** 2,
        terms='ODC Open Database License (ODbL) 1.0 per brca_metabric/LICENSE'
    )

# 2. METABRIC bounded mRNA expression via cBioPortal REST API (full data_mrna_illumina_microarray.txt is ~657 MiB > 500 MiB cap)
m.fetch_api_json(
    'https://www.cbioportal.org/api/studies/brca_metabric',
    'METABRIC/api/brca_metabric_study.json',
    terms='ODC Open Database License (ODbL) 1.0; cBioPortal REST API'
)
m.fetch_api_json(
    'https://www.cbioportal.org/api/molecular-profiles/brca_metabric_mrna',
    'METABRIC/api/brca_metabric_mrna_profile.json',
    terms='ODC Open Database License (ODbL) 1.0; cBioPortal REST API'
)
sids_path = m.fetch_api_json(
    'https://www.cbioportal.org/api/sample-lists/brca_metabric_mrna/sample-ids',
    'METABRIC/api/brca_metabric_mrna_sample_ids.json',
    terms='ODC Open Database License (ODbL) 1.0; cBioPortal REST API'
)
sample_ids = json.loads(sids_path.read_text(encoding='utf-8'))

gene_records = {}
for g in m.ALL_PANEL_GENES:
    gp = m.fetch_api_json(
        f'https://www.cbioportal.org/api/genes/{g}',
        f'METABRIC/api/genes/{g}.json',
        terms='cBioPortal REST API gene lookup'
    )
    gene_records[g] = json.loads(gp.read_text(encoding='utf-8'))

entrez_ids = [gene_records[g]['entrezGeneId'] for g in m.ALL_PANEL_GENES]
m.fetch_api_json(
    'https://www.cbioportal.org/api/molecular-profiles/brca_metabric_mrna/molecular-data/fetch',
    'METABRIC/api/brca_metabric_target_gep_mrna.json',
    payload={'entrezGeneIds': entrez_ids, 'sampleIds': sample_ids},
    limit=100 * 1024 ** 2,
    terms='ODC Open Database License (ODbL) 1.0; cBioPortal REST API bounded target+GEP genes fetch'
)

# 3. GEO series matrices & subseries / overlap metadata
for gse, prefix in [('GSE25066', 'GSE25nnn'), ('GSE21653', 'GSE21nnn')]:
    m.fetch(
        f'https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/matrix/{gse}_series_matrix.txt.gz',
        f'{gse}/{gse}_series_matrix.txt.gz',
        150 * 1024 ** 2,
        terms='NCBI GEO public FTP repository; unrestricted academic research use with original citation'
    )

for gse, folder in [
    ('GSE25055', 'GSE25066'),
    ('GSE25065', 'GSE25066'),
    ('GSE83937', 'GSE103091'),
    ('GSE194040', 'reference_geo_metadata'),
    ('GSE173839', 'reference_geo_metadata'),
    ('GSE96058', 'reference_geo_metadata'),
]:
    m.fetch(
        f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&targ=self&form=text&view=full',
        f'{folder}/{gse}_metadata.txt',
        10 * 1024 ** 2,
        terms='NCBI GEO metadata SOFT query'
    )

# 4. GEO Platform GPL annotations (compressed official .annot.gz from NCBI FTP)
for gpl in ['GPL96', 'GPL570']:
    m.fetch(
        f'https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/{gpl}/annot/{gpl}.annot.gz',
        f'platforms/{gpl}.annot.gz',
        50 * 1024 ** 2,
        terms='NCBI GEO official platform annotation archive'
    )
