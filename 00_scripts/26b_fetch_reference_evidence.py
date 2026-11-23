"""Public literature and terms/policies retrieval; metadata only."""
import importlib.util
import json
from pathlib import Path

spec = importlib.util.spec_from_file_location('a', Path(__file__).with_name('26b_external_cohort_eligibility.py'))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

PMIDS = [
    '22522925',  # Curtis et al. Nature 2012 (METABRIC)
    '27161491',  # Pereira et al. Nat Commun 2016 (METABRIC)
    '30867590',  # Rueda et al. Nature 2019 (METABRIC relapse)
    '21558518',  # Hatzis et al. JAMA 2011 (GSE25066)
    '20490655',  # Sabatier et al. BCRT 2011 (GSE21653)
    '25887482',  # Jezequel et al. BCR 2015 (GSE58812)
    '31101122',  # Jezequel et al. BCR 2019 (GSE103091)
    '28650338',  # Ayers et al. JCI 2017 (18-gene T-cell-inflamed GEP)
    '29929551',  # Danaher et al. JITC 2018 (TIS / 18-gene GEP review/TCGA)
]

for pmid in PMIDS:
    p = m.fetch(
        f'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:{pmid}%20AND%20SRC:MED&format=json',
        f'literature/PMID{pmid}.json',
        2 * 1024 ** 2,
        terms='EuropePMC public REST metadata API'
    )
    j = json.loads(p.read_text(encoding='utf-8'))
    for x in j.get('resultList', {}).get('result', []):
        if x.get('pmcid') and x.get('isOpenAccess') == 'Y':
            pmcid = x['pmcid']
            try:
                m.fetch(
                    f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id={pmcid}&retmode=xml',
                    f'literature/{pmcid}.xml',
                    20 * 1024 ** 2,
                    terms='NCBI PMC Open Access Subset via E-utilities'
                )
            except Exception as e:
                print(pmcid, e, flush=True)

for url, name in [
    ('https://docs.cbioportal.org/user-guide/faq/', 'terms/cbioportal_faq.html'),
    ('https://ftp.ncbi.nlm.nih.gov/geo/README.txt', 'terms/geo_ftp_readme.txt'),
    ('https://www.ncbi.nlm.nih.gov/home/about/policies/', 'terms/ncbi_policies.html'),
]:
    try:
        m.fetch(url, name, 10 * 1024 ** 2, terms='Public repository terms/documentation page')
    except Exception as e:
        print(name, e, flush=True)
