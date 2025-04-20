from bioservices import UniProt
import requests
import time


def map_enspid_to_uniprot(list_enspid, batch_size=500, sleep_time=1):
    u = UniProt()
    mapping = {}
    for i in range(0, len(list_enspid), batch_size):
        batch = list_enspid[i:i + batch_size]
        res = u.mapping(fr="Ensembl_Protein", to="UniProtKB", query=batch)
        for prot in res['results']:
            mapping[prot['from']] = prot['to']['primaryAccession']
        time.sleep(sleep_time)
    return mapping


def fetch_uniprot_molecular_weight(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    response = requests.get(url)
    if response.ok:
        return response.json()['sequence']['molWeight']
    return None


def fetch_protein_sequence(ensp_id, species=9606):
    server = "https://rest.ensembl.org"
    ext = f"/sequence/id/{ensp_id}?type=protein"
    headers = {"Content-Type": "application/json"}
    response = requests.get(server + ext, headers=headers)
    data = response.json()
    return data.get('seq')
