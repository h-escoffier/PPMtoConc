import requests
import time
import pandas as pd 
from bioservices import UniProt


def get_uniprot_from_enspid(list_enspid, batch_size=500, sleep_time=1):
    """
    Maps a list of Ensembl protein IDs to UniProt IDs in batches to avoid slowdowns and server throttling.
    
    Parameters:
        list_enspid (list): A list of Ensembl protein IDs.
        batch_size (int): Number of IDs to send per request.
        sleep_time (int): Pause between requests to avoid hitting server limits.
    
    Returns:
        dict: Mapping {ENSPID: UniProtID}
    """
    u = UniProt()
    mapping = {}

    for i in range(0, len(list_enspid), batch_size):
        batch = list_enspid[i:i + batch_size]
        res = u.mapping(fr="Ensembl_Protein", to="UniProtKB", query=batch)
        for prot in res['results']:
            ensp = prot['from']
            uniprot_id = prot['to']['primaryAccession']
            mapping[ensp] = uniprot_id
        time.sleep(sleep_time)  

    return mapping


def get_uniprot_molecular_weight(uniprot_id):
    """
    Retrieves the molecular weight of a protein from the UniProt database.

    Parameters:
        uniprot_id (str): The UniProt accession ID of the protein.

    Returns:
        int: The molecular weight of the protein in Dalton.
        None: If the UniProt ID is not found or the request fails.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    response = requests.get(url)
    if response.ok:
        data = response.json()
        return data['sequence']['molWeight']
    else:
        return None


# def run_api(input_data): 
#     """
#     Retrieves molecular weights for a list of Ensembl protein IDs.
#     Parameters:
#         input_data (list): A list of Ensembl protein IDs.
#     Returns:
#         pd.DataFrame: A DataFrame containing Ensembl protein IDs, UniProt accession IDs, and their molecular weights (Da).
#     """
#     lst_mw = []
#     lst_uniprot = get_uniprot_from_enspid(input_data)
#     for prot in lst_uniprot:
#         mw = get_uniprot_molecular_weight(prot)
#         lst_mw.append(mw)
#     df = pd.DataFrame(lst_mw, columns=['MolecularWeight'])
#     df['UniProtID'] = lst_uniprot
#     df['EnsemblProteinID'] = input_data
#     df = df[['EnsemblProteinID', 'UniProtID', 'MolecularWeight']]
#     return df
        

# Test
# output_df = run_api(input_data = ["ENSP00000340466", "ENSP00000352438"])
# print(output_df)
