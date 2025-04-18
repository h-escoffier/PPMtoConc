import pandas as pd
from tqdm import tqdm
from api_mw import get_uniprot_from_enspid, get_uniprot_molecular_weight

def format_proteomics_data(file_path, save_path=None): 
    """
    Reads human proteomics data from PaxDB (https://pax-db.org/species/9606), cleans the string_external_id column, 
    and converts abundance to mmol/gDW.

    Parameters:
        file_path (str): Path to the input file containing proteomics data.
        save_path (str, optional): Path to save the cleaned data. If None, the data will not be saved.
    
    Returns:
        pd.DataFrame: A DataFrame containing the cleaned proteomics data with columns:
            - EnsemblID
            - EntrezID
            - ENSPID
            - MW
            - abundance
    """
    # Load data
    df = pd.read_csv(file_path, comment='#', sep='\t', names=["string_external_id", "abundance"])
    df["ENSPID"] = df["string_external_id"].str.replace("9606.", "", regex=False)  # TODO: Add species number in the parameter

    # Get unique ENSPIDs
    unique_enspid = df["ENSPID"].unique().tolist()
    # Map to UniProt
    print('Mapping Ensembl protein IDs to UniProt accession IDs...')
    enspid_to_uniprot = get_uniprot_from_enspid(unique_enspid)
    print('Mapping completed.')
    # Retrieve MWs
    data = []
    for enspid, uniprot_id in tqdm(iterable=enspid_to_uniprot.items(), desc="Retrieving molecular weights.."): 
        mw = get_uniprot_molecular_weight(uniprot_id)
        data.append({"ENSPID": enspid, "UniProtID": uniprot_id, "MolecularWeight": mw})
    print('Molecular weights retrieved.')
    df_mw = pd.DataFrame(data)
    df_final = df.merge(df_mw, on="ENSPID", how="left")

    if save_path:
        df_final.to_csv(save_path, index=False)

    return df_final[["ENSPID", "UniProtID", "MolecularWeight"]]


test_file_path = "input/9606-iBAQ_MCF7_Geiger_2012_uniprot.txt"
save_path = "output/cleaned_proteomics_data.csv"
formatted_df = format_proteomics_data(test_file_path, save_path)
print(formatted_df.head())