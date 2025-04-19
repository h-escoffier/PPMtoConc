import pandas as pd
from tqdm import tqdm
from api_mw import get_uniprot_from_enspid, get_uniprot_molecular_weight, get_protein_sequence, calculate_molecular_weight


def format_proteomics_data(file_path, save_path=None): 
    """
    Reads human proteomics data from PaxDB (https://pax-db.org/species/9606), cleans the string_external_id column, 
    and converts abundance to mmol/gDW.

    Parameters:
        file_path (str): Path to the input file containing proteomics data.
        save_path (str, optional): Path to save the cleaned data. If None, the data will not be saved.
    
    Returns:
        pd.DataFrame: A DataFrame containing the cleaned proteomics data with columns:
            - ENSPID
            - UniProtID
            - MolecularWeight
            - Abundance
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
    for enspid in tqdm(unique_enspid, desc="Retrieving molecular weights.."): 
        uniprot_id = enspid_to_uniprot.get(enspid, None)
        mw = None

        if uniprot_id:  # If mapping exists, try UniProt method
            mw = get_uniprot_molecular_weight(uniprot_id)

        if mw is None or pd.isna(mw):  # Fallback: sequence from Ensembl
            seq = get_protein_sequence(enspid)
            if seq:
                mw = calculate_molecular_weight(seq)

        data.append({
            "ENSPID": enspid,
            "UniProtID": uniprot_id,
            "MolecularWeight": mw
        })

    # Create DataFrame
    df_mw = pd.DataFrame(data)
    df_final = df.merge(df_mw, on="ENSPID", how="left")
    df_final = df_final[["ENSPID", "UniProtID", "MolecularWeight"]]

    if save_path:
        df_final.to_csv(save_path, index=False)

    return df_final


test_file_path = "input/9606-iBAQ_MCF7_Geiger_2012_uniprot.txt"
save_path = "output/cleaned_proteomics_data_2.csv"
formatted_df = format_proteomics_data(test_file_path, save_path)
print(formatted_df.head())