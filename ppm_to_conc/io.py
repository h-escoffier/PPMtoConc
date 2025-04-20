import pandas as pd
from tqdm import tqdm
from .api import map_enspid_to_uniprot, fetch_uniprot_molecular_weight, fetch_protein_sequence
from .conversion import calculate_sequence_mw, ppm_to_grams


def load_paxdb(file_path, species=9606):
    df = pd.read_csv(file_path, comment='#', sep='\t', names=["string_external_id", "abundance"])
    df["ENSPID"] = df["string_external_id"].str.replace(f"{species}.", "", regex=False)
    return df


def enrich_protein_data(df):
    unique_enspid = df["ENSPID"].unique().tolist()
    mapping = map_enspid_to_uniprot(unique_enspid)

    enriched = []
    for enspid in tqdm(unique_enspid, desc="Retrieving molecular weights.."): 
        uniprot_id = mapping.get(enspid)
        mw = fetch_uniprot_molecular_weight(uniprot_id) if uniprot_id else None
        if mw is None:
            seq = fetch_protein_sequence(enspid)
            if seq:
                mw = calculate_sequence_mw(seq)
        enriched.append({
            "ENSPID": enspid,
            "UniProtID": uniprot_id,
            "MolecularWeight": mw
        })

    df_mw = pd.DataFrame(enriched)
    return df.merge(df_mw, on="ENSPID", how="left")


def convert_ppm_dataframe(df, total_prot_content=0.505717472):
    df['g'] = df.apply(lambda row: ppm_to_grams(row['abundance'], row['MolecularWeight']), axis=1)
    df['Mass_g_per_gDCW'] = (df['g'] / df['g'].sum()) * total_prot_content
    return df[['ENSPID', 'UniProtID', 'MolecularWeight', 'abundance', 'Mass_g_per_gDCW']]
