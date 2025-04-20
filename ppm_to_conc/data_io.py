import pandas as pd
from tqdm import tqdm

from .api import map_enspid_to_uniprot, fetch_uniprot_molecular_weight, fetch_protein_sequence
from .conversion import calculate_sequence_mw, ppm_to_grams


def load_paxdb(file_path, species=9606):
    df = pd.read_csv(file_path, comment='#', sep='\t', names=["string_external_id", "abundance"])
    df["ENSPID"] = df["string_external_id"].str.replace(f"{species}.", "", regex=False)
    return df


def load_tsv(file_path): 
    df = pd.read_csv(file_path, sep="\t", header=0)
    return df


def check_file_format(file_path):
    with open(file_path, 'r') as file:
        first_line = file.readline().strip()
        if first_line.startswith("#"):
            return "paxdb"
        else:
            return "tsv"


def load_file(file_path, species=9606):
    file_format = check_file_format(file_path)
    if file_format == "paxdb":
        return load_paxdb(file_path, species)
    elif file_format == "tsv":
        return load_tsv(file_path)
    else:
        raise ValueError("Unsupported file format. Please provide a valid PAXdb or TSV file.")


def get_protein_molecular_weight(df, fill_missing=None):
    unique_enspids = df["ENSPID"].unique().tolist()
    enspid_to_uniprot = map_enspid_to_uniprot(unique_enspids)

    mw_records = []
    for enspid in tqdm(unique_enspids, desc="Retrieving molecular weights.."): 
        uniprot_id = enspid_to_uniprot.get(enspid)
        mw = fetch_uniprot_molecular_weight(uniprot_id) if uniprot_id else None
        if mw is None:
            seq = fetch_protein_sequence(enspid)
            if seq:
                mw = calculate_sequence_mw(seq)
        mw_records.append({
            "ENSPID": enspid,
            "UniProtID": uniprot_id,
            "MolecularWeight": mw
        })

    mw_df = pd.DataFrame(mw_records)

    if fill_missing in ["mean", "median"]:
        fill_value = mw_df["MolecularWeight"].mean() if fill_missing == "mean" else mw_df["MolecularWeight"].median()
        mw_df["MolecularWeight"] = mw_df["MolecularWeight"].fillna(fill_value)


    return df.merge(mw_df, on="ENSPID", how="left")


def convert_ppm_dataframe(df, total_prot_content=0.505717472):
    df['g'] = df.apply(lambda row: ppm_to_grams(row['abundance'], row['MolecularWeight']), axis=1)
    df['Mass_g_per_gDCW'] = (df['g'] / df['g'].sum()) * total_prot_content
    return df[['ENSPID', 'UniProtID', 'MolecularWeight', 'abundance', 'Mass_g_per_gDCW']]
