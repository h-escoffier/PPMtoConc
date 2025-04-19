import pandas as pd


# TODO: Allow the possibility to give a df instead of a file path


def convert_data(file_path, save_path, total_prot_content=0.505717472): 
    """
    The resulting DataFrame contains the following columns:
        - ENSPID: Ensembl Protein ID
        - UniProtID: UniProt Accession ID
        - MolecularWeight: Molecular Weight in kDa
        - Abundance: Abundance in ppm
        - Mass_g_per_gDCW: Mass in grams per gram of dry cell weight (gDCW)
    
    Parameters:
        file_path (str): Path to the input file containing proteomics data.
        total_prot_content (float): Total protein content in gDCW. Default is 0.505717472 (Average across NCI60 cell lines). 
        save_path (str, optional): Path to save the converted data. If None, the data will not be saved.
    """
    df = pd.read_csv(file_path)
    df['g'] = df.apply(lambda row: conv_ppm_to_g(row['abundance'], row['MolecularWeight']),
                           axis=1
                           )
    df['g_gDCW'] = (df['g'] / df['g'].sum()) * total_prot_content
    df = df[['ENSPID', 'UniProtID', 'MolecularWeight', 'abundance', 'g_gDCW']]
    df = df.rename(columns={
        'ENSPID': 'ENSP',
        'UniProtID': 'UniProt',
        'MolecularWeight': 'MW_kDa',
        'abundance': 'Abundance',
        'g_gDCW': 'Mass_g_per_gDCW'
    })

    if save_path:
        df.to_csv(save_path, index=False)


def conv_ppm_to_g(ppm, molecular_weight):
    """
    
    """
    relative_abundance = ppm / 1e6
    relative_g = (relative_abundance / 6.022e23) * molecular_weight
    return relative_g

