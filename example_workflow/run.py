from ppm_to_conc.data_io import load_file, get_protein_molecular_weight, convert_ppm_dataframe


def run(input_file, output_file, species=9606, fill_missing=None, total_protein_content=0.505717472):
    df = load_file(input_file, species=species)
    df_with_mw = get_protein_molecular_weight(df, fill_missing=fill_missing)
    df_converted = convert_ppm_dataframe(df_with_mw, total_prot_content=total_protein_content)
    df_converted.to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":
    print('start')
    run('data/9606-iBAQ_MCF7_Geiger_2012_uniprot.txt', 'output/MCF7_test.txt')
    print('end')