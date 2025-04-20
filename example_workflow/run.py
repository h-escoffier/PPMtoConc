from ppm_to_conc.io import load_file, get_protein_molecular_weight, convert_ppm_dataframe


def run(input_file, output_file):
    df = load_file(input_file)
    df_with_mw = get_protein_molecular_weight(df)
    df_converted = convert_ppm_dataframe(df_with_mw)
    df_converted.to_csv(output_file, index=False)


if __name__ == "__main__":
    print('start')
    run('data/9606-iBAQ_MCF7_Geiger_2012_uniprot.txt', 'output/MCF7_test.txt')
    print('end')