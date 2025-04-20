import argparse
from ppm_to_conc.io import load_file, get_protein_molecular_weight, convert_ppm_dataframe


def main():
    parser = argparse.ArgumentParser(description="Convert PPM proteomics data to concentration units.")
    parser.add_argument("input_file", type=str, help="Path to input file.")
    parser.add_argument("output_file", type=str, help="Path to output file.")
    parser.add_argument("--species", type=int, default=9606, help="NCBI species identifier.")
    parser.add_argument("--fill_missing", type=str, default=None, help="Fill missing molecular weights with 'mean' or 'median'.")
    parser.add_argument("--total_protein_content", type=float, default=0.505717472, help="Total protein content in g/gDCW.")

    args = parser.parse_args()

    df = load_file(args.input_file, species=args.species)
    df_with_mw = get_protein_molecular_weight(df, fill_missing=args.fill_missing)
    df_converted = convert_ppm_dataframe(df_with_mw, total_prot_content=args.total_protein_content)
    df_converted.to_csv(args.output_file, index=False)


if __name__ == "__main__":
    main()