import argparse
from ppm_to_conc.io import load_paxdb, enrich_protein_data, convert_ppm_dataframe


def main():
    parser = argparse.ArgumentParser(description="Convert PPM proteomics data to concentration units.")
    parser.add_argument("input_file", type=str, help="Path to input file.")
    parser.add_argument("output_file", type=str, help="Path to output file.")
    parser.add_argument("--species", type=int, default=9606, help="NCBI species identifier.")
    parser.add_argument("--total_protein_content", type=float, default=0.505717472, help="Total protein content in g/gDCW.")

    args = parser.parse_args()

    df = load_paxdb(args.input_file, species=args.species)
    enriched = enrich_protein_data(df)
    converted = convert_ppm_dataframe(enriched, total_prot_content=args.total_protein_content)
    converted.to_csv(args.output_file, index=False)


if __name__ == "__main__":
    main()
