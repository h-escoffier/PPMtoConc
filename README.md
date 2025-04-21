# PPMtoConc

PPMtoConc is a tool for converting **proteomics data reported in parts per million (ppm)** into **cellular protein concentrations**, expressed as grams per gram of dry cell weight (g/gDCW). This conversion is useful for enzyme constrained metabolic modeling.

-----------------------

## Definition of PPM 

> The protein abundance is present in ppm which is short for parts per million. Abundance in 'ppm' is essentially describing each protein with reference to the entire expressed proteome.This means each protein entity is enumerated relative to all other protein molecules in the sample.
>
> From [PAXdb — Protein Abundance Database](https://pax-db.org/help).

-----------------------

## Conversion Formula

To integrate proteomics data into metabolic models, protein abundances measured in parts per million (PPM) are first converted into absolute mass fractions relative to the total cellular protein content.

The conversion process is as follows:

**1. PPM to Grams Conversion**

Each protein’s abundance in PPM is first converted into grams using its molecular weight:

```
mass_g = (abundance_ppm / 1e6) / Na * mw
```

Where:
* `abundance_ppm` is the protein's relative abundance (parts per million),
* `Na` is the Avogadro's number
* `mw` is the molecular weight of the protein given in Daltons.

This calculation gives the estimated mass contribution of each protein based on its abundance and molecular weight.

**2. Normalization to Total Protein Content**

The computed masses are normalized so that their sum matches the experimentally determined total protein fraction in the cell dry weight (g/gDCW). The default value corresponds to the average total protein content in the NCI60 cell lines.

```
Mass_g_per_gDCW = (mass_g / total_mass_g) * total_protein_content
```

Where:
* `total_mass_g` is the sum of all protein masses in the sample,
* `total_protein_content` is the expected total protein mass fraction in the biomass (default: 0.5057 g/gDCW).

-----------------------

## Usage 

### Installation 

```
git clone https://github.com/h-escoffier/PPMtoConc
cd PPMtoConc
pip install -r requirements.txt
```

### Testing

This repository includes unit tests using pytest.

```
pytest
```

This will automatically and run all tests under the tests/ directory.

### Command-Line Interface (CLI)

The tool can be used via the command line interface: 

```
python scripts/run_conversion.py <input_file> <output_file> [options]
``` 

#### CLI Arguments

* `input_file`: Path to the input file containing proteomics data (ppm values).
* `output_file`: Path to the output file to save the converted data.
* `--species`: Optional; NCBI species identifier (default: 9606).
* `--fill_missing`: Optional; Specify how to handle missing molecular weights. Use 'mean' or 'median' to fill missing values (default: None).
* `--total_protein_content`: Optional; Specify the total protein content in g/gDCW (default: 0.505717472).

#### Example Usage

```
python scripts/run_conversion.py data/9606-iBAQ_MCF7_Geiger_2012_uniprot.txt output_conc_data.tsv
```

-----------------------

### Function-Based Usage

The conversion can also be used in Python scripts by calling the `run()` function.

```
from ppm_to_conc.data_io import load_file, get_protein_molecular_weight, convert_ppm_dataframe


def run(input_file, output_file, species=9606, fill_missing=None, total_protein_content=0.505717472):
    df = load_file(input_file, species=species)
    df_with_mw = get_protein_molecular_weight(df, fill_missing=fill_missing)
    df_converted = convert_ppm_dataframe(df_with_mw, total_prot_content=total_protein_content)
    df_converted.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    print('start')
    run('data/9606-iBAQ_MCF7_Geiger_2012_uniprot.txt', 'output_conc_data.tsv')
    print('end')
```

-----------------------

## Input Format 

The input file should be either a `.txt` file download from the [PAXdb database](https://pax-db.org/download) or a `.tsv` file as follow : 

| ENSPID | Abundance_in_ppm |
| :--------------- | :--------------- |
| ENSP00000340466 | 459 |
| ENSP00000352438 | 454 |

-----------------------

## Output Format 

The output file will be a tab-delimited file with the same ENSPID and abundance values, with additional columns containing UniprotID, molecular weight and concentratioin in g/gDCW.

| ENSPID | UniProtID | MolecularWeight | abundance | Mass_g_per_gDCW |
| :--------------- | :--------------- | :--------------- | :--------------- | :--------------- | 
| ENSP00000340466 | Q14697 | 106874.0 | 459.0 | 0.00038304162818003304 |
| ENSP00000352438 | Q15366 | 38580.0 | 454.0 | 0.00013676636403379605 |

-----------------------

## Feedback & Improvements

Contributions, suggestions, and feedback are very welcome!
If you encounter any [issues](https://github.com/h-escoffier/PPMtoConc/issues), have ideas for new features, or notice room for improvement, feel free to open an issue or submit a pull request.

