
# Comparative evaluation of Olink Explore 3072 and mass spectrometry with peptide fractionation for plasma proteomics

## Overview

In this study, we performed a direct comparison of two plasma proteomics technologies, HiRIEF LC-MS/MS ("MS") and Olink Explore 3072 ("Olink"), across 88 plasma samples and 1,129 proteins analyzed with both methods. We evaluated and compared proteome coverage, technical precision, quantitative agreement at the protein and peptide level, and agreement and statistical power in detecting biological differences in protein abundance between groups, with sex as an example. In addition, we assessed the influence of various technical factors on quantitative agreement and compared the cross-platform correlations to those reported in previous studies. Finally, we developed an R Shiny app for exploring MS-Olink peptide-level correlations along the protein sequence and structure: https://peptaffinity.serve.scilifelab.se/app/peptaffinity. 

The results are published in Communications Chemistry:

Sissala, N., Babačić, H., Leo, I.R. et al. Comparative evaluation of Olink Explore 3072 and mass spectrometry with peptide fractionation for plasma proteomics. Commun Chem 8, 327 (2025). https://doi.org/10.1038/s42004-025-01753-2

## Content

The repository contains all R code necessary to reproduce the results presented in the paper. The scripts corresponding to each analysis can be found in the `/code` folder:

- `proteome_coverage.Rmd` – Comparison of the number of detected proteins, proportion of missing values, and proteome coverage.  
- `functional_annotation.R` – Comparison of Human Protein Atlas annotations, GO terms, and FDA-approved blood biomarkers among detected proteins.  
- `CV_analysis.Rmd` – Technical precision at the protein and peptide level.  
- `DAA.Rmd` – Comparison of differential abundance analysis results between platforms.  
- `correlation_analysis.Rmd` – MS–Olink correlations and association with technical factors and protein characteristics.  
- `correlation_comparisons.Rmd` – Comparison of correlations between MS, Olink, and SomaScan platforms across different studies.  
- `peptide_correlations.R` – Calculation of correlations between MS peptides and Olink assays.  
- `peptide_correlation_analysis.Rmd` – Visualization of peptide-level MS-Olink correlations mapped to protein sequences for select examples.  

These scripts produce results in the `/results` folder, with subfolders corresponding to each analysis:

```text
/results
├── proteome_coverage
├── functional_annotation
├── CV_analysis
├── DAA
├── correlation_analysis
├── correlation_comparisons
└── peptide_correlations
```


> **Note:**
> - The `DAA.Rmd` code cannot be run, as we cannot share participant data publicly. However, code and results are still provided.
> - Data and code for the PeptAffinity R Shiny application are deposited separately at https://github.com/isabelle-leo/PeptAffinity.  

Due to the large size of the data files, these are not stored in the github repository, but can be downloaded elsewhere. Below is a summary of files needed to re-run the analysis, and the structure of the `/data` folder. Details are provided under section "Running the code".

```text
/data
├── raw_data          # Raw MS and Olink data
├── processed_data    # Cleaned data files produced by `code/clean_data.R`
├── external_data     # External datasets used for cross-study comparisons (see Supplementary Data 10)
└── metadata          # Sample and protein metadata (participant metadata not included)
```

## Running the code

1. Clone the repository

	```bash
	git clone https://github.com/noorasissala/MS-Olink-comparison
	```

2. Download the raw data files from the PRIDE repository (https://www.ebi.ac.uk/pride/) and place them in `data/raw_data`.

   	- **HiRIEF LC-MS/MS**: PRIDE identifier `PXD061144`. Required files: `proteins_table.txt`, `peptides_table.txt`, and `target_psmtable.txt`.
   	- **Olink Explore 3072**: PRIDE identifier `PAD000006`. Required file: `VB-3207.npx.csv`

3. Download `sample_metadata.txt` (from `PXD061144`) and place in `data/metadata`

4. Run `clean_data.R`, which produces processed data files in `data/processed_data`

5. Run the analysis scripts.

	Suggested order to run analyses:

	- proteome_coverage.Rmd
	- functional_annotation.R
	- CV_analysis.Rmd
	- correlation_analysis.Rmd
	- correlation_comparisons.Rmd
	- peptide_correlations.R
	- peptide_correlation_analysis.Rmd (depends on output from peptide_correlations.R)

> **Note:**
> - Before running each script, ensure you have installed the required R packages (see `session_info.txt` in the corresponding `/results` subfolder).
> - `functional_annotation.R` optionally uses Supplementary Table 1 from "The Clinical Plasma Proteome: A Survey of Clinical Assays for Proteins in Plasma and Serum" by Anderson (Clinical Chemistry, 2010) to compare the detection of FDA approved blood biomarkers between platforms. The script can be run without it, but download it manually and place in `data/external_data` to fully reproduce the analysis.
> - To run `correlation_comparisons.Rmd`, you need to download external data from previous studies. A list of required files and sources is provided in **Supplementary Data 10** of the present publication. Place the files in `data/external_data`.

## Citation

If you use this code in your research, please cite the original publication and the code:

Article: Sissala, N., Babačić, H., Leo, I.R. et al. Comparative evaluation of Olink Explore 3072 and mass spectrometry with peptide fractionation for plasma proteomics. Commun Chem 8, 327 (2025). https://doi.org/10.1038/s42004-025-01753-2

Code: Sissala, N. MS-Olink-comparison. Zenodo https://doi.org/10.5281/zenodo.17246257 (2025).
