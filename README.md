# SMLR-Genotyping
This repository contains all analysis, figure, and table generation scripts for the article: \
'Enhanced SNP Genotyping with Symmetric Multinomial Logistic Regression' \
(preprint submitted to bioRxiv, DOI: https://doi.org/10.1101/2024.11.28.625807)

(submitted for review at Forensic Science International: Genetics) \
[DOI: to be updated if accepted]

The repository is published on GitHub and preserved on Zenodo. \
[https://github.com/maltebn/SMLR-Genotyping.git] \
[DOI: 10.5281/zenodo.14176720]

## Overview of scripts
### 00_functions_and_global_definitions.R
Defines functions and global objects to be used throughout the scripts.

### 00-set-up-data.Rmd
Prepares the raw data files for further analysis.

### 01_bootstrap.R
Generates objects with bootstrap results (and writes them to the disk).
Best to be run on a server.

### 01_crossval_custom.R
Generates objects with cross-validation results (and writes them to the disk).
Best to be run on a server.

### 02_Bootstrap_FSI_final.Rmd
Analyse bootstrap results and generates Supplementary Figure S1.

### 02_cross_validations_FSI_final.Rmd
Analyse cross-validation results and generates Supplementary Figure S2.
As part of this, the additional data objects are written to the disk.

### 03_accuracy-and-call-rate_FSI_final.Rmd
Analyse cross-validation results and generates Figure 3 and Figure 4.

### 03_plot_HID_all_relevant_dils_FSI_final.Rmd
Generates Figure 1 and Figure 2.

### 03-summary-statistics-and-tables_FSI_final.Rmd
Generates Table 1 and Supplementary Figure S3.

### 04-convert-bmp-to-jpeg.sh
Converts figures from the bmp-format to jpeg-format.
While R is capable of directly saving plots as jpeg-files, this approach unfortunately made the subscripts on some axis titles look weird.
This was also the case for the possible png- and pdf-format.
The bmp-format was found to produce the desired output, but LaTeX does not support this file format.
Therefore, the figures generated in the above scripts are first save as bmp-files and then converted to jpeg-files.
The Bash script relies on ImageMagick having been installed.
