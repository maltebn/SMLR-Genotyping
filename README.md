# SMLR-Genotyping
This repository contains all analysis, figure, and table generation scripts for the article: \
'Enhanced SNP genotyping with symmetric multinomial logistic regression' \
(preprint submitted to bioRxiv, DOI: https://doi.org/10.1101/2024.11.28.625807)

(peer-reviewed and accepted for publication at Forensic Science International: Genetics) \
[DOI: insert when received from publisher]

The repository is published on GitHub and preserved on Zenodo. \
[https://github.com/maltebn/SMLR-Genotyping.git] \
[DOI: 10.5281/zenodo.15341884]

## Overview of scripts
### 00-functions-and-global-definitions.R
Defines functions and global objects to be used throughout the scripts.

### 00-set-up-data.Rmd
Prepares the raw data files for further analysis.

### 01-MTK-models.Rmd
Exploration/selection of parameter settings for the observational model by Mostad, Tillmar, and Kling (MTK).
Using these to predict MTK genotypes within each dilution level.

### 02-bootstrap.R
Generates objects with bootstrap results.
Long run time: use server.

### 02-crossval-custom.R
Generates objects with cross-validation results.
Long run time: use server.

### 03-bootstrap-analysis.Rmd
Analyse bootstrap results and generates Supplementary Figure S1.

### 03-cross-validation-analysis.Rmd
Analyse cross-validation results and generates Supplementary Figure S2.
As part of this, the additional data objects are written to the disk.

### 04-accuracy-and-call-rate.Rmd
Analyse cross-validation results and generates Figure 3 and Figure 4.

### 04-plot-HSG-all-relevant-dils.Rmd
Generates Figure 1 and Figure 2.

### 04-summary-statistics-and-tables.Rmd
Generates Table 1, Supplementary Figure S3, and Supplementary Figure S4.

### 05-convert-bmp-to-jpeg.sh
Converts figures from the bmp-format to jpeg-format.
While R is capable of directly saving plots as jpeg-files, this approach unfortunately made the subscripts on some axis titles look weird.
This was also the case for the possible png- and pdf-format.
The bmp-format was found to produce the desired output, but LaTeX does not support this file format.
Therefore, the figures generated in the above scripts are first save as bmp-files and then converted to jpeg-files.
The Bash script relies on ImageMagick having been installed.
