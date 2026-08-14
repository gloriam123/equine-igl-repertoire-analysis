# Equine IGL Repertoire Analysis

This repository contains scripts and configuration files used for the analysis of the equine immunoglobulin lambda (IGL) repertoire.

## Repository contents

### `clonotyped_gene_frequency.R`

R script used to calculate V and J gene usage from clonotyped sequence files.

The script:

* processes all files ending in `YClon_clonotyped.tsv` (I used YClon: it's very good and fast for processing light chain)s;
* calculates absolute and relative frequencies of V and J gene usage;
* calculates the mean annotation identity for each V and J gene (important to know how similar is your sequences compared to the gemrline);
* exports the results to separate worksheets in an Excel file for each sample (it's good to keep it organized).

### `equcab2e3.ndm.imgt`

Custom NDM file used for equine immunoglobulin sequence annotation.

This file contains species-specific information required for the analysis of equine immunoglobulin rearrangements. It's a mix between Eq2 and Eq3

### preprocess.sh

### id_adjusted
---
