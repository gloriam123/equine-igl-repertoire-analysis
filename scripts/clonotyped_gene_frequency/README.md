# Clonotyped Gene Frequency Analysis

This directory contains the R script used to calculate V and J gene frequencies from clonotyped immunoglobulin sequence files.

## Script

`clonotyped_gene_frequency.R`

## What the script does

The script:

1. Searches the input directory for files ending in `YClon_clonotyped.tsv`.
2. Processes all matching files in a loop.
3. Calculates the absolute frequency of each V gene assignment.
4. Calculates the relative frequency of each V gene assignment.
5. Calculates the absolute frequency of each J gene assignment.
6. Calculates the relative frequency of each J gene assignment.
7. Calculates the mean V annotation identity for each V gene.
8. Calculates the mean J annotation identity for each J gene.
9. Exports the results to an Excel workbook for each sample.

Gene frequencies are calculated from the clonotype records contained in each input file.

## Requirements

* R
* `data.table`
* `tidyverse`
* `dplyr`
* `openxlsx`

Install the required packages in R if necessary:

```r
install.packages(c(
  "data.table",
  "tidyverse",
  "dplyr",
  "openxlsx"
))
```

## Input files

The script searches recursively for files whose names end with:

```text
YClon_clonotyped.tsv
```

The input files must contain at least the following columns:

```text
sequence_id
v_call
v_identity
j_call
j_identity
```

Before running the analysis, define the directory containing the clonotyped files:

```r
input_dir <- "path/to/clonotyped_files"
```

## Run

From R or RStudio:

```r
source("clonotyped_gene_frequency.R")
```

Alternatively, from the command line:

```bash
Rscript clonotyped_gene_frequency.R
```

## Output

For each input file, the script generates an Excel workbook containing four worksheets:

* `V_Gene_Frequency`
* `J_Gene_Frequency`
* `V_Gene_Identity`
* `J_Gene_Identity`

The frequency worksheets contain both absolute and relative gene frequencies.

The identity worksheets contain the mean annotation identity and the number of records associated with each V or J gene.

Relative frequencies are reported as proportions between 0 and 1.

