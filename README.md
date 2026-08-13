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

#### Requirements

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

#### Input

The script expects clonotyped TSV files whose names end with:

```text
YClon_clonotyped.tsv
```

Before running the script, edit the input directory:

```r
input_dir <- "path/to/clonotyped_files"
```

#### Run

From R or RStudio:

```r
source("clonotyped_gene_frequency.R")
```

#### Output

For each input sample, the script generates an Excel file containing:

* V gene absolute and relative frequencies;
* J gene absolute and relative frequencies;
* mean V gene annotation identity;
* mean J gene annotation identity.

---

### `equcab2e3.ndm.imgt`

Custom NDM file used for equine immunoglobulin sequence annotation.

This file contains species-specific information required for the analysis of equine immunoglobulin rearrangements. It's a mix between Eq2 and Eq3

---
