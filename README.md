# Equine IGL Repertoire Analysis
This repository contains scripts and configuration files used for the analysis of the equine immunoglobulin lambda (IGL) repertoire.

## Repository contents

### `clonotyped_gene_frequency.R`
R script used to calculate V and J gene usage from clonotyped sequence files.

### `equcab2e3.ndm.imgt`
Custom NDM file used for equine immunoglobulin sequence annotation.
This file contains species-specific information required for the analysis of equine immunoglobulin rearrangements. It's a mix between Eq2 and Eq3

### preprocess.sh
This file was based on pRESTO tool, according to described by Vander-Heiden et al., 2014

### id_adjusted.py
This py script was used to calculate the identity adjusted by coverage, based to described by Price and Arkin, 2019. 
Its a math formula which combines the identity (IgBLAST output) and the coverage (which is the proportion of the gremlin segment covered by the alignment)

---
