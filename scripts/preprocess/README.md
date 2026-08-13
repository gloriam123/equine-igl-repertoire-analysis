# Immunoglobulin Repertoire Pre-processing

This directory contains the pRESTO-based pipeline used to preprocess paired-end immunoglobulin repertoire sequencing reads before germline annotation.
It is important to mention that this tool was based on described by: 
Vander Heiden JA, Yaari G, Uduman M, Stern JN, O'Connor KC, Hafler DA, Vigneault F, Kleinstein SH. 
pRESTO: a toolkit for processing high-throughput sequencing raw reads of lymphocyte receptor repertoires. 
Bioinformatics. 2014 
Jul 1;30(13):1930-2. 
doi: 10.1093/bioinformatics/btu138. 
Epub 2014 Mar 10. 
PMID: 24618469; 
PMCID: PMC4071206.


## Script

`preprocess_repertoire_presto.sh`

## What the script does

The pipeline performs the following steps:

1. Assembles paired-end reads.
2. Filters sequences based on read quality (Q30).
3. Identifies forward and constant-region primers.
4. Organizes primer annotations in the sequence headers.
5. Filters sequences by length (minimum 100 nt).
6. Collapses duplicate sequences.
7. Separates sequences according to duplicate count.
8. Exports sequence information from the headers.
9. Converts FASTQ files to FASTA format for downstream germline analysis.

## Requirements

* pRESTO 0.6.2
* Bash
* Standard Unix command-line utilities (`cat`, `paste`, `sed`, `cut`, and `tr`)

pRESTO must be installed and its executables must be available in the system `PATH`.

## Input files

The script requires:

* paired-end FASTQ files (`R1` and `R2`);
* FASTA file containing the forward/V primers;
* FASTA file containing the constant-region primers.

Before running the script, define the input files and sample name:

```bash
SAMPLE="sample_name"

R1="sample_R1.fastq"
R2="sample_R2.fastq"

V_PRIMERS="V_primers.fasta"
C_PRIMERS="C_primers.fasta"
```

## Main parameters

The pipeline uses:

* quality threshold: Q30;
* minimum sequence length: 100 nt;
* maximum primer alignment error: 0.5;
* duplicate-count threshold: 2.

These parameters can be modified directly in the script if required.

## Run

```bash
bash preprocess_repertoire_presto.sh
```

## Output

The pipeline generates intermediate and final FASTQ files, processing logs, primer information, duplicate-count information, and FASTA files for downstream germline annotation.

Input sequencing files are not included in this repository.
