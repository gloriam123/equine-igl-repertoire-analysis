# Coverage-adjusted Identity Filtering

This directory contains the Python script used to calculate coverage-adjusted identity for V and J gene assignments and filter immunoglobulin sequences according to this metric.
This metric was modified based on:
Price MN, Arkin AP. Curated BLAST for Genomes. 
mSystems. 2019 Mar 26;4(2):e00072-19. 
doi: 10.1128/mSystems.00072-19. 
PMID: 30944879; 
PMCID: PMC6435814.

## Script

`id_adjusted.py`

## What the script does

The script:

1. Reads an AIRR-formatted TSV file.
2. Retains productive sequences.
3. Requires complete amino-acid annotation of FWR1, CDR1, FWR2, CDR2, FWR3, CDR3, and FWR4.
4. Uses the first V or J gene assignment when multiple gene calls are reported.
5. Retrieves the lengths of the assigned V and J germline genes from FASTA files.
6. Calculates the number of evaluated bases from the `M` operations in the V and J CIGAR strings.
7. Calculates V and J gene coverage.
8. Calculates coverage-adjusted identity for V and J.
9. Retains sequences with coverage-adjusted identity ≥ 80% for both V and J.
10. Exports the filtered sequences as a TSV file.

## Coverage-adjusted identity

Coverage-adjusted identity is calculated as:

```text
coverage-adjusted identity =
alignment identity × evaluated bases / germline gene length
```

Gene coverage is calculated as:

```text
gene coverage =
evaluated bases / germline gene length × 100
```

## Requirements

* Python 3
* pandas
* Biopython

Install the required Python packages with:

```bash
pip install pandas biopython
```

## Input files

The script requires:

* an AIRR-formatted TSV file containing sequence annotations;
* a V germline FASTA file;
* a J germline FASTA file.

Before running the script, define the input and output files:

```python
INPUT_FILE = "path/to/sample.airr.tsv"

V_FASTA = "path/to/V_germline.fasta"
J_FASTA = "path/to/J_germline.fasta"

OUTPUT_FILE = "sample_filtered_identity80.tsv"
```

The filtering threshold can also be modified if necessary:

```python
COVERAGE_ADJUSTED_IDENTITY_CUTOFF = 80.0
```

## Run

```bash
python id_adjusted.py
```

## Output

The output is a tab-separated file containing sequences that:

* are productive;
* have complete FR/CDR amino-acid annotation;
* have V coverage-adjusted identity ≥ 80%;
* have J coverage-adjusted identity ≥ 80%.

The output table also includes the calculated V and J gene coverage and coverage-adjusted identity values.

