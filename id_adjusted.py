#!/usr/bin/env python3

"""
Filter productive immunoglobulin sequences based on
coverage-adjusted V and J gene identity.

The script:
1. Reads an AIRR-formatted TSV file;
2. Retains productive sequences with complete FR/CDR annotation;
3. Retrieves V and J germline gene lengths from FASTA files;
4. Calculates the number of evaluated bases from the CIGAR strings;
5. Calculates V and J gene coverage;
6. Calculates coverage-adjusted identity;
7. Retains sequences with coverage-adjusted identity >= 80%
   for both V and J genes;
8. Exports the filtered sequences as a TSV file.

Requirements:
    pandas
    biopython
"""

import re
import pandas as pd
from Bio import SeqIO


# ============================================================
# CONFIGURATION
# ============================================================

INPUT_FILE = "path/to/sample.airr.tsv"

V_FASTA = "path/to/V_germline.fasta"
J_FASTA = "path/to/J_germline.fasta"

OUTPUT_FILE = "sample_filtered_identity80.tsv"

COVERAGE_ADJUSTED_IDENTITY_CUTOFF = 80.0


# AIRR amino-acid region columns required for complete annotation
REQUIRED_REGION_COLUMNS = [
    "fwr1_aa",
    "cdr1_aa",
    "fwr2_aa",
    "cdr2_aa",
    "fwr3_aa",
    "cdr3_aa",
    "fwr4_aa",
]


# ============================================================
# AUXILIARY FUNCTIONS
# ============================================================

def load_gene_lengths(
    fasta_path: str,
    gene_column: str,
    length_column: str
) -> pd.DataFrame:
    """
    Read a germline FASTA file and return the gene names
    and sequence lengths as a DataFrame.
    """

    records = []

    with open(fasta_path) as fasta_handle:

        for gene in SeqIO.parse(fasta_handle, "fasta"):

            records.append({
                gene_column: gene.id,
                length_column: len(gene.seq)
            })

    return pd.DataFrame(records)


def get_first_gene_call(call_value):
    """
    Return the first gene when multiple gene assignments
    are present.

    Example:
        IGLV8S1,IGLV8S2 -> IGLV8S1
    """

    if pd.isna(call_value):
        return pd.NA

    return str(call_value).split(",")[0].strip()


def sum_m_from_cigar(cigar):
    """
    Sum all values associated with the M operation
    in a CIGAR string.

    Examples:
        2S292M95S7N -> 292
        10M2I20M    -> 30

    This value is used here as the number of bases
    evaluated in the alignment.
    """

    if pd.isna(cigar):
        return pd.NA

    cigar = str(cigar)

    matches = re.findall(r"(\d+)(?=M)", cigar)

    if not matches:
        return 0

    return sum(int(value) for value in matches)


# ============================================================
# READ AIRR FILE
# ============================================================

print(f"Reading input file: {INPUT_FILE}")

sequences = pd.read_csv(
    INPUT_FILE,
    sep="\t",
    low_memory=False
)

print(f"Total input sequences: {len(sequences)}")


# ============================================================
# FILTER 1
# PRODUCTIVE SEQUENCES WITH COMPLETE FR/CDR ANNOTATION
# ============================================================

productive_mask = sequences["productive"] == "T"

complete_regions_mask = (
    sequences[REQUIRED_REGION_COLUMNS]
    .notna()
    .all(axis=1)
)

productive_complete = sequences[
    productive_mask & complete_regions_mask
].copy()

print(
    "After productive + complete-region filtering: "
    f"{len(productive_complete)}"
)


# ============================================================
# DEFINE UNIQUE V AND J GENE CALLS
# ============================================================

productive_complete["v_call_unique"] = (
    productive_complete["v_call"]
    .apply(get_first_gene_call)
)

productive_complete["j_call_unique"] = (
    productive_complete["j_call"]
    .apply(get_first_gene_call)
)


# ============================================================
# LOAD GERMLINE V AND J GENE LENGTHS
# ============================================================

print("Reading V germline gene lengths...")

v_gene_lengths = load_gene_lengths(
    fasta_path=V_FASTA,
    gene_column="v_gene",
    length_column="v_gene_length"
)


print("Reading J germline gene lengths...")

j_gene_lengths = load_gene_lengths(
    fasta_path=J_FASTA,
    gene_column="j_gene",
    length_column="j_gene_length"
)


# ============================================================
# ADD GERMLINE GENE LENGTHS TO THE AIRR TABLE
# ============================================================

sequences_with_lengths = pd.merge(
    productive_complete,
    v_gene_lengths,
    how="left",
    left_on="v_call_unique",
    right_on="v_gene"
)

sequences_with_lengths = pd.merge(
    sequences_with_lengths,
    j_gene_lengths,
    how="left",
    left_on="j_call_unique",
    right_on="j_gene"
)

print(
    "After merging germline gene lengths: "
    f"{len(sequences_with_lengths)}"
)


# ============================================================
# CALCULATE NUMBER OF EVALUATED BASES FROM CIGAR STRINGS
# ============================================================

sequences_with_lengths["v_evaluated_bases"] = (
    sequences_with_lengths["v_cigar"]
    .apply(sum_m_from_cigar)
)

sequences_with_lengths["j_evaluated_bases"] = (
    sequences_with_lengths["j_cigar"]
    .apply(sum_m_from_cigar)
)


# ============================================================
# ENSURE NUMERIC COLUMNS
# ============================================================

numeric_columns = [
    "v_identity",
    "j_identity",
    "v_evaluated_bases",
    "j_evaluated_bases",
    "v_gene_length",
    "j_gene_length",
]

for column in numeric_columns:

    sequences_with_lengths[column] = pd.to_numeric(
        sequences_with_lengths[column],
        errors="coerce"
    )


# ============================================================
# CALCULATE GERMLINE GENE COVERAGE
#
# coverage = evaluated bases / germline gene length * 100
# ============================================================

sequences_with_lengths["v_gene_coverage"] = (
    sequences_with_lengths["v_evaluated_bases"]
    / sequences_with_lengths["v_gene_length"]
    * 100
)

sequences_with_lengths["j_gene_coverage"] = (
    sequences_with_lengths["j_evaluated_bases"]
    / sequences_with_lengths["j_gene_length"]
    * 100
)


# ============================================================
# CALCULATE COVERAGE-ADJUSTED IDENTITY
#
# coverage-adjusted identity =
# alignment identity * evaluated bases / germline gene length
# ============================================================

sequences_with_lengths["v_coverage_adjusted_identity"] = (
    sequences_with_lengths["v_identity"]
    * sequences_with_lengths["v_evaluated_bases"]
    / sequences_with_lengths["v_gene_length"]
)

sequences_with_lengths["j_coverage_adjusted_identity"] = (
    sequences_with_lengths["j_identity"]
    * sequences_with_lengths["j_evaluated_bases"]
    / sequences_with_lengths["j_gene_length"]
)


# ============================================================
# FINAL FILTER
#
# Retain sequences with coverage-adjusted identity >= 80%
# for both V and J genes
# ============================================================

filtered_sequences = sequences_with_lengths[
    (
        sequences_with_lengths["v_coverage_adjusted_identity"]
        >= COVERAGE_ADJUSTED_IDENTITY_CUTOFF
    )
    &
    (
        sequences_with_lengths["j_coverage_adjusted_identity"]
        >= COVERAGE_ADJUSTED_IDENTITY_CUTOFF
    )
].copy()

print(
    "After coverage-adjusted identity filtering "
    f"(>= {COVERAGE_ADJUSTED_IDENTITY_CUTOFF}%): "
    f"{len(filtered_sequences)}"
)


# ============================================================
# SAVE OUTPUT
# ============================================================

filtered_sequences.to_csv(
    OUTPUT_FILE,
    sep="\t",
    index=False
)

print(f"Output saved to: {OUTPUT_FILE}")


# ============================================================
# FINAL SUMMARY
# ============================================================

print("\nSummary:")

print(
    f"Total input sequences: "
    f"{len(sequences)}"
)

print(
    f"Productive sequences with complete FR/CDR annotation: "
    f"{len(productive_complete)}"
)

print(
    f"Sequences with V and J coverage-adjusted identity "
    f">= {COVERAGE_ADJUSTED_IDENTITY_CUTOFF}%: "
    f"{len(filtered_sequences)}"
)