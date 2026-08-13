#!/bin/bash

# ============================================================
# Paired-end immunoglobulin repertoire preprocessing with pRESTO
# ============================================================
#
# This script:
# 1. Assembles paired-end reads
# 2. Filters reads by quality (Q30)
# 3. Identifies V and constant-region primers
# 4. Filters sequences by length
# 5. Collapses duplicate sequences
# 6. Separates sequences according to duplicate count
# 7. Converts FASTQ files to FASTA for downstream germline analysis
#
# Requirements:
# pRESTO must be installed and available in the system PATH.
#
# ============================================================


# ============================================================
# INPUT FILES
# ============================================================

SAMPLE="sample_name"

R1="sample_R1.fastq"
R2="sample_R2.fastq"

V_PRIMERS="V_primers.fasta"
C_PRIMERS="C_primers.fasta"


# ============================================================
# 1. ASSEMBLE PAIRED-END READS
# ============================================================

echo "********** Assembling paired-end reads **********"

AssemblePairs.py align \
    -1 "$R1" \
    -2 "$R2" \
    --coord illumina \
    --rc tail \
    --outname "$SAMPLE" \
    --log "${SAMPLE}_AP.log" \
    --failed


# Assembly statistics

ParseLog.py \
    -l "${SAMPLE}_AP.log" \
    -f ID LENGTH OVERLAP ERROR PVALUE


# ============================================================
# 2. QUALITY FILTERING
# ============================================================

echo "********** Filtering reads by quality (Q30) **********"

FilterSeq.py quality \
    -s "${SAMPLE}_assemble-pass.fastq" \
    -q 30 \
    --outname "${SAMPLE}_q30" \
    --log "${SAMPLE}_q30.log" \
    --failed


# Quality statistics

ParseLog.py \
    -l "${SAMPLE}_q30.log" \
    -f ID QUALITY


# ============================================================
# 3. PRIMER IDENTIFICATION
# ============================================================

echo "********** Identifying primers **********"

# Forward / V primer

MaskPrimers.py align \
    -s "${SAMPLE}_q30_quality-pass.fastq" \
    -p "$V_PRIMERS" \
    --mode tag \
    --outname "${SAMPLE}_FWD" \
    --log "${SAMPLE}_FWD.log" \
    --failed \
    --maxerror 0.5


# Reverse / constant-region primer

MaskPrimers.py align \
    -s "${SAMPLE}_FWD_primers-pass.fastq" \
    -p "$C_PRIMERS" \
    --mode tag \
    --outname "${SAMPLE}_REV" \
    --log "${SAMPLE}_REV.log" \
    --failed \
    --revpr \
    --maxerror 0.5


# Primer statistics

ParseLog.py \
    -l "${SAMPLE}_FWD.log" \
    -f ID PRIMER ERROR


# ============================================================
# 4. ORGANIZE PRIMER ANNOTATIONS
# ============================================================

echo "********** Organizing primer annotations **********"

ParseHeaders.py expand \
    -s "${SAMPLE}_REV_primers-pass.fastq" \
    -f PRIMER \
    -o "${SAMPLE}_primers_expanded.fastq"


ParseHeaders.py rename \
    -s "${SAMPLE}_primers_expanded.fastq" \
    -f PRIMER1 PRIMER2 \
    -k VPRIMER CPRIMER \
    -o "${SAMPLE}_primers.fastq"


# ============================================================
# 5. FILTER BY SEQUENCE LENGTH
# ============================================================

echo "********** Filtering reads by length **********"

FilterSeq.py length \
    -s "${SAMPLE}_primers.fastq" \
    -n 100 \
    --outname "${SAMPLE}_trim" \
    --log "${SAMPLE}_trim_length.log" \
    --failed


# ============================================================
# 6. COLLAPSE DUPLICATE SEQUENCES
# ============================================================

echo "********** Collapsing duplicate sequences **********"

CollapseSeq.py \
    -s "${SAMPLE}_trim_length-pass.fastq" \
    -n 10 \
    --inner \
    --uf "$C_PRIMERS" \
    --cf "$V_PRIMERS" \
    --act set \
    --outname "${SAMPLE}_trim"


# ============================================================
# 7. SPLIT SEQUENCES BY DUPLICATE COUNT
# ============================================================

echo "********** Splitting sequences by duplicate count **********"

SplitSeq.py group \
    -s "${SAMPLE}_trim_collapse-unique.fastq" \
    -f DUPCOUNT \
    --num 2 \
    --outname "${SAMPLE}_trim"


# ============================================================
# 8. EXPORT HEADER INFORMATION
# ============================================================

echo "********** Exporting sequence information **********"

ParseHeaders.py table \
    -s "${SAMPLE}_trim_atleast-2.fastq" \
    -f ID DUPCOUNT CPRIMER VPRIMER


# ============================================================
# 9. CONVERT FASTQ TO FASTA
# ============================================================

echo "********** Converting FASTQ files to FASTA **********"

cat "${SAMPLE}_trim_collapse-unique.fastq" \
    | paste - - - - \
    | sed 's/^@/>/g' \
    | cut -f1-2 \
    | tr '\t' '\n' \
    > "${SAMPLE}_collapse-unique.fasta"


cat "${SAMPLE}_trim_under-2.fastq" \
    | paste - - - - \
    | sed 's/^@/>/g' \
    | cut -f1-2 \
    | tr '\t' '\n' \
    > "${SAMPLE}_under-2.fasta"


cat "${SAMPLE}_trim_atleast-2.fastq" \
    | paste - - - - \
    | sed 's/^@/>/g' \
    | cut -f1-2 \
    | tr '\t' '\n' \
    > "${SAMPLE}_atleast-2.fasta"


echo "********** Preprocessing completed **********"