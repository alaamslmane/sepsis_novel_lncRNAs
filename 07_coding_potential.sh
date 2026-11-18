#!/bin/bash
# ==============================================================================
# Script: 07_coding_potential.sh
# Purpose: Extract transcript sequences and assess coding potential using
#          CNCI, CPC2, and CPAT through ezLncPred.
# ==============================================================================

set -e

# ------------------------------------------------------------------------------
# 1. Define Directories
# ------------------------------------------------------------------------------

PROJECT_DIR="/data/records_wts/03_hybrid_assembly"

INPUT_GTF="${PROJECT_DIR}/annotation/merged_pass2.gtf"
OUT_DIR="${PROJECT_DIR}/coding_potential"

GENOME_FASTA="/data/Ressources/gencode/GRCh38.primary_assembly.genome.fa"
TRANSCRIPT_FASTA="${OUT_DIR}/transcripts.fa"

mkdir -p "${OUT_DIR}"

# ------------------------------------------------------------------------------
# 2. Extract Transcript Sequences Using gffread
# ------------------------------------------------------------------------------

echo "Extracting transcript sequences..."

gffread "${INPUT_GTF}" \
    -g "${GENOME_FASTA}" \
    -w "${TRANSCRIPT_FASTA}"

# ------------------------------------------------------------------------------
# 3. Assess Coding Potential Using CNCI
# ------------------------------------------------------------------------------

echo "Running CNCI..."

mkdir -p "${OUT_DIR}/cnci"

ezLncPred \
    -i "${TRANSCRIPT_FASTA}" \
    -o "${OUT_DIR}/cnci/results" \
    CNCI \
    -p ve \
    --parallel

# ------------------------------------------------------------------------------
# 4. Assess Coding Potential Using CPC2
# ------------------------------------------------------------------------------

echo "Running CPC2..."

mkdir -p "${OUT_DIR}/cpc2"

ezLncPred \
    -i "${TRANSCRIPT_FASTA}" \
    -o "${OUT_DIR}/cpc2/results" \
    CPC2

# ------------------------------------------------------------------------------
# 5. Assess Coding Potential Using CPAT
# ------------------------------------------------------------------------------

echo "Running CPAT..."

mkdir -p "${OUT_DIR}/cpat"

ezLncPred \
    -i "${TRANSCRIPT_FASTA}" \
    -o "${OUT_DIR}/cpat/results" \
    CPAT \
    -p Human

# ------------------------------------------------------------------------------
# 6. Pipeline Completion
# ------------------------------------------------------------------------------

echo "=========================================="
echo "Coding potential assessment completed."
echo "Transcript FASTA: ${TRANSCRIPT_FASTA}"
echo "CNCI results: ${OUT_DIR}/cnci"
echo "CPC2 results: ${OUT_DIR}/cpc2"
echo "CPAT results: ${OUT_DIR}/cpat"
echo "=========================================="

