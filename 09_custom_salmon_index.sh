#!/bin/bash
# ==============================================================================
# Script: 09_build_custom_salmon_index.sh
# Purpose: Generate a custom decoy-aware Salmon index combining GENCODE v49
#          annotated transcripts and filtered novel transcript candidates.
# ==============================================================================

set -e
set -o pipefail

# ------------------------------------------------------------------------------
# 1. Define Directories and Parameters
# ------------------------------------------------------------------------------

PROJECT_DIR="/data/records_wts"
RESOURCES_DIR="/data/Ressources/gencode"

GENCODE_VERSION="v49"

# Reference annotation and genome
REFERENCE_GTF="${RESOURCES_DIR}/gencode.${GENCODE_VERSION}.primary_assembly.annotation.gtf"
GENOME_FASTA="${RESOURCES_DIR}/GRCh38.primary_assembly.genome.fa"
DECOYS="${RESOURCES_DIR}/GRCh38_decoys.txt"

# Filtered novel transcript annotation
NOVEL_GTF="${PROJECT_DIR}/03_hybrid_assembly/final_annotation/merged_pass4.gtf"

# Output directory
OUT_DIR="${PROJECT_DIR}/04_custom_salmon_index"

COMBINED_GTF="${OUT_DIR}/combined_annotation.gtf"
TRANSCRIPT_FASTA="${OUT_DIR}/combined_transcripts.fa"
GENTROME="${OUT_DIR}/gentrome_GRCh38_v49_seat.fa.gz"
SALMON_INDEX="${OUT_DIR}/salmon_index"

mkdir -p "${OUT_DIR}"

# Check required input files
for file in "${REFERENCE_GTF}" "${NOVEL_GTF}" \
            "${GENOME_FASTA}" "${DECOYS}"; do

    if [ ! -s "${file}" ]; then
        echo "ERROR: Missing or empty input file: ${file}"
        exit 1
    fi

done

# ------------------------------------------------------------------------------
# 2. Combine Reference and Novel Transcript Annotations
# ------------------------------------------------------------------------------

echo "Combining GENCODE and novel transcript annotations..."

cat "${REFERENCE_GTF}" "${NOVEL_GTF}" > "${COMBINED_GTF}"

echo "Combined annotation: ${COMBINED_GTF}"

# ------------------------------------------------------------------------------
# 3. Extract Transcript Sequences Using gffread
# ------------------------------------------------------------------------------

echo "Extracting transcript sequences..."

gffread "${COMBINED_GTF}" \
    -g "${GENOME_FASTA}" \
    -w "${TRANSCRIPT_FASTA}"

echo "Transcript FASTA generated."

# ------------------------------------------------------------------------------
# 4. Build the Decoy-Aware Transcriptome (Gentrome)
# ------------------------------------------------------------------------------

echo "Building custom gentrome..."

# Combine transcript sequences and the GRCh38 genome.
# Compress the resulting gentrome for Salmon indexing.

cat "${TRANSCRIPT_FASTA}" "${GENOME_FASTA}" | \
    gzip -c > "${GENTROME}"

echo "Gentrome generated: ${GENTROME}"

# ------------------------------------------------------------------------------
# 5. Build Custom Salmon Index
# ------------------------------------------------------------------------------

echo "Building Salmon index..."

salmon index \
    -t "${GENTROME}" \
    -d "${DECOYS}" \
    -i "${SALMON_INDEX}" \
    --gencode

echo "Salmon index generated."

# ------------------------------------------------------------------------------
# 6. Pipeline Completion
# ------------------------------------------------------------------------------

echo "=========================================================="
echo "Custom Salmon index generation completed successfully."
echo "GENCODE version:       ${GENCODE_VERSION}"
echo "Combined annotation:   ${COMBINED_GTF}"
echo "Transcript FASTA:      ${TRANSCRIPT_FASTA}"
echo "Gentrome FASTA:        ${GENTROME}"
echo "Decoy sequence list:   ${DECOYS}"
echo "Salmon index:          ${SALMON_INDEX}"
echo "=========================================================="

