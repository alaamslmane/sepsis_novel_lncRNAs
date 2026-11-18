#!/bin/bash
# ==============================================================================
# Script: 04_cramino_QC.sh
# Purpose: Assess the alignment quality of Oxford Nanopore long-read
#          RNA-seq data using Cramino.
# ==============================================================================

set -e

# ------------------------------------------------------------------------------
# 1. Define Directories and Parameters
# ------------------------------------------------------------------------------

PROJECT_DIR="/data/records/ONT/long_reads"

BAM_DIR="${PROJECT_DIR}/03_minimap2"
OUT_DIR="${PROJECT_DIR}/04_cramino"

THREADS=12

mkdir -p "${OUT_DIR}"

# ------------------------------------------------------------------------------
# 2. Assess Long-Read Alignment Quality Using Cramino
# ------------------------------------------------------------------------------

counter=0

for fn in "${BAM_DIR}"/*.bam; do

    # Extract sample identifier
    id=$(basename "${fn}" .bam)

    ((counter+=1))

    echo "=========================================="
    echo "Processing sample: ${id} (n=${counter})"
    echo "=========================================="

    cramino "${fn}" \
        --threads "${THREADS}" \
        --hist \
        > "${OUT_DIR}/${id}.cram.txt"

    echo "Cramino completed: ${id}"

done

# ------------------------------------------------------------------------------
# 3. Pipeline Completion
# ------------------------------------------------------------------------------

echo "=========================================="
echo "Long-read alignment quality control completed."
echo "Total samples processed: ${counter}"
echo "Cramino reports: ${OUT_DIR}"
echo "=========================================="

