#!/bin/bash
# ==============================================================================
# Script: 02_chopper_filtering.sh
# Purpose: Filter Oxford Nanopore long-read RNA-seq data based on read
#          quality and length using Chopper.
# ==============================================================================

set -o pipefail

# ------------------------------------------------------------------------------
# 1. Define Directories and Parameters
# ------------------------------------------------------------------------------

PROJECT_DIR="/data/records/ONT/long_reads"

FASTQ_DIR="${PROJECT_DIR}/01_fastq"
OUT_DIR="${PROJECT_DIR}/02_chopper"

MIN_QUALITY=10
MIN_LENGTH=200
THREADS=12

mkdir -p "${OUT_DIR}"

# ------------------------------------------------------------------------------
# 2. Filter Long-Read RNA-seq Data Using Chopper
# ------------------------------------------------------------------------------

counter=0

for fn in "${FASTQ_DIR}"/*.fastq.gz; do

    # Extract sample identifier
    id=$(basename "${fn}" .fastq.gz)

    output="${OUT_DIR}/${id}.chop.fastq.gz"

    ((counter+=1))

    echo "=========================================="
    echo "Processing sample: ${id} (n=${counter})"
    echo "=========================================="

    chopper \
        -q "${MIN_QUALITY}" \
        -l "${MIN_LENGTH}" \
        -i "${fn}" \
        -t "${THREADS}" | \
    gzip > "${output}"

    if [ $? -ne 0 ]; then
        echo "ERROR: Chopper filtering failed for ${id}"
        exit 1
    fi

    echo "Filtering completed: ${output}"

done

# ------------------------------------------------------------------------------
# 3. Pipeline Completion
# ------------------------------------------------------------------------------

echo "=========================================="
echo "Long-read RNA-seq filtering completed."
echo "Total samples processed: ${counter}"
echo "Filtered FASTQ files: ${OUT_DIR}"
echo "=========================================="

