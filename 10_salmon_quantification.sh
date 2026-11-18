#!/bin/bash
# ==============================================================================
# Script: 10_salmon_quantification.sh
# Purpose: Quantify transcript abundance in paired-end short-read RNA-seq
#          samples using Salmon and the custom decoy-aware transcriptome index.
# ==============================================================================

set -e

# ------------------------------------------------------------------------------
# 1. Define Directories and Parameters
# ------------------------------------------------------------------------------

PROJECT_DIR="/data/records_wts"

# Input: fastp-processed paired-end reads
FASTQ_DIR="${PROJECT_DIR}/02_fastp"

# Custom Salmon index generated from GENCODE v49 and novel transcripts
SALMON_INDEX="${PROJECT_DIR}/04_custom_salmon_index/salmon_index"

# Output: sample-level Salmon quantification results
OUT_DIR="${PROJECT_DIR}/06_salmon"

THREADS=12
BOOTSTRAPS=100

mkdir -p "${OUT_DIR}"

echo "=========================================="
echo "Salmon Transcript Quantification"
echo "Reference: Custom GENCODE v49 transcriptome"
echo "Index: ${SALMON_INDEX}"
echo "=========================================="

# ------------------------------------------------------------------------------
# 2. Quantify Transcript Expression Using Salmon
# ------------------------------------------------------------------------------

counter=0

for fn1 in "${FASTQ_DIR}"/*_R1.fp.fastq.gz; do

    # Extract sample identifier
    id=$(basename "${fn1}" _R1.fp.fastq.gz)

    # Corresponding paired-end read file
    fn2="${FASTQ_DIR}/${id}_R2.fp.fastq.gz"

    if [ ! -f "${fn2}" ]; then
        echo "ERROR: Missing R2 file for ${id}"
        exit 1
    fi

    ((counter+=1))

    echo "=========================================="
    echo "Processing sample: ${id} (n=${counter})"
    echo "=========================================="

    salmon quant \
        -i "${SALMON_INDEX}" \
        -l A \
        -1 "${fn1}" \
        -2 "${fn2}" \
        --numBootstraps "${BOOTSTRAPS}" \
        --validateMappings \
        -p "${THREADS}" \
        -o "${OUT_DIR}/${id}"

    echo "Quantification completed: ${id}"

done

# ------------------------------------------------------------------------------
# 3. Pipeline Completion
# ------------------------------------------------------------------------------

echo "=========================================="
echo "Salmon quantification completed successfully."
echo "Total samples processed: ${counter}"
echo "Custom Salmon index: ${SALMON_INDEX}"
echo "Quantification results: ${OUT_DIR}"
echo "=========================================="
