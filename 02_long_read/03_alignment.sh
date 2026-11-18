#!/bin/bash
# ==============================================================================
# Script: 03_minimap2_alignment.sh
# Purpose: Align Chopper-filtered ONT long reads to the human reference genome
#          using minimap2 and convert SAM alignments to BAM using SAMtools.
# ==============================================================================

set -o pipefail

# ------------------------------------------------------------------------------
# 1. Define Directories and Parameters
# ------------------------------------------------------------------------------

PROJECT_DIR="/data/records/ONT/long_reads"

FASTQ_DIR="${PROJECT_DIR}/03_chopper"
OUT_DIR="${PROJECT_DIR}/04_minimap2"

# Pre-built minimap2 genome index
GENOME_INDEX="/data/Ressources/GRCh38.p14/GRCh38.mmi"

THREADS=12

mkdir -p "${OUT_DIR}"

# ------------------------------------------------------------------------------
# 2. Align Long-Read RNA-seq Data Using minimap2
# ------------------------------------------------------------------------------

counter=0

for fn in "${FASTQ_DIR}"/*.chop.fastq.gz; do

    # Extract sample identifier
    id=$(basename "${fn}" .chop.fastq.gz)

    bam_out="${OUT_DIR}/${id}.bam"

    ((counter+=1))

    echo "=========================================="
    echo "Processing sample: ${id} (n=${counter})"
    echo "=========================================="

    minimap2 \
        -ax splice \
        -t "${THREADS}" \
        "${GENOME_INDEX}" "${fn}" | \
    samtools view \
        -@ "${THREADS}" \
        -Sb \
        -o "${bam_out}" -

    if [ $? -ne 0 ]; then
        echo "ERROR: Alignment failed for ${id}"
        exit 1
    fi

    echo "Alignment completed: ${bam_out}"

done

# ------------------------------------------------------------------------------
# 3. Pipeline Completion
# ------------------------------------------------------------------------------

echo "=========================================="
echo "Long-read RNA-seq alignment completed."
echo "Total samples processed: ${counter}"
echo "Output BAM files: ${OUT_DIR}"
echo "=========================================="
