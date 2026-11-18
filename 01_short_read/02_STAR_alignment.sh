
#!/bin/bash
# ==============================================================================
# Script: 02_STAR_alignment.sh
# Purpose: Align paired-end short-read RNA-seq data to the GRCh38 human
#          reference genome using STAR and a pre-built GENCODE v49 index.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Define Directories and Parameters
# ------------------------------------------------------------------------------

PROJECT_DIR="/data/records_wts"
RESOURCES_DIR="/data/Ressources"

GENCODE_VERSION="v49"
SJDB_OVERHANG=100   # 101 bp reads: read length - 1
THREADS=12

# Input: FASTQ files processed using fastp
FASTQ_DIR="${PROJECT_DIR}/02_fastp"

# Output: STAR alignments
OUT_DIR="${PROJECT_DIR}/04_star"

# Pre-built STAR genome index
INDEX_DIR="${RESOURCES_DIR}/STAR/index_gencode_${GENCODE_VERSION}_L${SJDB_OVERHANG}"

mkdir -p "${OUT_DIR}"

echo "=========================================="
echo "STAR alignment - GENCODE ${GENCODE_VERSION}"
echo "Reference genome: GRCh38"
echo "STAR index: ${INDEX_DIR}"
echo "=========================================="

# ------------------------------------------------------------------------------
# 2. Align Paired-End Reads Using STAR
# ------------------------------------------------------------------------------

counter=0

for fn1 in "${FASTQ_DIR}"/*_R1.fp.fastq.gz; do

    # Extract sample ID and identify corresponding R2 file
    id=$(basename "${fn1}" _R1.fp.fastq.gz)
    fn2="${FASTQ_DIR}/${id}_R2.fp.fastq.gz"

    ((counter++))

    echo "=========================================="
    echo "Processing sample: ${id} (n=${counter})"
    echo "=========================================="

    STAR \
        --runThreadN "${THREADS}" \
        --genomeDir "${INDEX_DIR}" \
        --readFilesIn "${fn1}" "${fn2}" \
        --readFilesCommand zcat \
        --outSAMattrRGline ID:"${id}" SM:"${id}" LB:"${id}" PL:illumina \
        --outFileNamePrefix "${OUT_DIR}/${id}." \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts TranscriptomeSAM \
        --outSAMattributes All

    echo "Completed sample: ${id}"

done

# ------------------------------------------------------------------------------
# 3. Pipeline Completion
# ------------------------------------------------------------------------------

echo "=========================================="
echo "STAR alignment completed."
echo "Total samples processed: ${counter}"
echo "Output directory: ${OUT_DIR}"
echo "=========================================="
