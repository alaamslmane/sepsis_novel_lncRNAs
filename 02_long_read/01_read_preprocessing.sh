
#!/bin/bash
# ==============================================================================
# Script: 01_long_read_preprocessing.sh
# Purpose: Download and verify ONT FASTQ files, merge two sequencing runs
#          per sample, and assess read quality using NanoPlot.
# ==============================================================================

set -e

# ------------------------------------------------------------------------------
# 1. Define Directories and Parameters
# ------------------------------------------------------------------------------

PROJECT_DIR="/data/records/ONT/long_reads"

DOWNLOAD_DIR="${PROJECT_DIR}/00_download"
FASTQ_DIR="${PROJECT_DIR}/01_fastq"
QC_DIR="${PROJECT_DIR}/02_nanoplot"

URL_LIST="${PROJECT_DIR}/fastq_urls.txt"
SAMPLE_LIST="${PROJECT_DIR}/samples_runs.tsv"
MD5_FILE="${DOWNLOAD_DIR}/MD5SUMS"

mkdir -p "${DOWNLOAD_DIR}" "${FASTQ_DIR}" "${QC_DIR}"

# ------------------------------------------------------------------------------
# 2. Download Raw Sequencing Files
# ------------------------------------------------------------------------------

echo "Downloading ONT sequencing files..."

# Each line of fastq_urls.txt contains a direct download URL.
wget -c -i "${URL_LIST}" -P "${DOWNLOAD_DIR}"

echo "Download completed."

# ------------------------------------------------------------------------------
# 3. Verify File Integrity
# ------------------------------------------------------------------------------

echo "Checking MD5 checksums..."

(
    cd "${DOWNLOAD_DIR}"
    md5sum -c "${MD5_FILE}"
)

echo "MD5 verification completed."

# ------------------------------------------------------------------------------
# 4. Merge the Two Sequencing Runs per Sample
# ------------------------------------------------------------------------------

echo "Merging sequencing runs..."

# samples_runs.tsv contains three tab-separated columns:
# sample_id, run1_file, run2_file
# No header line is expected.
#
# Run files can be .fastq or .fastq.gz.
# The original files are preserved in 00_download.

read_fastq() {
    case "$1" in
        *.fastq.gz) gzip -cd "$1" ;;
        *.fastq)    cat "$1" ;;
        *) echo "Unsupported FASTQ format: $1" >&2; return 1 ;;
    esac
}

while IFS=$'\t' read -r sample run1 run2; do

    [[ -z "${sample}" ]] && continue

    FILE1="${DOWNLOAD_DIR}/${run1}"
    FILE2="${DOWNLOAD_DIR}/${run2}"
    OUTPUT="${FASTQ_DIR}/${sample}.fastq.gz"

    echo "Merging runs for sample: ${sample}"

    if [[ ! -f "${FILE1}" || ! -f "${FILE2}" ]]; then
        echo "ERROR: Missing sequencing run for ${sample}"
        exit 1
    fi

    # Decompress each run and combine reads into one FASTQ.gz file.
    (
        read_fastq "${FILE1}"
        read_fastq "${FILE2}"
    ) | gzip -c > "${OUTPUT}"

    gzip -t "${OUTPUT}"

    echo "Merged FASTQ: ${OUTPUT}"

done < "${SAMPLE_LIST}"

# ------------------------------------------------------------------------------
# 5. Raw Long-Read Quality Assessment Using NanoPlot
# ------------------------------------------------------------------------------

echo "Starting NanoPlot quality control..."

for fn in "${FASTQ_DIR}"/*.fastq.gz; do

    id=$(basename "${fn}" .fastq.gz)

    echo "Processing sample: ${id}"

    NanoPlot \
        --fastq "${fn}" \
        -o "${QC_DIR}/${id}_NanoPlot" \
        --downsample 10000 \
        --loglength \
        -t 12 \
        --N50

    echo "NanoPlot completed: ${id}"

done

# ------------------------------------------------------------------------------
# 6. Pipeline Completion
# ------------------------------------------------------------------------------

echo "=========================================================="
echo "Long-read RNA-seq preprocessing completed."
echo "Raw sequencing files: ${DOWNLOAD_DIR}"
echo "Merged FASTQ files:   ${FASTQ_DIR}"
echo "NanoPlot reports:     ${QC_DIR}"
echo "=========================================================="

