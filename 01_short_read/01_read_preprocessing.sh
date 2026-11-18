
#!/usr/bin/env bash
# ==============================================================================
# Script: 01_read_preprocessing.sh
# Purpose: Download and verify raw paired-end FASTQ files, merge sequencing
#          lanes, validate read pairing, and perform preprocessing with fastp.
#
# Input:   Illumina paired-end FASTQ files (lanes L001-L004)
# Output:  Merged FASTQ files and fastp-processed reads with QC reports
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# 0. Define Project Directories and Parameters
# ------------------------------------------------------------------------------

PROJECT_DIR="/data/records_wts"

DOWNLOAD_DIR="${PROJECT_DIR}/00_download"
FASTQ_DIR="${PROJECT_DIR}/01_fastq"
FASTP_DIR="${PROJECT_DIR}/02_fastp"

# One download URL per line
URL_LIST="${PROJECT_DIR}/fastq_urls.txt"

# Original MD5 checksum manifest supplied by the sequencing provider
MD5_FILE="${DOWNLOAD_DIR}/MD5SUMS"

# Expected Illumina sequencing lanes
LANES=(001 002 003 004)

# Preserve the original fastp configuration
THREADS=12

mkdir -p "${DOWNLOAD_DIR}" "${FASTQ_DIR}" "${FASTP_DIR}"

# Check required software
for tool in wget md5sum gzip fastp awk cmp sort; do
    command -v "${tool}" >/dev/null 2>&1 || {
        echo "ERROR: ${tool} is not installed."
        exit 1
    }
done

# Temporary files are removed when the script exits
TMP_DIR=$(mktemp -d "${PROJECT_DIR}/.readprep.XXXXXX")
trap 'rm -rf -- "${TMP_DIR}"' EXIT

echo "Starting short-read RNA-seq preprocessing..."

# ------------------------------------------------------------------------------
# 1. Download Raw FASTQ Files
# ------------------------------------------------------------------------------

echo "Step 1: Downloading raw FASTQ files..."

# The URL list must contain the actual provider download links.
# wget -c resumes interrupted downloads when supported by the server.

if [[ ! -s "${URL_LIST}" ]]; then
    echo "ERROR: Missing FASTQ download URL list: ${URL_LIST}"
    exit 1
fi

wget -c \
     -i "${URL_LIST}" \
     -P "${DOWNLOAD_DIR}"

echo "Download completed."

# ------------------------------------------------------------------------------
# 2. Verify File Integrity Using MD5 Checksums
# ------------------------------------------------------------------------------

echo "Step 2: Verifying MD5 checksums..."

# The original MD5 manifest must be provided separately.
# Each entry must correspond to a downloaded file.

if [[ ! -s "${MD5_FILE}" ]]; then
    echo "ERROR: MD5 checksum manifest not found."
    exit 1
fi

(
    cd "${DOWNLOAD_DIR}"
    md5sum --check --strict "${MD5_FILE}" \
        > "${DOWNLOAD_DIR}/md5_verification.log"
) || {
    echo "ERROR: MD5 verification failed."
    echo "See ${DOWNLOAD_DIR}/md5_verification.log"
    exit 1
}

echo "MD5 verification successful."

# ------------------------------------------------------------------------------
# 3. Identify Samples and Verify Sequencing Lanes
# ------------------------------------------------------------------------------

echo "Step 3: Checking sample files and sequencing lanes..."

# Identify samples without truncating names containing underscores.
# Expected naming:
# SAMPLE_L001_R1_001.fastq.gz
# SAMPLE_L001_R2_001.fastq.gz

shopt -s nullglob

INPUT_FILES=("${DOWNLOAD_DIR}"/*_L???_R?_001.fastq.gz)

if (( ${#INPUT_FILES[@]} == 0 )); then
    echo "ERROR: No paired-end FASTQ files found."
    exit 1
fi

declare -A SAMPLE_IDS

for file in "${INPUT_FILES[@]}"; do

    filename="${file##*/}"

    if [[ "${filename}" =~ ^(.+)_L([0-9]{3})_R([12])_001\.fastq\.gz$ ]]; then
        sample="${BASH_REMATCH[1]}"
        SAMPLE_IDS["${sample}"]=1
    else
        echo "ERROR: Unexpected FASTQ filename: ${filename}"
        exit 1
    fi

    # Confirm that every input FASTQ is included in the MD5 manifest.
    if ! awk -v target="${filename}" '
        {
            name = substr($0, 35)
            sub(/^\.\//, "", name)
            if (name == target) found = 1
        }
        END { exit !found }
    ' "${MD5_FILE}"; then
        echo "ERROR: No MD5 entry found for ${filename}"
        exit 1
    fi
done

mapfile -t SAMPLES < <(
    printf '%s\n' "${!SAMPLE_IDS[@]}" | LC_ALL=C sort
)

echo "Identified ${#SAMPLES[@]} samples."

# Check that every sample contains exactly the expected R1/R2 lanes.
for sample in "${SAMPLES[@]}"; do

    for lane in "${LANES[@]}"; do
        for read in 1 2; do

            file="${DOWNLOAD_DIR}/${sample}_L${lane}_R${read}_001.fastq.gz"

            if [[ ! -f "${file}" ]]; then
                echo "ERROR: Missing file: ${file}"
                exit 1
            fi

        done
    done

    # Reject unexpected or additional lanes.
    for file in "${DOWNLOAD_DIR}/${sample}"_L???_R?_001.fastq.gz; do

        name="${file##*/}"

        if [[ "${name}" =~ ^(.+)_L([0-9]{3})_R([12])_001\.fastq\.gz$ ]]; then
            lane="${BASH_REMATCH[2]}"
        fi

        valid=false
        for expected in "${LANES[@]}"; do
            if [[ "${lane}" == "${expected}" ]]; then
                valid=true
                break
            fi
        done

        if [[ "${valid}" == false ]]; then
            echo "ERROR: Unexpected lane ${lane} for ${sample}"
            exit 1
        fi
    done

    echo "All expected lanes found for ${sample}"

done

# ------------------------------------------------------------------------------
# 4. Merge FASTQ Files Across Sequencing Lanes
# ------------------------------------------------------------------------------

echo "Step 4: Merging sequencing lanes..."

for sample in "${SAMPLES[@]}"; do

    echo "Merging sample: ${sample}"

    R1_FILES=()
    R2_FILES=()

    for lane in "${LANES[@]}"; do

        R1="${DOWNLOAD_DIR}/${sample}_L${lane}_R1_001.fastq.gz"
        R2="${DOWNLOAD_DIR}/${sample}_L${lane}_R2_001.fastq.gz"

        # Verify compressed file integrity before merging.
        gzip -t "${R1}" "${R2}"

        R1_FILES+=("${R1}")
        R2_FILES+=("${R2}")

    done

    # Concatenation preserves lane order and gzip compression.
    cat "${R1_FILES[@]}" \
        > "${TMP_DIR}/${sample}_R1.fastq.gz"

    cat "${R2_FILES[@]}" \
        > "${TMP_DIR}/${sample}_R2.fastq.gz"

    echo "Merged ${#LANES[@]} lanes for ${sample}"

done

# ------------------------------------------------------------------------------
# 5. Verify R1/R2 Read Concordance
# ------------------------------------------------------------------------------

echo "Step 5: Checking paired-end read concordance..."

# Extract read identifiers and validate standard four-line FASTQ records.
# Unlike filtering on a specific sequencer prefix, all reads are checked.

extract_read_ids() {

    local input_file="$1"
    local output_file="$2"

    gzip -cd "${input_file}" | awk '

        NR % 4 == 1 {
            if ($0 !~ /^@/) exit 1

            id = $1
            sub(/^@/, "", id)
            sub(/\/[12]$/, "", id)

            print id
        }

        NR % 4 == 2 {
            sequence_length = length($0)
        }

        NR % 4 == 3 {
            if ($0 !~ /^\+/) exit 1
        }

        NR % 4 == 0 {
            if (length($0) != sequence_length) exit 1
        }

        END {
            if (NR == 0 || NR % 4 != 0) exit 1
        }

    ' > "${output_file}"
}

for sample in "${SAMPLES[@]}"; do

    echo "Checking sample: ${sample}"

    R1="${TMP_DIR}/${sample}_R1.fastq.gz"
    R2="${TMP_DIR}/${sample}_R2.fastq.gz"

    ID_R1="${TMP_DIR}/${sample}.R1.ids"
    ID_R2="${TMP_DIR}/${sample}.R2.ids"

    extract_read_ids "${R1}" "${ID_R1}"
    extract_read_ids "${R2}" "${ID_R2}"

    if ! cmp -s "${ID_R1}" "${ID_R2}"; then
        echo "ERROR: R1/R2 read identifiers do not match: ${sample}"
        exit 1
    fi

    echo "R1/R2 concordance verified: ${sample}"

    # Save validated merged FASTQ files.
    mv "${R1}" "${FASTQ_DIR}/${sample}_R1.fastq.gz"
    mv "${R2}" "${FASTQ_DIR}/${sample}_R2.fastq.gz"

    rm -f "${ID_R1}" "${ID_R2}"

done

# ------------------------------------------------------------------------------
# 6. Adapter Trimming and Quality Filtering Using fastp
# ------------------------------------------------------------------------------

echo "Step 6: Running fastp..."

for sample in "${SAMPLES[@]}"; do

    echo "Processing sample: ${sample}"

    R1="${FASTQ_DIR}/${sample}_R1.fastq.gz"
    R2="${FASTQ_DIR}/${sample}_R2.fastq.gz"

    fastp \
        -i "${R1}" \
        -I "${R2}" \
        -o "${FASTP_DIR}/${sample}_R1.fp.fastq.gz" \
        -O "${FASTP_DIR}/${sample}_R2.fp.fastq.gz" \
        -w "${THREADS}" \
        -j "${FASTP_DIR}/${sample}.fp.json" \
        -h "${FASTP_DIR}/${sample}.fp.html" \
        -R "${sample} fastp report"

    echo "fastp completed: ${sample}"

done

# ------------------------------------------------------------------------------
# 7. Pipeline Completion
# ------------------------------------------------------------------------------

echo "=========================================================="
echo "Short-read RNA-seq preprocessing completed successfully."
echo "Number of samples: ${#SAMPLES[@]}"
echo "Raw FASTQ files:    ${DOWNLOAD_DIR}"
echo "Merged FASTQ files: ${FASTQ_DIR}"
echo "Processed FASTQ:    ${FASTP_DIR}"
echo "=========================================================="
