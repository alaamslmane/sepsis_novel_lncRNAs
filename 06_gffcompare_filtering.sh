#!/bin/bash
# ==============================================================================
# Script: 06_gffcompare_filtering.sh
# Purpose: Classify assembled transcripts against GENCODE v49 and retain
#          transcripts with gffcompare class codes u, i, x and  y
# ==============================================================================

set -e

# ------------------------------------------------------------------------------
# 1. Define Directories and Parameters
# ------------------------------------------------------------------------------

PROJECT_DIR="/data/records_wts/03_hybrid_assembly"

INPUT_GTF="${PROJECT_DIR}/stringtie/merged_pass1.gtf"
OUT_DIR="${PROJECT_DIR}/annotation"

REFERENCE_GTF="/data/Ressources/gencode/gencode.v49.primary_assembly.annotation.gtf"

mkdir -p "${OUT_DIR}"

# ------------------------------------------------------------------------------
# 2. Compare Assembled Transcripts with GENCODE v49
# ------------------------------------------------------------------------------

echo "Running gffcompare..."

gffcompare \
    -r "${REFERENCE_GTF}" \
    -o "${OUT_DIR}/compare" \
    "${INPUT_GTF}"

# ------------------------------------------------------------------------------
# 3. Filter Transcripts by Genomic Class
# ------------------------------------------------------------------------------

echo "Selecting transcripts classified as u, i, x and y"

ANNOTATED_GTF="${OUT_DIR}/compare.annotated.gtf"
FILTERED_GTF="${OUT_DIR}/merged_pass2.gtf"

# Extract a quoted attribute from GTF column 9.
# First collect transcript IDs belonging to the selected classes.
# Then retain all GTF records for those transcript IDs.

awk '
function attribute(text, key,   n, fields, i, field) {
    n = split(text, fields, ";")
    for (i = 1; i <= n; i++) {
        field = fields[i]
        sub(/^[[:space:]]+/, "", field)
        if (index(field, key " \"") == 1) {
            sub(/^[^"]*"/, "", field)
            sub(/".*$/, "", field)
            return field
        }
    }
    return ""
}

FNR == NR {
    if ($3 == "transcript") {
        tid = attribute($9, "transcript_id")
        cls = attribute($9, "class_code")

        if (tid != "" && (cls == "u" || cls == "i" || cls == "x" || cls == "y"))
            keep[tid] = 1
    }
    next
}

{
    tid = attribute($9, "transcript_id")
    if (tid in keep)
        print
}
' "${ANNOTATED_GTF}" "${ANNOTATED_GTF}" > "${FILTERED_GTF}"

# ------------------------------------------------------------------------------
# 4. Summarize Selected Transcripts
# ------------------------------------------------------------------------------

echo "=========================================="
echo "gffcompare classification completed."
echo "Selected classes: u, i, x", "y"
echo "Filtered annotation: ${FILTERED_GTF}"
echo "Number of retained transcripts:"
awk '$3 == "transcript"' "${FILTERED_GTF}" | wc -l
echo "=========================================="

