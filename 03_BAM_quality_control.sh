
#!/bin/bash
# ==============================================================================
# Script: 03_BAM_quality_control.sh
# Purpose: Assess the quality of STAR-aligned RNA-seq reads using SAMtools
#          and Qualimap, including alignment statistics, BAM quality
#          metrics, and RNA-seq-specific quality assessment.
# ==============================================================================

set -e

# ------------------------------------------------------------------------------
# 1. Define Directories and Parameters
# ------------------------------------------------------------------------------

PROJECT_DIR="/data/records_wts"
RESOURCES_DIR="/data/Ressources"

GENCODE_VERSION="v49"

BAM_DIR="${PROJECT_DIR}/04_star"
QC_DIR="${PROJECT_DIR}/05_qc_bam"

GTF="${RESOURCES_DIR}/gencode/gencode.${GENCODE_VERSION}.primary_assembly.annotation.gtf"

# Create output directories
mkdir -p "${QC_DIR}/samtools_bamstats"
mkdir -p "${QC_DIR}/qualimap_bamqc"
mkdir -p "${QC_DIR}/qualimap_rnaseq"

# Run Qualimap without graphical display
unset DISPLAY

echo "=========================================="
echo "BAM Quality Control - GENCODE ${GENCODE_VERSION}"
echo "Input BAM directory: ${BAM_DIR}"
echo "Reference annotation: ${GTF}"
echo "=========================================="

# ------------------------------------------------------------------------------
# 2. Alignment Statistics Using SAMtools
# ------------------------------------------------------------------------------

echo "Running SAMtools alignment statistics..."

for bam in "${BAM_DIR}"/*.Aligned.sortedByCoord.out.bam; do

    id=$(basename "${bam}" .Aligned.sortedByCoord.out.bam)

    echo "Processing sample: ${id}"

    # Generate detailed alignment statistics
    samtools stats "${bam}" \
        > "${QC_DIR}/samtools_bamstats/${id}.samtools.stats"

    # Summarize mapping and alignment flag statistics
    samtools flagstat "${bam}" \
        > "${QC_DIR}/samtools_bamstats/${id}.samtools.flagstats"

done

# ------------------------------------------------------------------------------
# 3. BAM Quality Assessment Using Qualimap
# ------------------------------------------------------------------------------

echo "Running Qualimap BamQC..."

for bam in "${BAM_DIR}"/*.Aligned.sortedByCoord.out.bam; do

    id=$(basename "${bam}" .Aligned.sortedByCoord.out.bam)

    echo "Processing sample: ${id}"

    qualimap bamqc \
        -bam "${bam}" \
        -gff "${GTF}" \
        -c \
        --genome-gc-distr HUMAN \
        --outside-stats \
        -outdir "${QC_DIR}/qualimap_bamqc/${id}" \
        -p strand-specific-reverse

done

# ------------------------------------------------------------------------------
# 4. RNA-seq Quality Assessment Using Qualimap
# ------------------------------------------------------------------------------

echo "Running Qualimap RNA-seq QC..."

for bam in "${BAM_DIR}"/*.Aligned.sortedByCoord.out.bam; do

    id=$(basename "${bam}" .Aligned.sortedByCoord.out.bam)

    echo "Processing sample: ${id}"

    qualimap rnaseq \
        -bam "${bam}" \
        -gtf "${GTF}" \
        -oc "${QC_DIR}/qualimap_rnaseq/${id}.counts" \
        -outdir "${QC_DIR}/qualimap_rnaseq/${id}" \
        -p strand-specific-reverse \
        --paired \
        --sorted

done

# ------------------------------------------------------------------------------
# 5. Pipeline Completion
# ------------------------------------------------------------------------------

echo "=========================================="
echo "BAM quality control completed successfully."
echo "Input BAM files:     ${BAM_DIR}"
echo "SAMtools reports:    ${QC_DIR}/samtools_bamstats"
echo "Qualimap BamQC:      ${QC_DIR}/qualimap_bamqc"
echo "Qualimap RNA-seq QC: ${QC_DIR}/qualimap_rnaseq"
echo "=========================================="
