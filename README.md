# RNA-seq Pipeline for Novel long noncoding RNA (lncRNAs) Identification and Quantification in sepsis

This repository contains the Linux-based bioinformatics pipeline developed to identify and quantify novel long noncoding sepsis-associated transcripts (SEATs) by integrating long-read and short-read RNA sequencing data. It combines long-read RNA-seq (n = 10 samples) and short-read RNA-seq (n = 10 samples) for hybrid transcriptome assembly. Both sequencing datasets are first processed and aligned independently before being jointly used for transcript reconstruction. The resulting transcript assemblies are merged to generate a unified transcriptome annotation, including previously unannotated transcripts. This annotation then undergoes a series of filtering steps to identify and retain high-confidence novel long noncoding SEATs candidates. The resulting custom annotation is subsequently used to quantify transcript expression across n = 812 additional short-read RNA-seq samples, enabling the characterisation of newly identified transcripts in a larger cohort.

## Pipeline Architecture
The pipeline is organized into three main phases:

**Phase 1: RNA-seq Pre-processing and Alignment**

**A. Short-read RNA-seq Processing**

- **`01_read_preprocessing.sh`**: Downloads raw paired-end FASTQ files, verifies file integrity using MD5 checksums, and merges sequencing lanes for each sample. Checks R1/R2 read concordance and performs adapter trimming and quality filtering using `fastp`. Sequencing quality is assessed using `FastQC` and `FastQ Screen`, with QC metrics consolidated into a single MultiQC report.
- **`02_Alignment.sh`**: Aligns fastp-processed paired-end reads to the human reference genome (GRCh38) using `STAR`, with GENCODE v49 annotations.
- **`03_BAM_quality_control.sh`**: Assesses STAR-aligned BAM files using SAMtools and Qualimap. Generates alignment statistics, evaluates BAM quality metrics, and performs RNA-seq-specific quality assessment using GENCODE v49 annotations and a reverse-stranded paired-end library configuration.

  
**B. Long-read RNA-seq Processing**

- **`01_read_preprocessing.sh`**: Downloads raw FASTQ files, verifies file integrity using MD5 checksums, and assesses read length distributions and sequencing quality using `NanoPlot`. Performs read trimming and quality filtering using `Chopper`.
- **`02_Alignment.sh`**: Aligns quality-filtered long RNA-seq reads to the human reference genome (GRCh38) using `minimap2`, generating genome-aligned reads for downstream transcriptome assembly.
- **`03_BAM_quality_control.sh`**: Evaluates long-read alignment quality using `Cramino`, generating alignment statistics and quality metrics.
