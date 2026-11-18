# RNA-seq Pipeline for Novel long noncoding RNA (lncRNAs) Identification and Quantification in sepsis

This repository contains the Linux-based bioinformatics pipeline developed to identify and quantify novel long noncoding sepsis-associated transcripts (SEATs) by integrating long-read and short-read RNA sequencing data. It combines long-read RNA-seq (n = 10 samples) and short-read RNA-seq (n = 10 samples) for hybrid transcriptome assembly. Both sequencing datasets are first processed and aligned independently before being jointly used for transcript reconstruction. The resulting transcript assemblies are merged to generate a unified transcriptome annotation, including previously unannotated transcripts. This annotation then undergoes a series of filtering steps to identify and retain high-confidence novel long noncoding SEATs candidates. The resulting custom annotation is subsequently used to quantify transcript expression across n = 812 additional short-read RNA-seq samples, enabling the characterisation of newly identified transcripts in a larger cohort.

## Pipeline Architecture
The pipeline is organized into three main phases:

**Phase 1: RNA-seq Pre-processing and Alignment**

**A. Short-read RNA-seq Processing**

- **`01_read_preprocessing.sh`**: Downloads raw FASTQ files, verifies MD5 checksums, and performs adapter trimming and quality filtering using `fastp`. Sequencing quality is assessed using `FastQC` and `FastQ Screen`, with reports summarized by MultiQC.
- **`02_STAR_alignment.sh`**: Aligns fastp-processed paired-end reads to the human reference genome (GRCh38) using `STAR`, with GENCODE v49 annotations.
- **`03_BAM_quality_control.sh`**: Evaluates alignment quality and RNA-seq-specific metrics using `SAMtools` and `Qualimap`.

  
**B. Long-read RNA-seq Processing**
