# RNA-seq Pipeline for Novel long noncoding RNA (lncRNAs) Identification and Quantification in sepsis

This repository contains the Linux-based bioinformatics pipeline developed to identify and quantify novel long noncoding sepsis-associated transcripts (SEATs) by integrating long-read and short-read RNA sequencing data. It combines long-read RNA-seq (n = 10 samples) and short-read RNA-seq (n = 10 samples) for hybrid transcriptome assembly. Both sequencing datasets are first processed and aligned independently before being jointly used for transcript reconstruction. The resulting transcript assemblies are merged to generate a unified transcriptome annotation, including previously unannotated transcripts. This annotation then undergoes a series of filtering steps to identify and retain high-confidence novel long noncoding SEATs candidates. The resulting custom annotation is subsequently used to quantify transcript expression across n = 812 additional short-read RNA-seq samples, enabling the characterisation of newly identified transcripts in a larger cohort.

## Pipeline Architecture
The pipeline is organized into three main phases:

**Phase 1: RNA-seq Pre-processing and Alignment**

**A. Short-read RNA-seq**

- **`01_read_QC.sh`**: Downloads raw paired-end FASTQ files, verifies file integrity using MD5 checksums, and merges sequencing lanes for each sample. Checks R1/R2 read concordance and performs adapter trimming and quality filtering using `fastp`. Sequencing quality is assessed using `FastQC` and `FastQ Screen`, with QC metrics consolidated into a single MultiQC report.
- **`02_alignment.sh`**: Aligns fastp-processed paired-end reads to the human reference genome (GRCh38) using `STAR`, with GENCODE v49 annotations.
- **`03_BAM_QC.sh`**: Assesses STAR-aligned BAM files using SAMtools and Qualimap. Generates alignment statistics, evaluates BAM quality metrics, and performs RNA-seq-specific quality assessment using GENCODE v49 annotations and a reverse-stranded paired-end library configuration.

  
**B. Long-read RNA-seq**

- **`01_read_QC.sh`**: Downloads raw FASTQ files, verifies file integrity using MD5 checksums, and assesses read length distributions and sequencing quality using `NanoPlot`.
- **`02_filtering.sh`**: Filters Oxford Nanopore long-read RNA-seq data using `Chopper`, and generates compressed, quality-filtered FASTQ files for downstream splice-aware alignment.
- **`03_alignment.sh`**: Aligns quality-filtered long RNA-seq reads to the human reference genome (GRCh38) using `minimap2`, generating genome-aligned reads for downstream transcriptome assembly.
- **`04_BAM_QC.sh`**: Evaluates long-read alignment quality using `Cramino`, generating alignment statistics and quality metrics.


**Phase 2: Hybrid Transcriptome Assembly**

- **`05_hybrid_assembly.sh`**: Performs reference-guided hybrid transcriptome assembly by integrating matched short-read and long-read RNA-seq alignments using `StringTie`.

**Phase 3: Novel Transcript Annotation and Coding Potential Assessment**

- **`06_gffcompare_filtering.sh`**: Compares the merged hybrid transcriptome assembly against GENCODE v49 using `gffcompare`.
- **`07_coding_potential.sh`**: Extracts transcript sequences using gffread and evaluates protein-coding potential using three complementary prediction methods: `CNCI`, `CPC2`, and `CPAT`, implemented through `ezLncPred`.
- **`08_noncoding_transcript_filtering.R`**: Integrates coding-potential predictions.

**Phase 4: Transcript Quantification**

- **`09_custom_salmon_index.sh`**: Combines GENCODE v49 reference annotations with filtered novel transcript candidates to generate a custom transcriptome annotation. Extracts transcript sequences using `gffread` and constructs a decoy-aware Salmon index by incorporating the GRCh38 reference genome, enabling subsequent transcript-level quantification.
- **`10_salmon_quantification.sh`**: Quantifies transcript abundance in paired-end short-read RNA-seq samples using Salmon and the custom decoy-aware transcriptome index combining GENCODE v49 annotations with newly identified transcript candidates.
- **`11_generate_expression_matrices.R`**: Integrates sample-level Salmon quantification results and the custom transcript annotation to generate transcript- and gene-level expression matrices.

### Dependencies

**fastp**

**FastQC, FastQ Screen, MultiQC**

**STAR**

**NanoPlot, Chopper**

**minimap2**

**SAMtools**

**Qualimap, Cramino**

**StringTie**

**gffcompare, gffread**

**ezLncPred (CNCI, CPC2, CPAT)**

**Salmon**

**GNU Utilities (wget, md5sum, gzip, awk, sed, sort)**

**`R packages: dplyr, readr, stringr, tibble, magrittr, rtracklayer`**
