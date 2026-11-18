#!/usr/bin/env Rscript

# ==============================================================================
# Script: 11_generate_expression_matrices.R
# Purpose: Generate transcript- and gene-level expression matrices from
#          Salmon quantification results using the custom GTF annotation.
#
# Inputs:
#   1. Salmon quantification directory
#   2. Combined GENCODE v49 + novel transcript annotation (GTF)
#   3. Output directory
#
# Outputs:
#   Transcript estimated counts, transcript TPM, gene estimated counts,
#   gene TPM, and gene TPM with gene names.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Define Input Arguments and Output Directory
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    "Usage: Rscript 06_generate_expression_matrices.R ",
    "<salmon_dir> <combined_annotation.gtf> <output_dir>"
  )
}

salmon_dir <- args[1]
gtf_file <- args[2]
output_dir <- args[3]

project_name <- basename(normalizePath(salmon_dir, mustWork = TRUE))

if (!file.exists(gtf_file)) {
  stop("GTF annotation not found: ", gtf_file)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save processing messages in a log file
sink(file.path(output_dir, "output.log"), split = TRUE)

cat("Start time:", date(), "\n")
cat("Project:", project_name, "\n")
cat("Salmon directory:", salmon_dir, "\n")
cat("Annotation:", gtf_file, "\n")
cat("Output directory:", output_dir, "\n")

# ------------------------------------------------------------------------------
# 2. Load Required Libraries
# ------------------------------------------------------------------------------

library(dplyr)
library(readr)
library(rtracklayer)

# ------------------------------------------------------------------------------
# 3. Import Transcript Reference from Combined GTF
# ------------------------------------------------------------------------------

cat("\nImporting transcript annotation...\n")

gtf <- rtracklayer::import(gtf_file)
gtf_df <- as.data.frame(gtf)

# Select transcript records
gtf_df <- gtf_df %>%
  filter(type == "transcript")

# Add optional annotation fields when missing
if (!"gene_name" %in% names(gtf_df)) {
  gtf_df$gene_name <- gtf_df$gene_id
}

# Build transcript-to-gene reference table
transcripts_info <- gtf_df %>%
  transmute(
    transcript_id = as.character(transcript_id),
    gene_id = as.character(gene_id),
    gene_name = as.character(gene_name)
  )

# Replace missing gene names with gene IDs
transcripts_info <- transcripts_info %>%
  mutate(
    gene_name = if_else(
      is.na(gene_name) | gene_name == "",
      gene_id,
      gene_name
    )
  )

# Verify transcript identifiers and gene assignments
if (anyNA(transcripts_info$transcript_id) ||
    anyNA(transcripts_info$gene_id)) {
  stop("Missing transcript_id or gene_id in the GTF.")
}

if (anyDuplicated(transcripts_info$transcript_id)) {
  stop("Duplicate transcript IDs found in the combined annotation.")
}

cat("Annotated transcripts:", nrow(transcripts_info), "\n")

# ------------------------------------------------------------------------------
# 4. Identify Salmon Quantification Files
# ------------------------------------------------------------------------------

cat("\nIdentifying Salmon samples...\n")

sample_list <- sort(list.files(
  salmon_dir,
  pattern = "^quant\\.sf$",
  recursive = TRUE,
  full.names = TRUE
))

if (length(sample_list) == 0) {
  stop("No quant.sf files found in: ", salmon_dir)
}

sample_ids <- basename(dirname(sample_list))

if (anyDuplicated(sample_ids)) {
  stop("Duplicate sample IDs detected in Salmon directories.")
}

cat("Number of samples:", length(sample_list), "\n")

# ------------------------------------------------------------------------------
# 5. Generate Transcript-Level Expression Matrices
# ------------------------------------------------------------------------------

cat("\nBuilding transcript expression matrices...\n")

# Import a single Salmon quantification file
read_salmon <- function(file) {

  tab <- readr::read_tsv(
    file,
    show_col_types = FALSE,
    progress = FALSE
  )

  required <- c("Name", "TPM", "NumReads")

  if (!all(required %in% names(tab))) {
    stop("Missing required Salmon columns in: ", file)
  }

  if (anyDuplicated(tab$Name)) {
    stop("Duplicate transcript IDs in: ", file)
  }

  tab %>%
    transmute(
      transcript_id = as.character(Name),
      NumReads = as.numeric(NumReads),
      TPM = as.numeric(TPM)
    )
}

# Read all sample quantifications
quant_list <- lapply(sample_list, read_salmon)
names(quant_list) <- sample_ids

# Check that samples contain the same transcript reference
reference_ids <- quant_list[[1]]$transcript_id

for (i in seq_along(quant_list)) {

  ids <- quant_list[[i]]$transcript_id

  if (!setequal(reference_ids, ids)) {
    stop(
      "Transcript ID mismatch in sample: ",
      sample_ids[i]
    )
  }
}

# Check correspondence with the GTF annotation
missing_in_gtf <- setdiff(
  reference_ids,
  transcripts_info$transcript_id
)

missing_in_salmon <- setdiff(
  transcripts_info$transcript_id,
  reference_ids
)

cat("Salmon transcripts absent from GTF:",
    length(missing_in_gtf), "\n")

cat("GTF transcripts absent from Salmon:",
    length(missing_in_salmon), "\n")

if (length(missing_in_gtf) > 0 ||
    length(missing_in_salmon) > 0) {
  stop("Transcript IDs do not fully match the combined GTF.")
}

# Build matrices from each Salmon quantification
build_matrix <- function(column) {

  result <- data.frame(
    transcript_id = reference_ids,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(quant_list)) {

    tab <- quant_list[[i]]

    # Align transcript rows to the same reference order
    idx <- match(reference_ids, tab$transcript_id)

    result[[sample_ids[i]]] <- tab[[column]][idx]
  }

  result
}

numreads <- build_matrix("NumReads")
tpm <- build_matrix("TPM")

cat("Transcript count matrix:", dim(numreads), "\n")
cat("Transcript TPM matrix:", dim(tpm), "\n")

# ------------------------------------------------------------------------------
# 6. Generate Gene-Level Estimated Counts
# ------------------------------------------------------------------------------

cat("\nBuilding gene-level estimated counts...\n")

gene_mapping <- transcripts_info %>%
  select(transcript_id, gene_id)

genecounts <- numreads %>%
  left_join(gene_mapping, by = "transcript_id") %>%
  select(-transcript_id) %>%
  group_by(gene_id) %>%
  summarise(
    across(
      all_of(sample_ids),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

cat("Gene count matrix:", dim(genecounts), "\n")

# ------------------------------------------------------------------------------
# 7. Generate Gene-Level TPM
# ------------------------------------------------------------------------------

cat("\nBuilding gene-level TPM matrix...\n")

gtpm <- tpm %>%
  left_join(gene_mapping, by = "transcript_id") %>%
  select(-transcript_id) %>%
  group_by(gene_id) %>%
  summarise(
    across(
      all_of(sample_ids),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

cat("Gene TPM matrix:", dim(gtpm), "\n")

# ------------------------------------------------------------------------------
# 8. Generate Gene TPM Matrix with Gene Names
# ------------------------------------------------------------------------------

cat("\nAdding gene names to gene TPM matrix...\n")

gene_names <- transcripts_info %>%
  select(gene_id, gene_name) %>%
  distinct()

# Check consistent gene names per gene identifier
if (anyDuplicated(gene_names$gene_id)) {
  stop("Conflicting gene names for the same gene_id.")
}

# Gene names can be shared by different gene IDs.
# Warn because the gene-symbol output may have duplicate labels.
if (anyDuplicated(gene_names$gene_name)) {
  warning(
    "Some gene names are shared by multiple gene IDs. ",
    "The gene_id-based TPM table remains the unambiguous reference."
  )
}

gtpm_symbols <- gtpm %>%
  left_join(gene_names, by = "gene_id") %>%
  select(gene_name, all_of(sample_ids))

cat("Gene-name TPM matrix:", dim(gtpm_symbols), "\n")

# ------------------------------------------------------------------------------
# 9. Export Expression Tables
# ------------------------------------------------------------------------------

cat("\nWriting expression matrices...\n")

write_tsv(
  numreads,
  file.path(output_dir, paste0(project_name, "_tcounts.tsv"))
)

write_tsv(
  tpm,
  file.path(output_dir, paste0(project_name, "_tpm.tsv"))
)

write_tsv(
  genecounts,
  file.path(output_dir, paste0(project_name, "_gcounts.tsv"))
)

write_tsv(
  gtpm,
  file.path(output_dir, paste0(project_name, "_gtpm.tsv"))
)

write_tsv(
  gtpm_symbols,
  file.path(output_dir, paste0(project_name, "_gtpm_symbols.tsv"))
)

# ------------------------------------------------------------------------------
# 10. Pipeline Completion
# ------------------------------------------------------------------------------

cat("\n==========================================\n")
cat("Expression matrix generation completed.\n")
cat("Samples processed:", length(sample_ids), "\n")
cat("Transcripts:", nrow(numreads), "\n")
cat("Genes:", nrow(gtpm), "\n")
cat("Output directory:", output_dir, "\n")
cat("Completion time:", date(), "\n")
cat("==========================================\n")

sink()
