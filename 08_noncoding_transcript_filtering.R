
#!/usr/bin/env Rscript

# ==============================================================================
# Script: 08_noncoding_transcript_filtering.R
# Purpose: Integrate CNCI, CPC2, and CPAT predictions and retain candidate
#          noncoding transcripts satisfying all three filtering criteria.
# ==============================================================================

library(dplyr)
library(tibble)
library(stringr)

# ------------------------------------------------------------------------------
# 1. Define Directories and Filtering Thresholds
# ------------------------------------------------------------------------------

project_dir <- "/data/records_wts/03_hybrid_assembly"

annotation_dir <- file.path(project_dir, "annotation")
coding_dir <- file.path(project_dir, "coding_potential")

gtf_file <- file.path(annotation_dir, "merged_pass2.gtf")

cnci_file <- file.path(coding_dir, "cnci/results/CNCI.index")
cpc2_file <- file.path(coding_dir, "cpc2/results")
cpat_file <- file.path(coding_dir, "cpat/results")

out_dir <- file.path(project_dir, "final_annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

CPAT_THRESHOLD <- 0.364
CPC2_THRESHOLD <- 0.5
CNCI_THRESHOLD <- 0

# ------------------------------------------------------------------------------
# 2. Import Coding Potential Predictions
# ------------------------------------------------------------------------------

cnci <- read.table(
  cnci_file,
  sep = "\t",
  header = TRUE
) %>%
  dplyr::select(Transcript.ID, score)

cpc2 <- read.table(
  cpc2_file,
  sep = "\t",
  header = TRUE,
  comment.char = ""
) %>%
  dplyr::rename(Transcript.ID = X.ID) %>%
  dplyr::select(Transcript.ID, coding_probability)

# Preserve your original CPAT ID-processing procedure.
# This requires verification against the actual CPAT output format.

cpat <- read.table(
  cpat_file,
  sep = "\t",
  header = TRUE
) %>%
  rownames_to_column(var = "Transcript.ID") %>%
  mutate(
    Transcript.ID = sub("^b", "", Transcript.ID)
  ) %>%
  dplyr::select(Transcript.ID, coding_prob)

# ------------------------------------------------------------------------------
# 3. Integrate Prediction Results
# ------------------------------------------------------------------------------

# Avoid unintended many-to-many joins.
for (df in list(cnci, cpc2, cpat)) {
  if (anyDuplicated(df$Transcript.ID) > 0) {
    stop("Duplicate transcript identifiers detected.")
  }
}

coding_potential <- cnci %>%
  inner_join(cpat, by = "Transcript.ID") %>%
  inner_join(cpc2, by = "Transcript.ID")

write.table(
  coding_potential,
  file.path(out_dir, "coding_potential_all.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 4. Apply Noncoding Potential Thresholds
# ------------------------------------------------------------------------------

coding_potential_f <- coding_potential %>%
  filter(
    coding_prob < CPAT_THRESHOLD,
    coding_probability < CPC2_THRESHOLD,
    score < CNCI_THRESHOLD
  ) %>%
  rename(transcript_id = Transcript.ID)

write.table(
  coding_potential_f,
  file.path(out_dir, "coding_potential_filtered.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 5. Filter the Transcript Annotation
# ------------------------------------------------------------------------------

# Read GTF as text to preserve the original annotation records.
gtf_lines <- readLines(gtf_file)

# Extract transcript IDs from GTF column 9.
# The expected format is: transcript_id "ID";

get_transcript_id <- function(line) {
  fields <- strsplit(line, "\t", fixed = TRUE)[[1]]

  if (length(fields) < 9) return(NA_character_)

  str_match(
    fields[9],
    '(?:^|;)\\s*transcript_id "([^"]+)"'
  )[, 2]
}

gtf_ids <- vapply(
  gtf_lines,
  get_transcript_id,
  character(1)
)

# Retain transcript and exon records for selected IDs.
keep <- !is.na(gtf_ids) &
  gtf_ids %in% coding_potential_f$transcript_id

merged_pass4 <- gtf_lines[keep]

writeLines(
  merged_pass4,
  file.path(out_dir, "merged_pass4.gtf")
)

# ------------------------------------------------------------------------------
# 6. Summary
# ------------------------------------------------------------------------------

cat("==========================================\n")
cat("Coding potential filtering completed.\n")
cat("Transcripts with predictions from all three tools:",
    nrow(coding_potential), "\n")
cat("Transcripts passing all thresholds:",
    nrow(coding_potential_f), "\n")
cat("Final GTF:", file.path(out_dir, "merged_pass4.gtf"), "\n")
cat("==========================================\n")
