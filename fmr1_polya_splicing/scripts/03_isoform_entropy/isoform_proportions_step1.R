#!/usr/bin/env Rscript
# =============================================================================
# Isoform proportion stacked-bar figure -- STEP 1 of 4
#
# Goal: for a set of genes, show WT vs KO isoform usage as two stacked
# bars per gene, one segment per isoform, dropping isoforms with TPM == 0
# in that condition (a gene can show a different number of segments in
# each bar).
#
# This step: averages each condition's 3 replicates, applies the TPM>0
# filter per condition, and saves an intermediate table. Grouping by gene
# and computing proportions happens in step 2.
# =============================================================================

BASE_DIR   <- "/home/grangel/fmr1_polya_splicing"
TPM_PATH   <- file.path(BASE_DIR, "results", "normalisation", "TPM_all_samples.tpm")
OUTPUT_DIR <- file.path(BASE_DIR, "results", "isoform_proportions")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# LOAD TPM MATRIX
# =============================================================================
cat("Loading TPM matrix...\n")
tpm <- read.table(TPM_PATH, header = TRUE, sep = "\t", row.names = 1,
                  check.names = FALSE)
cat("Transcripts loaded:", nrow(tpm), "\n")
cat("Columns:", paste(colnames(tpm), collapse = ", "), "\n")

# =============================================================================
# AVERAGE EACH CONDITION'S 3 REPLICATES
# =============================================================================
cat("\nAveraging replicates per condition...\n")
wt_cols <- c("WT_Rep1", "WT_Rep2", "WT_Rep3")
ko_cols <- c("KO_Rep1", "KO_Rep2", "KO_Rep3")

avg_tpm <- data.frame(
  transcript_id = rownames(tpm),
  WT            = rowMeans(tpm[, wt_cols]),
  KO            = rowMeans(tpm[, ko_cols]),
  row.names     = NULL,
  stringsAsFactors = FALSE
)

# =============================================================================
# APPLY TPM > 0 FILTER, PER CONDITION SEPARATELY
# (a transcript can be "present" in WT but "absent" in KO, or vice versa --
#  that's the whole point, don't filter both conditions with the same rule
#  at once, keep the presence/absence flag separate per condition)
# =============================================================================
avg_tpm$present_in_WT <- avg_tpm$WT > 0
avg_tpm$present_in_KO <- avg_tpm$KO > 0

cat("\nTranscripts with TPM > 0 in WT:", sum(avg_tpm$present_in_WT), "\n")
cat("Transcripts with TPM > 0 in KO:", sum(avg_tpm$present_in_KO), "\n")
cat("Transcripts present in BOTH   :", sum(avg_tpm$present_in_WT & avg_tpm$present_in_KO), "\n")
cat("Transcripts present in NEITHER:", sum(!avg_tpm$present_in_WT & !avg_tpm$present_in_KO), "\n")

# =============================================================================
# SAVE INTERMEDIATE TABLE FOR INSPECTION
# =============================================================================
write.table(avg_tpm,
            file      = file.path(OUTPUT_DIR, "avg_tpm_per_condition.tsv"),
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

cat("\n=== STEP 1 DONE ===\n")
cat("Saved:", file.path(OUTPUT_DIR, "avg_tpm_per_condition.tsv"), "\n")
