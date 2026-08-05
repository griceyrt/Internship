#!/usr/bin/env Rscript
# =============================================================================
# Isoform proportion stacked-bar figure -- STEP 4 of N
# Author: Gricey
#
# Step 3 confirmed the low-expression-noise concern: 92/100 top genes by
# entropy shift had total TPM < 10 in at least one condition, several as
# low as 0.33-1.0 TPM (almost certainly 1-2 stray reads, not real isoform
# usage). This step adds a minimum total-expression filter BEFORE ranking
# by entropy shift, tested at a few threshold levels so we can see the
# tradeoff (stricter threshold = fewer but more trustworthy genes) rather
# than picking one blindly.
# =============================================================================

BASE_DIR   <- "/home/grangel/fmr1_polya_splicing"
RANKED_PATH <- file.path(BASE_DIR, "results", "isoform_proportions", "entropy_shift_ranked.tsv")
OUTPUT_DIR  <- file.path(BASE_DIR, "results", "isoform_proportions")

cat("Loading entropy-ranked table from Step 3...\n")
entropy_ranked <- read.table(RANKED_PATH, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Total genes:", nrow(entropy_ranked), "\n\n")

# A gene needs decent expression in BOTH conditions to trust its entropy
# calculation in either one -- require the MINIMUM of the two totals to
# clear the threshold, not just one of them
entropy_ranked$min_total_TPM <- pmin(entropy_ranked$total_TPM_WT, entropy_ranked$total_TPM_KO)

cat("=== HOW MANY GENES SURVIVE AT EACH MINIMUM-EXPRESSION THRESHOLD ===\n")
for (thresh in c(0, 5, 10, 20, 50, 100)) {
  n <- sum(entropy_ranked$min_total_TPM > thresh)
  cat("  min(total_TPM_WT, total_TPM_KO) > ", thresh, ": ", n, "genes\n", sep = "")
}

# =============================================================================
# CHECK MULTIPLE THRESHOLDS IN ONE RUN, per Gricey's request 2026-08-05, to
# compare 10 vs 20 (already looked at) vs 50 TPM sensitivity before picking
# =============================================================================
for (THRESHOLD in c(10, 20, 50)) {
  cat("\n\n########################################\n")
  cat("=== TOP 15 GENES BY ENTROPY SHIFT, AFTER REQUIRING min_total_TPM >", THRESHOLD, "===\n")
  cat("########################################\n")
  filtered <- entropy_ranked[entropy_ranked$min_total_TPM > THRESHOLD, ]
  filtered <- filtered[order(-filtered$abs_entropy_diff), ]
  cat("Genes remaining after filter:", nrow(filtered), "of", nrow(entropy_ranked), "\n\n")
  print(head(filtered[, c("gene_id","gene_name","entropy_WT","entropy_KO","entropy_diff",
                           "total_TPM_WT","total_TPM_KO","in_overlap1")], 15))

  write.table(filtered,
              file      = file.path(OUTPUT_DIR, paste0("entropy_shift_ranked_min", THRESHOLD, "TPM.tsv")),
              sep       = "\t",
              quote     = FALSE,
              row.names = FALSE)
}

cat("\n\n=== STEP 4 DONE ===\n")
cat("Saved one ranked file per threshold (10, 20, 50 TPM) in", OUTPUT_DIR, "\n")
cat("Compare the three top-15 lists above: how much does the actual gene\n")
cat("list change between 10 and 50? If it's fairly stable, the choice of\n")
cat("exact threshold doesn't matter much. If it shifts a lot, that's worth\n")
cat("discussing with Kiran before finalizing.\n")
