#!/usr/bin/env Rscript
# =============================================================================
# Isoform proportion stacked-bar figure -- STEP 4 of 4
#
# Adds a minimum total-expression filter before ranking by entropy shift
# (protects against low-read-count genes producing an unstable entropy
# value by chance), tested at 10/20/50 TPM thresholds to compare the
# tradeoff between stricter filtering and gene count retained.
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
# CHECK MULTIPLE THRESHOLDS IN ONE RUN
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
