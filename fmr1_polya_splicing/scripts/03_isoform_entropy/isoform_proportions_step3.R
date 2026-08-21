#!/usr/bin/env Rscript
# =============================================================================
# Isoform proportion stacked-bar figure -- STEP 3 of 4
#
# Computes Shannon entropy of isoform usage per gene per condition, from
# step 2's isoform_proportions_long.tsv, then ranks genes by how much that
# entropy shifted between WT and KO. Independent of DESeq2/SUPPA gene lists
# by design -- entropy shift can happen with no change in total expression
# at all, which those tools wouldn't catch. Overlap #1's gene list is kept
# only as an annotation/cross-reference here, not used as a filter.
#
# Entropy formula: H = -sum(p * log2(p)) over a gene's isoform proportions
# in one condition. H = 0 means one isoform totally dominates (minimum
# diversity); higher H means usage is spread more evenly across isoforms.
# =============================================================================

BASE_DIR    <- "/home/grangel/fmr1_polya_splicing"
LONG_PATH   <- file.path(BASE_DIR, "results", "isoform_proportions", "isoform_proportions_long.tsv")
OVERLAP_PATH <- file.path(BASE_DIR, "results", "overlap1", "overlap1_gene_list.tsv")
OUTPUT_DIR  <- file.path(BASE_DIR, "results", "isoform_proportions")

# =============================================================================
# LOAD STEP 2 OUTPUT
# =============================================================================
cat("Loading isoform proportions (long format)...\n")
long_table <- read.table(LONG_PATH, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Rows loaded:", nrow(long_table), "\n")

# =============================================================================
# COMPUTE SHANNON ENTROPY PER (gene, condition)
# =============================================================================
cat("\nComputing Shannon entropy per gene per condition...\n")
shannon_entropy <- function(p) {
  p <- p[p > 0]   # safety net; proportions should already all be >0 here
  -sum(p * log2(p))
}

entropy_by_gene_cond <- aggregate(proportion ~ gene_id + condition,
                                   data = long_table, FUN = shannon_entropy)
colnames(entropy_by_gene_cond)[3] <- "entropy"

# Also compute each gene's TOTAL TPM per condition (sum across its present
# isoforms) -- needed to check whether extreme entropy shifts are backed
# by real expression or are noise from very low read counts (a gene with
# only a handful of reads can show an unstable entropy value by chance).
total_tpm_by_gene_cond <- aggregate(TPM ~ gene_id + condition,
                                     data = long_table, FUN = sum)
colnames(total_tpm_by_gene_cond)[3] <- "total_TPM"

# =============================================================================
# RESHAPE TO WIDE (one row per gene, entropy_WT and entropy_KO columns)
# =============================================================================
cat("Reshaping to wide format...\n")
entropy_WT <- entropy_by_gene_cond[entropy_by_gene_cond$condition == "WT", c("gene_id", "entropy")]
colnames(entropy_WT)[2] <- "entropy_WT"
entropy_KO <- entropy_by_gene_cond[entropy_by_gene_cond$condition == "KO", c("gene_id", "entropy")]
colnames(entropy_KO)[2] <- "entropy_KO"

entropy_wide <- merge(entropy_WT, entropy_KO, by = "gene_id", all = TRUE)

# Merge in total TPM per condition too
tpm_WT <- total_tpm_by_gene_cond[total_tpm_by_gene_cond$condition == "WT", c("gene_id", "total_TPM")]
colnames(tpm_WT)[2] <- "total_TPM_WT"
tpm_KO <- total_tpm_by_gene_cond[total_tpm_by_gene_cond$condition == "KO", c("gene_id", "total_TPM")]
colnames(tpm_KO)[2] <- "total_TPM_KO"
entropy_wide <- merge(entropy_wide, tpm_WT, by = "gene_id", all.x = TRUE)
entropy_wide <- merge(entropy_wide, tpm_KO, by = "gene_id", all.x = TRUE)

# Gene names, one lookup per gene (from the long table, which already has them)
gene_names <- unique(long_table[, c("gene_id", "gene_name")])
entropy_wide <- merge(entropy_wide, gene_names, by = "gene_id", all.x = TRUE)

# =============================================================================
# FLAG GENES PRESENT IN ONLY ONE CONDITION -- a different category (the
# whole gene switching on/off), not a comparable entropy shift, so these
# don't get an entropy_diff
# =============================================================================
only_one_condition <- is.na(entropy_wide$entropy_WT) | is.na(entropy_wide$entropy_KO)
cat("\nGenes present in only one condition (excluded from entropy-diff ranking):",
    sum(only_one_condition), "\n")
cat("  (of these, WT-only:", sum(is.na(entropy_wide$entropy_KO) & !is.na(entropy_wide$entropy_WT)),
    ", KO-only:", sum(is.na(entropy_wide$entropy_WT) & !is.na(entropy_wide$entropy_KO)), ")\n")

write.table(entropy_wide[only_one_condition, ],
            file      = file.path(OUTPUT_DIR, "entropy_single_condition_only.tsv"),
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

# =============================================================================
# COMPUTE ENTROPY SHIFT FOR GENES PRESENT IN BOTH CONDITIONS
# =============================================================================
entropy_both <- entropy_wide[!only_one_condition, ]
entropy_both$entropy_diff     <- entropy_both$entropy_KO - entropy_both$entropy_WT
entropy_both$abs_entropy_diff <- abs(entropy_both$entropy_diff)

# =============================================================================
# ANNOTATE (not filter) with Overlap #1 membership, for cross-reference
# =============================================================================
if (file.exists(OVERLAP_PATH)) {
  overlap1 <- read.table(OVERLAP_PATH, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  overlap1_genes <- unique(overlap1$gene_id)
  entropy_both$in_overlap1 <- entropy_both$gene_id %in% overlap1_genes
  cat("\nGenes flagged as also being in Overlap #1's list:", sum(entropy_both$in_overlap1), "\n")
} else {
  cat("\nOverlap #1 file not found, skipping cross-reference annotation.\n")
  entropy_both$in_overlap1 <- NA
}

# =============================================================================
# RANK BY ABSOLUTE ENTROPY SHIFT
# =============================================================================
entropy_ranked <- entropy_both[order(-entropy_both$abs_entropy_diff), ]

write.table(entropy_ranked,
            file      = file.path(OUTPUT_DIR, "entropy_shift_ranked.tsv"),
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

# =============================================================================
# LOOK AT THE DISTRIBUTION
# =============================================================================
cat("\n=== DISTRIBUTION OF ENTROPY SHIFT (genes present in both conditions) ===\n")
cat("Total genes ranked:", nrow(entropy_ranked), "\n\n")
cat("Summary of entropy_diff (signed, KO minus WT):\n")
print(summary(entropy_ranked$entropy_diff))
cat("\nSummary of abs_entropy_diff:\n")
print(summary(entropy_ranked$abs_entropy_diff))

cat("\nHow many genes at increasing abs_entropy_diff thresholds:\n")
for (thresh in c(0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0)) {
  cat("  > ", thresh, ":", sum(entropy_ranked$abs_entropy_diff > thresh), "genes\n")
}

cat("\nTop 15 genes by absolute entropy shift (now WITH total TPM, to check\n")
cat("whether these are backed by real expression or just low-N noise):\n")
print(head(entropy_ranked[, c("gene_id","gene_name","entropy_WT","entropy_KO",
                               "entropy_diff","total_TPM_WT","total_TPM_KO",
                               "in_overlap1")], 15))

cat("\n=== CHECKING FOR THE LOW-EXPRESSION-NOISE PATTERN ===\n")
cat("Among the top 100 genes by absolute entropy shift, how many have low\n")
cat("total expression (<10 TPM) in at least one condition -- these are the\n")
cat("ones most likely to be noise rather than a real biological shift:\n")
top100 <- head(entropy_ranked, 100)
low_expr_top100 <- sum(top100$total_TPM_WT < 10 | top100$total_TPM_KO < 10, na.rm = TRUE)
cat("  ", low_expr_top100, "out of 100\n")

cat("\n=== STEP 3 DONE ===\n")
cat("Saved:", file.path(OUTPUT_DIR, "entropy_shift_ranked.tsv"), "\n")
cat("Saved:", file.path(OUTPUT_DIR, "entropy_single_condition_only.tsv"), "\n")
