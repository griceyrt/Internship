#!/usr/bin/env Rscript
# =============================================================================
# Tests 9-12: adds a THIRD axis to the original 2x2x2 grid (level x padj/rawp
# x fc-cutoff-or-not) -- HOW the 2.5-fold cutoff is applied.
#
# Tests 5-8 filtered post-hoc: test against fold-change-0, then throw away
# anything with |log2FC| < log2(2.5) afterward. DESeq2 has a built-in,
# statistically different way to do this: results(..., lfcThreshold=log2(2.5),
# altHypothesis="greaterAbs") runs the Wald test AGAINST the threshold
# boundary directly, rather than against zero. This is generally more
# conservative and is very plausibly closer to what a paper means when it
# says "significant at >=2.5-fold" without spelling out post-hoc filtering.
#
# Gene-level dds is reused from qc_replicate_outlier_check.R (same 3+3
# samples, same design, already fit -- no need to refit). Isoform-level dds
# has to be refit (deseq2_isoform_level_test.R never saved it, only the
# results table) -- saved this time so it doesn't need refitting again.
# =============================================================================

library(tximport)
library(DESeq2)

BASE_DIR    <- "/home/grangel/fmr1_polya_splicing"
OARFISH_DIR <- file.path(BASE_DIR, "results", "gse188840_oarfish_2026-08-04")
QC_DIR      <- file.path(BASE_DIR, "results", "qc_replicate_check")
OUTPUT_DIR  <- file.path(BASE_DIR, "results", "de_lfcthreshold_test")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

FC_CUTOFF <- log2(2.5)

sample_names <- c("WT_Rep1", "WT_Rep2", "WT_Rep3", "KO_Rep1", "KO_Rep2", "KO_Rep3")
conditions   <- c("WT", "WT", "WT", "KO", "KO", "KO")
col_data <- data.frame(sample = sample_names,
                        condition = factor(conditions, levels = c("WT", "KO")),
                        row.names = sample_names)

# ---- gene level: reuse already-fit dds -------------------------------------
cat("Loading existing gene-level dds...\n")
dds_gene <- readRDS(file.path(QC_DIR, "dds_gene_level.rds"))

# ---- isoform level: refit (not saved anywhere previously) ------------------
cat("Rebuilding isoform-level dds (not previously saved)...\n")
quant_files <- file.path(OARFISH_DIR, paste0(sample_names, ".quant"))
names(quant_files) <- sample_names
txi_tx <- tximport(quant_files, type = "oarfish", txOut = TRUE, ignoreTxVersion = TRUE)
rownames(txi_tx$counts)    <- sub("\\.[0-9]+$", "", rownames(txi_tx$counts))
rownames(txi_tx$abundance) <- sub("\\.[0-9]+$", "", rownames(txi_tx$abundance))
rownames(txi_tx$length)    <- sub("\\.[0-9]+$", "", rownames(txi_tx$length))

dds_iso <- DESeqDataSetFromTximport(txi_tx, colData = col_data, design = ~condition)
dds_iso <- DESeq(dds_iso)
saveRDS(dds_iso, file.path(OUTPUT_DIR, "dds_isoform_level.rds"))

# ---- proper lfcThreshold Wald test, both levels -----------------------------
run_lfcthreshold_test <- function(dds, level_name) {
  res <- results(dds, contrast = c("condition", "KO", "WT"),
                  lfcThreshold = FC_CUTOFF, altHypothesis = "greaterAbs")
  df <- as.data.frame(res)
  write.table(df, file.path(OUTPUT_DIR, paste0(level_name, "_lfcThreshold_KO_vs_WT.tsv")),
              sep = "\t", quote = FALSE, row.names = TRUE)

  n_up_padj   <- sum(df$padj < 0.05 & df$log2FoldChange > 0, na.rm = TRUE)
  n_down_padj <- sum(df$padj < 0.05 & df$log2FoldChange < 0, na.rm = TRUE)
  n_up_p      <- sum(df$pvalue < 0.05 & df$log2FoldChange > 0, na.rm = TRUE)
  n_down_p    <- sum(df$pvalue < 0.05 & df$log2FoldChange < 0, na.rm = TRUE)

  data.frame(
    level     = level_name,
    threshold = c("padj", "rawp"),
    n_up      = c(n_up_padj, n_up_p),
    n_down    = c(n_down_padj, n_down_p),
    n_total   = c(n_up_padj + n_down_padj, n_up_p + n_down_p)
  )
}

cat("Running proper lfcThreshold test, gene level...\n")
summary_gene <- run_lfcthreshold_test(dds_gene, "gene")
cat("Running proper lfcThreshold test, isoform level...\n")
summary_iso  <- run_lfcthreshold_test(dds_iso, "isoform")

summary_all <- rbind(summary_gene, summary_iso)
summary_all$test <- 9:12
summary_all <- summary_all[, c("test", "level", "threshold", "n_up", "n_down", "n_total")]

write.table(summary_all, file.path(OUTPUT_DIR, "lfcthreshold_summary_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== Tests 9-12 (proper lfcThreshold Wald test, all at |log2FC| boundary = log2(2.5)) ===\n")
print(summary_all)
cat("\nFor direct comparison, tests 5-8 (post-hoc filtering, from the original grid):\n")
cat("  5  gene     padj  fc2.5    0 up    2 down    2 total\n")
cat("  6  gene     rawp  fc2.5    8 up   48 down   56 total\n")
cat("  7  isoform  padj  fc2.5   21 up   40 down   61 total\n")
cat("  8  isoform  rawp  fc2.5  455 up  715 down 1170 total\n")
cat("\nPaper's stated target: ~100 RNAs altered by at least 2.5-fold.\n")
