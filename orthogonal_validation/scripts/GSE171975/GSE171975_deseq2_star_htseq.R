#!/usr/bin/env Rscript
# GSE171975 (CT20 window) -- DESeq2 on STAR + htseq-count gene-level counts.
# Author: Gricey
# Second independent Asher dataset (WT n=4, PerDKO n=3, unbalanced by design).
# Same pipeline as GSE130613's deseq2_star_htseq.R.
# Run from: orthogonal_validation/

library(DESeq2)

BASE_DIR       <- "/Users/gricey/Desktop/Internship/orthogonal_validation"
RESULTS_FOLDER <- "GSE171975_CT20_2026-08-17"   # <-- UPDATE if different
COUNTS_DIR     <- file.path(BASE_DIR, "results", "GSE171975", RESULTS_FOLDER)
COUNTS_TSV     <- file.path(COUNTS_DIR, "combined_gene_counts.tsv")
OUTPUT_DIR     <- file.path(BASE_DIR, "results", "GSE171975", "normalisation_GSE171975_star_htseq")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

LFC_THRESH  <- 0.5
PADJ_THRESH <- 0.05

cat("Loading combined gene count matrix...\n")
raw <- read.table(COUNTS_TSV, sep = "\t", header = TRUE, row.names = 1,
                   check.names = FALSE)

gene_name <- raw$gene_name
count_mat <- as.matrix(raw[, !colnames(raw) %in% "gene_name"])
storage.mode(count_mat) <- "integer"
cat("Genes:", nrow(count_mat), " Samples:", ncol(count_mat), "\n")
cat("Sample columns:", paste(colnames(count_mat), collapse = ", "), "\n")

wt_samples <- c("WT_A", "WT_B", "WT_C", "WT_D")
ko_samples <- c("PerDKO_A", "PerDKO_B", "PerDKO_C")

condition <- ifelse(colnames(count_mat) %in% wt_samples, "WT",
             ifelse(colnames(count_mat) %in% ko_samples, "KO", NA))
stopifnot(!any(is.na(condition)))
col_data <- data.frame(condition = factor(condition, levels = c("WT", "KO")),
                        row.names = colnames(count_mat))
cat("Group sizes -- WT:", sum(condition == "WT"), " KO:", sum(condition == "KO"), "\n")

cat("\nRunning DESeq2...\n")
dds <- DESeqDataSetFromMatrix(countData = count_mat,
                               colData   = col_data,
                               design    = ~condition)

dds <- dds[rowSums(counts(dds)) >= 10, ]
cat("Genes after >=10 total count filter:", nrow(dds), "\n")

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "KO", "WT"))

out_df <- as.data.frame(res)
out_df$gene_id   <- rownames(out_df)
out_df$gene_name <- gene_name[match(out_df$gene_id, rownames(count_mat))]
out_df <- out_df[order(out_df$padj, na.last = TRUE), ]

write.table(out_df,
            file      = file.path(OUTPUT_DIR, "deseq2_GSE171975_KO_vs_WT.tsv"),
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = FALSE)

cat("\n--- DESeq2 results diagnostics ---\n")
cat("Total genes tested       :", nrow(out_df), "\n")
cat("log2FC range (all genes, incl. NaN-padj/unreliable):",
    round(min(out_df$log2FoldChange, na.rm = TRUE), 2), "to",
    round(max(out_df$log2FoldChange, na.rm = TRUE), 2), "\n")

reliable <- out_df[!is.na(out_df$padj), ]
cat("Genes with a real (non-NaN) padj (reliable)         :", nrow(reliable), "\n")
cat("log2FC range (reliable genes only)                  :",
    round(min(reliable$log2FoldChange), 2), "to",
    round(max(reliable$log2FoldChange), 2), "\n")

sig <- out_df[!is.na(out_df$padj) &
              out_df$padj < PADJ_THRESH &
              abs(out_df$log2FoldChange) > LFC_THRESH, ]
cat("Significant (padj<0.05, |log2FC|>0.5):", nrow(sig),
    "(", sum(sig$log2FoldChange > 0), "up,",
    sum(sig$log2FoldChange < 0), "down )\n")
cat("Sig log2FC range          :",
    round(min(sig$log2FoldChange), 2), "to",
    round(max(sig$log2FoldChange), 2), "\n")
cat("-----------------------------------\n")

cat("\n=== DONE ===\n")
cat("DESeq2 results:", file.path(OUTPUT_DIR, "deseq2_GSE171975_KO_vs_WT.tsv"), "\n")
