#!/usr/bin/env Rscript
# =============================================================================
# Orthogonal Validation — Steps 5 & 6
# Author: Gricey
# Description:
#   Step 5: tximport → filter → normalize → (optional limma) → SUPPA input
#   Step 6: DESeq2 differential expression (WT vs PerDKO)
#
# Dataset: GSE130613 — WT and Per1/2 KO mouse liver, constant darkness, CT16-20
#   WT  (n=4): SRR9002567, SRR9002568, SRR9002569, SRR9002570
#   KO  (n=4): SRR9002583, SRR9002584, SRR9002585, SRR9002586
#
# Run from: orthogonal_validation/
# =============================================================================

# =============================================================================
# INSTALL PACKAGES (run once if needed)
# =============================================================================
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(c("tximport", "DESeq2", "limma"))

# =============================================================================
# LIBRARIES
# =============================================================================
library(tximport)
library(DESeq2)
library(limma)

# =============================================================================
# PATHS — update SALMON_DATE to match your results folder
# =============================================================================
BASE_DIR    <- "/Users/gricey/Desktop/Internship/orthogonal_validation"
SALMON_DATE <- "2026-06-16"   # <-- UPDATE if different
SALMON_DIR  <- file.path(BASE_DIR, "results", paste0("SRP194523_salmon_", SALMON_DATE))
GTF_PATH    <- "/Users/gricey/Desktop/Internship/boundary_analysis/data/transcriptome_productivity.gtf"
OUTPUT_DIR  <- file.path(BASE_DIR, "results", "normalisation")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# SAMPLE INFO
# =============================================================================
sample_ids  <- c("SRR9002567", "SRR9002568", "SRR9002569", "SRR9002570",
                 "SRR9002583", "SRR9002584", "SRR9002585", "SRR9002586")
sample_names <- c("WT_a", "WT_b", "WT_c", "WT_d",
                  "KO_a", "KO_b", "KO_c", "KO_d")
conditions  <- c("WT", "WT", "WT", "WT", "KO", "KO", "KO", "KO")

col_data <- data.frame(
  sample    = sample_names,
  condition = factor(conditions, levels = c("WT", "KO")),
  row.names = sample_names
)

# quant.sf paths
quant_files <- file.path(SALMON_DIR, sample_ids, "quant.sf")
names(quant_files) <- sample_names

# Check all files exist
missing <- quant_files[!file.exists(quant_files)]
if (length(missing) > 0) stop("Missing quant.sf files:\n", paste(missing, collapse="\n"))
cat("All 8 quant.sf files found.\n")

# =============================================================================
# BUILD tx2gene FROM GTF
# (maps transcript_id -> gene_id, needed for DESeq2 gene-level aggregation)
# =============================================================================
cat("Building tx2gene from GTF...\n")
gtf <- read.table(GTF_PATH, sep="\t", quote="", comment.char="#",
                  col.names=c("chr","source","feature","start","end",
                              "score","strand","frame","attributes"))
gtf_tx <- gtf[gtf$feature == "transcript", ]

extract_attr <- function(attr_string, key) {
  pattern <- paste0(key, ' "([^"]+)"')
  m <- regmatches(attr_string, regexpr(pattern, attr_string))
  if (length(m) == 0) return(NA)
  sub(paste0(key, ' "([^"]+)"'), "\\1", m)
}

tx2gene <- data.frame(
  tx_id   = sapply(gtf_tx$attributes, extract_attr, key = "transcript_id"),
  gene_id = sapply(gtf_tx$attributes, extract_attr, key = "gene_id"),
  stringsAsFactors = FALSE,
  row.names = NULL
)
tx2gene <- tx2gene[!is.na(tx2gene$tx_id) & !is.na(tx2gene$gene_id), ]
cat("tx2gene built:", nrow(tx2gene), "transcripts.\n")

# =============================================================================
# STEP 5a — tximport (transcript-level, for SUPPA)
# =============================================================================
cat("Running tximport (transcript level)...\n")
txi_tx <- tximport(quant_files,
                   type   = "salmon",
                   txOut  = TRUE)   # transcript-level, no gene aggregation

raw_counts_tx <- txi_tx$counts
cat("Transcripts imported:", nrow(raw_counts_tx), "\n")

# =============================================================================
# STEP 5b — Filter: keep transcripts with >5 counts in >=2 samples
# =============================================================================
cat("Filtering low-expressed transcripts...\n")
keep <- rowSums(raw_counts_tx > 5) >= 2
counts_filtered <- raw_counts_tx[keep, ]
cat("Transcripts after filtering:", nrow(counts_filtered),
    "(removed:", sum(!keep), ")\n")

# =============================================================================
# STEP 5c — Normalize: counts / total reads per sample
# =============================================================================
cat("Normalizing to total reads per sample...\n")
total_reads     <- colSums(counts_filtered)
counts_norm     <- sweep(counts_filtered, 2, total_reads, FUN = "/")

# =============================================================================
# STEP 5d — limma removeBatchEffect (OPTIONAL)
# =============================================================================
# NOTE: The SI Appendix of the hypoxia paper (GSE130613) does not provide
# explicit batch information for the DD samples we are using.
# Options:
#   A) Skip batch correction (set APPLY_BATCH_CORRECTION <- FALSE)
#   B) If Bharath/Kiran provide batch info, define batch_vector and set to TRUE
#
# Example batch vector if a/b = batch1, c/d = batch2:
#   batch_vector <- c(1, 1, 2, 2, 1, 1, 2, 2)

APPLY_BATCH_CORRECTION <- FALSE  # <-- change to TRUE if batch info available
# batch_vector <- c(...)          # <-- define if APPLY_BATCH_CORRECTION = TRUE

if (APPLY_BATCH_CORRECTION) {
  cat("Applying limma removeBatchEffect...\n")
  counts_suppa <- removeBatchEffect(counts_norm, batch = batch_vector)
} else {
  cat("Skipping batch correction (no batch info available).\n")
  counts_suppa <- counts_norm
}

# =============================================================================
# STEP 8 — Export expression table for SUPPA
# Format: rows = transcript IDs, columns = sample names, tab-separated
# =============================================================================
cat("Exporting SUPPA expression table...\n")
suppa_table <- as.data.frame(counts_suppa)
colnames(suppa_table) <- sample_names

write.table(suppa_table,
            file      = file.path(OUTPUT_DIR, "combined_norm.tab"),
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = TRUE)
cat("Saved: combined_norm.tab\n")

# =============================================================================
# STEP 6 — DESeq2 differential expression (gene-level, parallel track)
# Uses raw estimated counts from tximport aggregated to gene level
# =============================================================================
cat("\nRunning DESeq2 (gene-level)...\n")

# tximport at gene level for DESeq2
txi_gene <- tximport(quant_files,
                     type    = "salmon",
                     tx2gene = tx2gene,
                     txOut   = FALSE)

# Build DESeq2 object
dds <- DESeqDataSetFromTximport(txi_gene,
                                colData = col_data,
                                design  = ~ condition)

# Run DESeq2
dds <- DESeq(dds)

# Results: KO vs WT
res <- results(dds, contrast = c("condition", "KO", "WT"))
res_df <- as.data.frame(res)
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

# Save results
write.table(res_df,
            file      = file.path(OUTPUT_DIR, "deseq2_KO_vs_WT.tsv"),
            sep       = "\t",
            quote     = FALSE,
            col.names = NA)

# Summary
sig <- sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1.5, na.rm = TRUE)
cat("DESeq2 significant genes (|log2FC|>1.5, padj<0.05):", sig, "\n")
cat("Saved: deseq2_KO_vs_WT.tsv\n")

cat("\n=== DONE ===\n")
cat("SUPPA input table : results/normalisation/combined_norm.tab\n")
cat("DESeq2 results    : results/normalisation/deseq2_KO_vs_WT.tsv\n")
