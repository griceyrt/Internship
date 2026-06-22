#!/usr/bin/env Rscript
# =============================================================================
# Orthogonal Validation — Steps 5 & 6
# Author: Gricey
# Description:
#   Step 5: tximport → filter → normalize → (optional limma) → SUPPA input
#   Step 6: DESeq2 at TRANSCRIPT level + stageR two-stage adjustment
#           Stage 1 (screen):  gene-level q-values (perGeneQValueHack)
#           Stage 2 (confirm): transcript-level p-values within significant genes
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
# BiocManager::install(c("tximport", "DESeq2", "limma", "stageR"))

# =============================================================================
# LIBRARIES
# =============================================================================
library(tximport)
library(DESeq2)
library(limma)
library(stageR)

# =============================================================================
# PATHS — update SALMON_DATE to match your results folder
# =============================================================================
BASE_DIR    <- "/Users/gricey/Desktop/Internship/orthogonal_validation"
SALMON_DATE <- "2026-06-20"   # <-- UPDATE if different
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
# STEP 6 — Transcript-level DE with stageR two-stage adjustment
#
# Why transcript level?
#   A gene can have multiple transcripts. We want to know not just IF a gene
#   changes, but WHICH specific transcript changes. stageR does this in two stages:
#     Stage 1 (screen):  gene-level — does anything change in this gene at all?
#     Stage 2 (confirm): transcript-level — which specific transcript changed?
#
# Functions from Bharath's code (perGeneQValueHack + perGeneQValueExact):
#   These compute gene-level q-values by summarising transcript p-values per gene
#   (taking the minimum p-value per gene), then applying BH correction.
# =============================================================================

# --- Helper functions (from Bharath) -----------------------------------------

perGeneQValueExact <- function(pGene, theta, geneSplit) {
  stopifnot(length(pGene) == length(geneSplit))
  numExons <- listLen(geneSplit)
  tab      <- tabulate(numExons)
  notZero  <- (tab > 0)
  numerator <- mapply(function(m, n) m * (1 - (1 - theta)^n),
                      m = tab[notZero], n = which(notZero))
  numerator <- rowSums(numerator)
  bins  <- cut(pGene, breaks = c(-Inf, as.vector(theta)), right = TRUE,
               include.lowest = TRUE)
  counts <- tabulate(bins, nbins = nlevels(bins))
  denom  <- cumsum(counts)
  return(numerator / denom)
}

perGeneQValueHack <- function(object, p = "pvalue",
                               method = perGeneQValueExact) {
  wTest    <- which(!is.na(object$padj))
  pvals    <- object[[p]][wTest]
  geneID   <- factor(object[["groupID"]][wTest])
  geneSplit <- split(seq_along(geneID), geneID)
  pGene    <- sapply(geneSplit, function(i) min(pvals[i]))
  stopifnot(all(is.finite(pGene)))
  theta <- unique(sort(pGene))
  q     <- method(pGene, theta, geneSplit)
  res        <- rep(NA_real_, length(pGene))
  res        <- q[match(pGene, theta)]
  res        <- pmin(1, res)
  names(res) <- names(geneSplit)
  stopifnot(!any(is.na(res)))
  return(res)
}

# --- DESeq2 at transcript level -----------------------------------------------

cat("\nRunning DESeq2 (transcript level)...\n")

# Use the same txi_tx from Step 5a (transcript-level, txOut=TRUE)
# Build DESeq2 object directly from transcript-level tximport
dds_tx <- DESeqDataSetFromTximport(txi_tx,
                                   colData = col_data,
                                   design  = ~ condition)

# Add gene_id to rowData so stageR can group transcripts by gene
tx_to_gene              <- setNames(tx2gene$gene_id, tx2gene$tx_id)
rowData(dds_tx)$groupID <- tx_to_gene[rownames(dds_tx)]

# Remove transcripts not found in GTF (no gene_id mapping — NA groupID)
# These exist in transcriptome_ext.fa but not in transcriptome_productivity.gtf
no_mapping <- is.na(rowData(dds_tx)$groupID)
cat("Transcripts with no GTF gene mapping (removed):", sum(no_mapping), "\n")
dds_tx <- dds_tx[!no_mapping, ]

# Filter: keep transcripts with total counts > 20 across all samples
keep_tx <- rowSums(counts(dds_tx)) > 20
dds_tx  <- dds_tx[keep_tx, ]
cat("Transcripts after DESeq2 filter (rowSums > 20):", nrow(dds_tx), "\n")

# Run DESeq2
dds_tx <- DESeq(dds_tx)

# Extract results: KO vs WT
res_tx <- results(dds_tx,
                  contrast            = c("condition", "KO", "WT"),
                  independentFiltering = TRUE,
                  alpha               = 0.05)

# Add groupID (gene_id) to results for stageR
res_tx$groupID <- rowData(dds_tx)$groupID

# --- Diagnostics before stageR -----------------------------------------------
cat("\n--- DESeq2 raw results diagnostics ---\n")
cat("Total transcripts tested          :", nrow(res_tx), "\n")
cat("Transcripts with non-NA pvalue    :", sum(!is.na(res_tx$pvalue)), "\n")
cat("Transcripts with non-NA padj      :", sum(!is.na(res_tx$padj)), "\n")
cat("Transcripts with pvalue < 0.05    :", sum(res_tx$pvalue < 0.05, na.rm=TRUE), "\n")
cat("Transcripts with padj < 0.05      :", sum(res_tx$padj < 0.05, na.rm=TRUE), "\n")
cat("Transcripts with non-NA groupID   :", sum(!is.na(res_tx$groupID)), "\n")
cat("--------------------------------------\n\n")

# --- Stage 1: gene-level screening -------------------------------------------

cat("Running stageR stage 1 (gene-level screening)...\n")
pScreen <- perGeneQValueHack(res_tx)   # named vector: gene_id -> q-value
cat("pScreen summary:\n")
print(summary(pScreen))
cat("Genes with pScreen < 0.05:", sum(pScreen < 0.05, na.rm=TRUE), "\n")

# --- Stage 2: transcript-level confirmation -----------------------------------

cat("Running stageR stage 2 (transcript-level confirmation)...\n")
pConfirmation <- matrix(res_tx$pvalue, ncol = 1)
dimnames(pConfirmation) <- list(rownames(res_tx), "transcript")

# tx2gene subset: only transcripts in the analysis
tx2gene_sub <- tx2gene[tx2gene$tx_id %in% rownames(res_tx), ]

stageRObj <- stageRTx(pScreen      = pScreen,
                       pConfirmation = pConfirmation,
                       pScreenAdjusted = TRUE,
                       tx2gene      = tx2gene_sub)

stageRObj <- stageWiseAdjustment(stageRObj,
                                  method  = "dte",
                                  alpha   = 0.05,
                                  allowNA = TRUE)

suppressWarnings({
  padj_stageR <- getAdjustedPValues(stageRObj,
                                     order                 = FALSE,
                                     onlySignificantGenes  = FALSE)
})
cat("padj_stageR class:", class(padj_stageR), "\n")
cat("padj_stageR dim:", dim(padj_stageR), "\n")
cat("padj_stageR colnames:", colnames(padj_stageR), "\n")
cat("padj_stageR head:\n")
print(head(padj_stageR))
cat("Values < 0.05 in gene col:", sum(padj_stageR[,"gene"] < 0.05, na.rm=TRUE), "\n")
cat("Values < 0.05 in tx col:  ", sum(padj_stageR[,"transcript"] < 0.05, na.rm=TRUE), "\n")

# --- Build output table -------------------------------------------------------

# gene_name map (gene_id -> gene_name)
gene_name_map_by_gene <- setNames(
  sapply(gtf_tx$attributes, extract_attr, key = "gene_name"),
  sapply(gtf_tx$attributes, extract_attr, key = "gene_id")
)

# Base table from DESeq2 results
res_df               <- as.data.frame(res_tx)
res_df$transcript_id <- rownames(res_df)
res_df$gene_id       <- res_df$groupID
res_df$gene_name     <- gene_name_map_by_gene[res_df$gene_id]
res_df$gene_name[is.na(res_df$gene_name)] <- "unannotated"

# Merge stageR output — use txID column (NOT rownames)
padj_df <- padj_stageR
colnames(padj_df)[colnames(padj_df) == "txID"]        <- "transcript_id"
colnames(padj_df)[colnames(padj_df) == "gene"]        <- "padj_gene_stageR"
colnames(padj_df)[colnames(padj_df) == "transcript"]  <- "padj_tx_stageR"

res_df <- merge(res_df,
                padj_df[, c("transcript_id", "padj_gene_stageR", "padj_tx_stageR")],
                by = "transcript_id", all.x = TRUE)

# Select and order final columns
out_df <- res_df[, c("transcript_id", "gene_id", "gene_name",
                      "log2FoldChange", "pvalue",
                      "padj_gene_stageR", "padj_tx_stageR")]
out_df <- out_df[order(out_df$padj_gene_stageR, out_df$padj_tx_stageR, na.last = TRUE), ]

write.table(out_df,
            file      = file.path(OUTPUT_DIR, "deseq2_KO_vs_WT.tsv"),
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

# Summary
sig_genes <- length(unique(out_df$gene_id[!is.na(out_df$padj_gene_stageR) & out_df$padj_gene_stageR < 0.05]))
sig_tx    <- sum(!is.na(out_df$padj_tx_stageR) & out_df$padj_tx_stageR < 0.05)
cat("Significant genes (stageR stage 1, padj<0.05)      :", sig_genes, "\n")
cat("Significant transcripts (stageR stage 2, padj<0.05):", sig_tx, "\n")
cat("Saved: deseq2_KO_vs_WT.tsv\n")

cat("\n=== DONE ===\n")
cat("SUPPA input table : results/normalisation/combined_norm.tab\n")
cat("DESeq2 results    : results/normalisation/deseq2_KO_vs_WT.tsv\n")
