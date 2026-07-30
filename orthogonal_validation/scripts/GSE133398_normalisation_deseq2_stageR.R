#!/usr/bin/env Rscript
# =============================================================================
# GSE133398 — Option B: transcript-level DESeq2 + stageR
# Author: Gricey
# Adapted from: normalisation_deseq2.R (GSE130613) — same stageR pattern,
#   including Bharath's perGeneQValueHack/perGeneQValueExact functions verbatim.
#
# Why this exists: Bharath's answer on gene-level vs transcript-level was
# conditional ("depends on how we're using the result") — this is Option B
# (transcript-level, stageR two-stage), to sit alongside the existing
# Option A (gene-level, GSE133398_normalisation_deseq2.R) so he/Kiran can
# pick based on whether they want widespread-effects framing (A) or
# splicing-emphasis framing (B).
#
# Dataset: GSE133398 — SRSF3 KO vs control hepatocytes (mm39, 2026-07-08 run)
#   WT / GFP control (n=4): SRR9606472, SRR9606473, SRR9606474, SRR9606475
#   KO / CRE SRSF3 KO (n=4): SRR9606468, SRR9606469, SRR9606470, SRR9606471
#
# Run from: orthogonal_validation/
# =============================================================================

# =============================================================================
# LIBRARIES — stageR is new for this script, check it's installed first:
#   BiocManager::install("stageR")
# =============================================================================
library(tximport)
library(DESeq2)
library(limma)
library(stageR)

# =============================================================================
# PATHS
# =============================================================================
BASE_DIR    <- "/Users/gricey/Desktop/Internship/orthogonal_validation"
SALMON_DATE <- "2026-07-08"
SALMON_DIR  <- file.path(BASE_DIR, "results", paste0("SRP212166_salmon_", SALMON_DATE))
GTF_PATH    <- "/Users/gricey/Desktop/Internship/boundary_analysis/data/transcriptome_productivity.gtf"
OUTPUT_DIR  <- file.path(BASE_DIR, "results", "normalisation_GSE133398_stageR")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# SAMPLE INFO — same convention as Option A
# =============================================================================
sample_ids <- c(
  "SRR9606472", "SRR9606473", "SRR9606474", "SRR9606475",  # GFP control = "WT"
  "SRR9606468", "SRR9606469", "SRR9606470", "SRR9606471"   # CRE SRSF3 KO = "KO"
)
sample_names <- c("WT_a", "WT_b", "WT_c", "WT_d", "KO_a", "KO_b", "KO_c", "KO_d")
conditions   <- c("WT", "WT", "WT", "WT", "KO", "KO", "KO", "KO")

col_data <- data.frame(
  sample    = sample_names,
  condition = factor(conditions, levels = c("WT", "KO")),
  row.names = sample_names
)

quant_files <- file.path(SALMON_DIR, sample_ids, "quant.sf")
names(quant_files) <- sample_names
missing <- quant_files[!file.exists(quant_files)]
if (length(missing) > 0) stop("Missing quant.sf files:\n", paste(missing, collapse="\n"))
cat("All 8 quant.sf files found.\n")

# =============================================================================
# BUILD tx2gene FROM GTF
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
# tximport (transcript level)
# =============================================================================
cat("Running tximport (transcript level)...\n")
txi_tx <- tximport(quant_files, type = "salmon", txOut = TRUE)
cat("Transcripts imported:", nrow(txi_tx$counts), "\n")

# =============================================================================
# DESeq2 at transcript level
# =============================================================================
cat("\nRunning DESeq2 (transcript level)...\n")
dds_tx <- DESeqDataSetFromTximport(txi_tx, colData = col_data, design = ~ condition)

tx_to_gene              <- setNames(tx2gene$gene_id, tx2gene$tx_id)
rowData(dds_tx)$groupID <- tx_to_gene[rownames(dds_tx)]

no_mapping <- is.na(rowData(dds_tx)$groupID)
cat("Transcripts with no GTF gene mapping (removed):", sum(no_mapping), "\n")
dds_tx <- dds_tx[!no_mapping, ]

keep_tx <- rowSums(counts(dds_tx)) > 20
dds_tx  <- dds_tx[keep_tx, ]
cat("Transcripts after DESeq2 filter (rowSums > 20):", nrow(dds_tx), "\n")

dds_tx <- DESeq(dds_tx)
res_tx <- results(dds_tx, contrast = c("condition", "KO", "WT"),
                  independentFiltering = TRUE, alpha = 0.05)
res_tx$groupID <- rowData(dds_tx)$groupID

cat("\n--- DESeq2 raw results diagnostics ---\n")
cat("Total transcripts tested       :", nrow(res_tx), "\n")
cat("Transcripts with non-NA padj   :", sum(!is.na(res_tx$padj)), "\n")
cat("Transcripts with padj < 0.05   :", sum(res_tx$padj < 0.05, na.rm=TRUE), "\n")
cat("---------------------------------\n\n")

# =============================================================================
# stageR — Bharath's helper functions, verbatim from GSE130613 script
# =============================================================================
perGeneQValueExact <- function(pGene, theta, geneSplit) {
  stopifnot(length(pGene) == length(geneSplit))
  numExons <- listLen(geneSplit)
  tab      <- tabulate(numExons)
  notZero  <- (tab > 0)
  numerator <- mapply(function(m, n) m * (1 - (1 - theta)^n),
                      m = tab[notZero], n = which(notZero))
  numerator <- rowSums(numerator)
  bins  <- cut(pGene, breaks = c(-Inf, as.vector(theta)), right = TRUE, include.lowest = TRUE)
  counts <- tabulate(bins, nbins = nlevels(bins))
  denom  <- cumsum(counts)
  return(numerator / denom)
}

perGeneQValueHack <- function(object, p = "pvalue", method = perGeneQValueExact) {
  wTest    <- which(!is.na(object$padj))
  pvals    <- object[[p]][wTest]
  geneID   <- factor(object[["groupID"]][wTest])
  geneSplit <- split(seq_along(geneID), geneID)
  pGene    <- sapply(geneSplit, function(i) min(pvals[i]))
  stopifnot(all(is.finite(pGene)))
  theta <- unique(sort(pGene))
  q     <- method(pGene, theta, geneSplit)
  res   <- q[match(pGene, theta)]
  res   <- pmin(1, res)
  names(res) <- names(geneSplit)
  stopifnot(!any(is.na(res)))
  return(res)
}

cat("Running stageR stage 1 (gene-level screening)...\n")
pScreen <- perGeneQValueHack(res_tx)
cat("Genes with pScreen < 0.05:", sum(pScreen < 0.05, na.rm=TRUE), "\n")

cat("Running stageR stage 2 (transcript-level confirmation)...\n")
pConfirmation <- matrix(res_tx$pvalue, ncol = 1)
dimnames(pConfirmation) <- list(rownames(res_tx), "transcript")
tx2gene_sub <- tx2gene[tx2gene$tx_id %in% rownames(res_tx), ]

stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,
                       pScreenAdjusted = TRUE, tx2gene = tx2gene_sub)
stageRObj <- stageWiseAdjustment(stageRObj, method = "dte", alpha = 0.05, allowNA = TRUE)

suppressWarnings({
  padj_stageR <- getAdjustedPValues(stageRObj, order = FALSE, onlySignificantGenes = FALSE)
})
cat("Values < 0.05 in gene col:", sum(padj_stageR[,"gene"] < 0.05, na.rm=TRUE), "\n")
cat("Values < 0.05 in tx col:  ", sum(padj_stageR[,"transcript"] < 0.05, na.rm=TRUE), "\n")

# =============================================================================
# Build output table
# =============================================================================
gene_name_map_by_gene <- setNames(
  sapply(gtf_tx$attributes, extract_attr, key = "gene_name"),
  sapply(gtf_tx$attributes, extract_attr, key = "gene_id")
)

res_df               <- as.data.frame(res_tx)
res_df$transcript_id <- rownames(res_df)
res_df$gene_id       <- res_df$groupID
res_df$gene_name     <- gene_name_map_by_gene[res_df$gene_id]
res_df$gene_name[is.na(res_df$gene_name)] <- "unannotated"

padj_df <- padj_stageR
colnames(padj_df)[colnames(padj_df) == "txID"]       <- "transcript_id"
colnames(padj_df)[colnames(padj_df) == "gene"]        <- "padj_gene_stageR"
colnames(padj_df)[colnames(padj_df) == "transcript"]  <- "padj_tx_stageR"

res_df <- merge(res_df, padj_df[, c("transcript_id", "padj_gene_stageR", "padj_tx_stageR")],
                by = "transcript_id", all.x = TRUE)

out_df <- res_df[, c("transcript_id", "gene_id", "gene_name",
                      "log2FoldChange", "pvalue", "padj_gene_stageR", "padj_tx_stageR")]
out_df <- out_df[order(out_df$padj_gene_stageR, out_df$padj_tx_stageR, na.last = TRUE), ]

write.table(out_df, file = file.path(OUTPUT_DIR, "deseq2_stageR_KO_vs_WT.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

sig_genes <- length(unique(out_df$gene_id[!is.na(out_df$padj_gene_stageR) & out_df$padj_gene_stageR < 0.05]))
sig_tx    <- sum(!is.na(out_df$padj_tx_stageR) & out_df$padj_tx_stageR < 0.05)
cat("\nSignificant genes (stageR stage 1, padj<0.05)      :", sig_genes, "\n")
cat("Significant transcripts (stageR stage 2, padj<0.05):", sig_tx, "\n")
cat("Saved: deseq2_stageR_KO_vs_WT.tsv\n")
cat("\n=== DONE (Option B) ===\n")
