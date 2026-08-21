#!/usr/bin/env Rscript
# =============================================================================
# 2x2 test: isoform-level DESeq2, to compare against Figure 3F of Shin/
# Chikhaoui et al. 2022 (the GSE188840 source paper, which used isoform-
# level DESeq2, vs. our main pipeline's gene-level). Holds mm39 + oarfish
# fixed so only level (and, in downstream plotting, threshold type) varies.
#
# Reuses the transcript-level tximport import already built in
# normalisation_deseq2_oarfish.R (that script only ever exported it as a
# TPM matrix for SUPPA) and adds the missing DESeq2 run on top of it.
# =============================================================================

library(tximport)
library(DESeq2)

BASE_DIR    <- "/home/grangel/fmr1_polya_splicing"
OARFISH_DIR <- file.path(BASE_DIR, "results", "gse188840_oarfish_2026-08-04")
GTF_PATH    <- file.path(BASE_DIR, "data", "reference", "Mus_musculus.GRCm39.116.gtf.gz")
OUTPUT_DIR  <- file.path(BASE_DIR, "results", "de_2x2_test")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

sample_names <- c("WT_Rep1", "WT_Rep2", "WT_Rep3", "KO_Rep1", "KO_Rep2", "KO_Rep3")
conditions   <- c("WT", "WT", "WT", "KO", "KO", "KO")

col_data <- data.frame(
  sample    = sample_names,
  condition = factor(conditions, levels = c("WT", "KO")),
  row.names = sample_names
)

quant_files <- file.path(OARFISH_DIR, paste0(sample_names, ".quant"))
names(quant_files) <- sample_names
missing <- quant_files[!file.exists(quant_files)]
if (length(missing) > 0) stop("Missing .quant files:\n", paste(missing, collapse = "\n"))
cat("All 6 .quant files found.\n")

# =============================================================================
# tx2gene from mm39 GTF, also saved to file here (unlike the main pipeline
# script) since the isoform-level Panel-A bar chart needs it to join each
# transcript back to its parent gene's biotype.
# =============================================================================
cat("Building tx2gene from GTF...\n")
gtf <- read.table(gzfile(GTF_PATH), sep = "\t", quote = "", comment.char = "#",
                   col.names = c("chr", "source", "feature", "start", "end",
                                 "score", "strand", "frame", "attributes"))
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

write.table(tx2gene, file = file.path(OUTPUT_DIR, "tx2gene_mm39.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# =============================================================================
# Isoform-level tximport + DESeq2 (the new part -- gene-level equivalent
# already exists in results/normalisation/deseq2_KO_vs_WT.tsv)
# =============================================================================
cat("Running tximport (transcript level, oarfish)...\n")
txi_tx <- tximport(quant_files, type = "oarfish", txOut = TRUE, ignoreTxVersion = TRUE)
# same version-suffix fix as the gene-level script, needed here too since
# DESeq2 will use these rownames directly as transcript IDs
rownames(txi_tx$counts) <- sub("\\.[0-9]+$", "", rownames(txi_tx$counts))
rownames(txi_tx$abundance) <- sub("\\.[0-9]+$", "", rownames(txi_tx$abundance))
rownames(txi_tx$length) <- sub("\\.[0-9]+$", "", rownames(txi_tx$length))

cat("Running DESeq2 (isoform/transcript-level)...\n")
dds_tx <- DESeqDataSetFromTximport(txi_tx, colData = col_data, design = ~ condition)
dds_tx <- DESeq(dds_tx)

res_tx <- results(dds_tx, contrast = c("condition", "KO", "WT"))
res_tx_df <- as.data.frame(res_tx)
res_tx_df$transcript_id <- rownames(res_tx_df)
res_tx_df <- merge(res_tx_df, tx2gene, by.x = "transcript_id", by.y = "tx_id", all.x = TRUE)
res_tx_df <- res_tx_df[order(res_tx_df$padj, na.last = TRUE), ]

write.table(res_tx_df,
            file = file.path(OUTPUT_DIR, "deseq2_isoform_level_KO_vs_WT.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== DONE ===\n")
cat("Isoform-level DESeq2 results:", file.path(OUTPUT_DIR, "deseq2_isoform_level_KO_vs_WT.tsv"), "\n")
cat("tx2gene table               :", file.path(OUTPUT_DIR, "tx2gene_mm39.tsv"), "\n")
cat("\nSummary (padj<0.05, no log2FC cutoff, matching the gene-level Panel A convention):\n")
cat("  padj<0.05  :", sum(res_tx_df$padj < 0.05, na.rm = TRUE), "isoforms\n")
cat("  raw p<0.05 :", sum(res_tx_df$pvalue < 0.05, na.rm = TRUE), "isoforms\n")
