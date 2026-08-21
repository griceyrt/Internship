#!/usr/bin/env Rscript
# =============================================================================
# 2+3 replicate test: gene-level DESeq2 with WT_Rep1 excluded (2 WT + 3 KO),
# matching the paper's replicate structure.
#
# Built after qc_replicate_outlier_check.R confirmed WT_Rep1 as an outlier
# on all three checks (PCA, sample-to-sample distance, Cook's distance) --
# see results/qc_replicate_check/. Tests whether excluding it explains part
# of the gap to the paper's numbers.
#
# Kept in its own folder (de_2plus3_test/), separate from de_2x2_test/ --
# a different axis (replicate count) from the level x threshold x fold-
# change grid tested there. Same tximport/DESeq2 pattern as
# normalisation_deseq2_oarfish.R, with WT_Rep1 dropped from the sample list.
# =============================================================================

library(tximport)
library(DESeq2)

BASE_DIR    <- "/home/grangel/fmr1_polya_splicing"
OARFISH_DIR <- file.path(BASE_DIR, "results", "gse188840_oarfish_2026-08-04")
GTF_PATH    <- file.path(BASE_DIR, "data", "reference", "Mus_musculus.GRCm39.116.gtf.gz")
OUTPUT_DIR  <- file.path(BASE_DIR, "results", "de_2plus3_test")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# WT_Rep1 EXCLUDED here -- the only difference from the main pipeline's
# sample list (normalisation_deseq2_oarfish.R uses all 6)
sample_names <- c("WT_Rep2", "WT_Rep3", "KO_Rep1", "KO_Rep2", "KO_Rep3")
conditions   <- c("WT", "WT", "KO", "KO", "KO")

col_data <- data.frame(
  sample    = sample_names,
  condition = factor(conditions, levels = c("WT", "KO")),
  row.names = sample_names
)

quant_files <- file.path(OARFISH_DIR, paste0(sample_names, ".quant"))
names(quant_files) <- sample_names
missing <- quant_files[!file.exists(quant_files)]
if (length(missing) > 0) stop("Missing .quant files:\n", paste(missing, collapse = "\n"))
cat("Using", length(sample_names), "samples (WT_Rep1 excluded):",
    paste(sample_names, collapse = ", "), "\n")

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
tx2gene <- data.frame(tx_id = sapply(gtf_tx$attributes, extract_attr, key = "transcript_id"),
                       gene_id = sapply(gtf_tx$attributes, extract_attr, key = "gene_id"),
                       stringsAsFactors = FALSE, row.names = NULL)
tx2gene <- tx2gene[!is.na(tx2gene$tx_id) & !is.na(tx2gene$gene_id), ]

cat("Running tximport (gene level, oarfish, 2+3 samples)...\n")
txi_gene <- tximport(quant_files, type = "oarfish", tx2gene = tx2gene,
                      txOut = FALSE, ignoreTxVersion = TRUE)

cat("Running DESeq2 (gene-level, WT_Rep1 excluded)...\n")
dds <- DESeqDataSetFromTximport(txi_gene, colData = col_data, design = ~condition)
dds <- DESeq(dds)
saveRDS(dds, file.path(OUTPUT_DIR, "dds_gene_level_2plus3.rds"))

res <- results(dds, contrast = c("condition", "KO", "WT"))
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

write.table(res_df,
            file      = file.path(OUTPUT_DIR, "deseq2_gene_level_2plus3_KO_vs_WT.tsv"),
            sep       = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== DONE ===\n")
cat("padj<0.05  :", sum(res_df$padj < 0.05, na.rm = TRUE), "genes\n")
cat("raw p<0.05 :", sum(res_df$pvalue < 0.05, na.rm = TRUE), "genes\n")
cat("Saved:", file.path(OUTPUT_DIR, "deseq2_gene_level_2plus3_KO_vs_WT.tsv"), "\n")
