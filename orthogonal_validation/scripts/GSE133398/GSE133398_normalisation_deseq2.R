#!/usr/bin/env Rscript
# GSE133398 -- SRSF3 KO vs control hepatocytes: tximport/filter/normalize for
# SUPPA, plus gene-level DESeq2. Adapted from normalisation_deseq2.R (GSE130613).
# Author: Gricey
# Run from: orthogonal_validation/

library(tximport)
library(DESeq2)
library(limma)

BASE_DIR    <- "/Users/gricey/Desktop/Internship/orthogonal_validation"
SALMON_DATE <- "2026-07-08"   # <-- UPDATE to match the date in your results folder name
SALMON_DIR  <- file.path(BASE_DIR, "results", "GSE133398", paste0("SRP212166_salmon_", SALMON_DATE))
GTF_PATH    <- "/Users/gricey/Desktop/Internship/boundary_analysis/data/transcriptome_productivity.gtf"
OUTPUT_DIR  <- file.path(BASE_DIR, "results", "GSE133398", "normalisation_GSE133398")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# WT = GFP control, KO = CRE / SRSF3 KO -- order matches run_suppa.sh's split step
sample_ids <- c(
  "SRR9606472", "SRR9606473", "SRR9606474", "SRR9606475",
  "SRR9606468", "SRR9606469", "SRR9606470", "SRR9606471"
)
sample_names <- c("WT_a", "WT_b", "WT_c", "WT_d",
                  "KO_a", "KO_b", "KO_c", "KO_d")
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

cat("Running tximport (transcript level)...\n")
txi_tx <- tximport(quant_files, type = "salmon", txOut = TRUE)
raw_counts_tx <- txi_tx$counts
cat("Transcripts imported:", nrow(raw_counts_tx), "\n")

cat("Filtering low-expressed transcripts...\n")
keep <- rowSums(raw_counts_tx > 5) >= 2
counts_filtered <- raw_counts_tx[keep, ]
cat("Transcripts after filtering:", nrow(counts_filtered),
    "(removed:", sum(!keep), ")\n")

cat("Normalizing to total reads per sample...\n")
total_reads <- colSums(counts_filtered)
counts_norm <- sweep(counts_filtered, 2, total_reads, FUN = "/")

APPLY_BATCH_CORRECTION <- FALSE
# batch_vector <- c(...)  # define if needed

if (APPLY_BATCH_CORRECTION) {
  cat("Applying limma removeBatchEffect...\n")
  counts_suppa <- removeBatchEffect(counts_norm, batch = batch_vector)
} else {
  cat("Skipping batch correction (no batch info available).\n")
  counts_suppa <- counts_norm
}

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

cat("\nRunning DESeq2 (gene-level)...\n")

txi_gene <- tximport(quant_files,
                     type    = "salmon",
                     tx2gene = tx2gene,
                     txOut   = FALSE)

dds <- DESeqDataSetFromTximport(txi_gene,
                                colData = col_data,
                                design  = ~ condition)
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "KO", "WT"))
res_df <- as.data.frame(res)
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

write.table(res_df,
            file      = file.path(OUTPUT_DIR, "deseq2_KO_vs_WT.tsv"),
            sep       = "\t",
            quote     = FALSE,
            col.names = NA)

sig <- sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0.5, na.rm = TRUE)
cat("DESeq2 significant genes (|log2FC|>0.5, padj<0.05):", sig, "\n")
cat("Saved: deseq2_KO_vs_WT.tsv\n")

cat("\n=== DONE ===\n")
cat("SUPPA input table : results/GSE133398/normalisation_GSE133398/combined_norm.tab\n")
cat("DESeq2 results    : results/GSE133398/normalisation_GSE133398/deseq2_KO_vs_WT.tsv\n")
cat("\nNOTE: WT = GFP control (SRR9606472-75), KO = CRE / SRSF3 KO (SRR9606468-71)\n")
