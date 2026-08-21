#!/usr/bin/env Rscript
# =============================================================================
# GSE188840 -- Normalisation & DESeq2 (gene-level, oarfish quantification)
#
# Dataset: WT vs Fmr1 KO mouse brain cortex, ONT direct RNA-seq, 3 reps each.
# DESeq2 runs at gene-level (SUPPA's own TPM input is unaffected by this
# choice, since it comes from the quantifier directly, not from DESeq2).
#
# Step 5: tximport (transcript-level) -> SUPPA TPM export
# Step 6: tximport (gene-level)       -> DESeq2 (KO vs WT)
# =============================================================================

# =============================================================================
# LIBRARIES
# =============================================================================
library(tximport)
library(DESeq2)

# =============================================================================
# PATHS
# =============================================================================
BASE_DIR    <- "/home/grangel/fmr1_polya_splicing"
OARFISH_DIR <- file.path(BASE_DIR, "results", "gse188840_oarfish_2026-08-04")
GTF_PATH    <- file.path(BASE_DIR, "data", "reference", "Mus_musculus.GRCm39.116.gtf.gz")
OUTPUT_DIR  <- file.path(BASE_DIR, "results", "normalisation")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# SAMPLE INFO
# =============================================================================
sample_names <- c("WT_Rep1", "WT_Rep2", "WT_Rep3",
                   "KO_Rep1", "KO_Rep2", "KO_Rep3")
conditions   <- c("WT", "WT", "WT", "KO", "KO", "KO")

col_data <- data.frame(
  sample    = sample_names,
  condition = factor(conditions, levels = c("WT", "KO")),
  row.names = sample_names
)

# oarfish quant file paths
quant_files <- file.path(OARFISH_DIR, paste0(sample_names, ".quant"))
names(quant_files) <- sample_names

missing <- quant_files[!file.exists(quant_files)]
if (length(missing) > 0) stop("Missing .quant files:\n", paste(missing, collapse="\n"))
cat("All 6 .quant files found.\n")

# =============================================================================
# BUILD tx2gene FROM mm39 GTF (gzipped)
# =============================================================================
cat("Building tx2gene from GTF...\n")
gtf <- read.table(gzfile(GTF_PATH), sep="\t", quote="", comment.char="#",
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
# STEP 5 -- tximport (transcript-level, for SUPPA)
# =============================================================================
cat("Running tximport (transcript level, oarfish)...\n")
# ignoreTxVersion only strips versions when matching against a tx2gene
# table -- with txOut=TRUE and no tx2gene here, it has nothing to strip
# against, so version suffixes are stripped manually below instead (SUPPA's
# GTF-based .ioe files use non-versioned transcript IDs).
txi_tx <- tximport(quant_files,
                   type  = "oarfish",
                   txOut = TRUE,
                   ignoreTxVersion = TRUE)

# SUPPA needs TPM (comparable across samples), not raw counts -- tximport's
# $abundance is TPM regardless of which quantifier produced the input
tpm_tx <- txi_tx$abundance

# Explicitly strip Ensembl version suffixes (e.g. ".1", ".2") from the
# transcript IDs so they match the GTF's non-versioned transcript_id
# convention that SUPPA's .ioe event files use.
rownames(tpm_tx) <- sub("\\.[0-9]+$", "", rownames(tpm_tx))

cat("Transcripts imported:", nrow(tpm_tx), "\n")

# =============================================================================
# STEP 5b -- Export TPM matrix for SUPPA
# =============================================================================
cat("Exporting SUPPA TPM matrix...\n")
suppa_table <- as.data.frame(tpm_tx)
colnames(suppa_table) <- sample_names

write.table(suppa_table,
            file      = file.path(OUTPUT_DIR, "TPM_all_samples.tpm"),
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = TRUE)
cat("Saved: TPM_all_samples.tpm\n")

# =============================================================================
# STEP 6 -- DESeq2 (gene-level, per Lies's call)
# Contrast: KO vs WT
# =============================================================================
cat("\nRunning tximport (gene level, oarfish)...\n")
txi_gene <- tximport(quant_files,
                     type    = "oarfish",
                     tx2gene = tx2gene,
                     txOut   = FALSE,
                     ignoreTxVersion = TRUE)

cat("Running DESeq2 (gene-level)...\n")
dds <- DESeqDataSetFromTximport(txi_gene,
                                colData = col_data,
                                design  = ~ condition)
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "KO", "WT"))
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

write.table(res_df,
            file      = file.path(OUTPUT_DIR, "deseq2_KO_vs_WT.tsv"),
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

sig <- sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0.5, na.rm = TRUE)
cat("DESeq2 significant genes (|log2FC|>0.5, padj<0.05):", sig, "\n")
cat("Saved: deseq2_KO_vs_WT.tsv\n")

cat("\n=== DONE ===\n")
cat("SUPPA input table :", file.path(OUTPUT_DIR, "TPM_all_samples.tpm"), "\n")
cat("DESeq2 results    :", file.path(OUTPUT_DIR, "deseq2_KO_vs_WT.tsv"), "\n")
