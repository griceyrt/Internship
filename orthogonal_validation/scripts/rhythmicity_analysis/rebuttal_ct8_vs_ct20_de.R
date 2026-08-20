#!/usr/bin/env Rscript
# Rebuttal point (d) -- WT CT8 vs CT20 differential expression, real R DESeq2.
# Reproduces fig4_1_2024_12_04.R's method (transcript-level tximport, global
# expressed-transcript filter across all 16 samples, plain DESeq2 Wald test,
# fitType="local"), swapping the genotype contrast for a timepoint contrast
# (WT-CT8 vs WT-CT20). Significance: padj<0.05 & |log2FoldChange|>1.5.
# Author: Gricey
# Needs R + Bioconductor (tximport, DESeq2) -- run on the PSMN cluster.
# USAGE: Rscript rebuttal_ct8_vs_ct20_de.R
# Run from: orthogonal_validation/ (or adjust PROJECT_ROOT below)

suppressMessages({
  library(tximport)
  library(DESeq2)
  library(stringr)
  library(dplyr)
  library(tibble)
})

PROJECT_ROOT <- "."  # change to an absolute path if running from elsewhere on the cluster
SALMON_DIR <- file.path(PROJECT_ROOT, "data", "salmon_2024-09-25")
GTF_PATH <- file.path(PROJECT_ROOT, "data", "transcriptome_productivity.gtf")
OUT_DIR <- file.path(PROJECT_ROOT, "results", "rebuttal_ct8_vs_ct20")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

COUNT_THRESHOLD <- 5
MIN_SAMPLES <- 2
FDR_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 1.5

cat("Building tx2gene map from", GTF_PATH, "\n")
gtf_lines <- readLines(GTF_PATH)
tx_lines <- gtf_lines[grepl("\ttranscript\t", gtf_lines)]
gene_id <- str_match(tx_lines, 'gene_id "([^"]+)"')[, 2]
transcript_id <- str_match(tx_lines, 'transcript_id "([^"]+)"')[, 2]
tx2gene <- data.frame(transcript_id = transcript_id, gene_id = gene_id) %>% distinct()
cat(" ", nrow(tx2gene), "transcript -> gene mappings\n")

cat("\nLoading all 16 samples for the global expression filter...\n")
all_files <- list.files(SALMON_DIR, pattern = "\\.sf$", full.names = TRUE, recursive = TRUE)
stopifnot(length(all_files) == 16)
names(all_files) <- str_extract(all_files, "([CKTO]+\\d+_\\d)(?=\\.sf$)")

txi_all <- tximport(all_files, type = "salmon", txOut = TRUE,
                     countsFromAbundance = "no", ignoreTxVersion = TRUE,
                     dropInfReps = TRUE,
                     importer = function(x) read.delim(x))
counts_all <- txi_all$counts
mode(counts_all) <- "integer"

keep <- rowSums(counts_all > COUNT_THRESHOLD) >= MIN_SAMPLES
expressed_tx <- rownames(counts_all)[keep]
cat(" ", length(expressed_tx), "/", nrow(counts_all), "transcripts pass the filter\n")
write.csv(data.frame(transcript_id = expressed_tx),
          file.path(OUT_DIR, "expressed_tx_list.csv"), row.names = FALSE)

cat("\nLoading CT8 + CT20 transcript counts...\n")
ct_files <- all_files[c("CT8_1", "CT8_2", "CT20_1", "CT20_2")]
txi <- tximport(ct_files, type = "salmon", txOut = TRUE,
                 countsFromAbundance = "no", ignoreTxVersion = TRUE,
                 dropInfReps = TRUE,
                 importer = function(x) read.delim(x))

exp_design <- data.frame(
  sample_id = colnames(txi$abundance),
  condition = factor(ifelse(grepl("^CT8", colnames(txi$abundance)), "CT8", "CT20"),
                      levels = c("CT8", "CT20"))
)
rownames(exp_design) <- exp_design$sample_id

counts <- txi$counts
mode(counts) <- "integer"

dds <- DESeqDataSetFromMatrix(countData = counts, colData = exp_design, design = ~condition)
keep <- rownames(counts) %in% expressed_tx
dds <- dds[keep, ]
cat(" ", sum(keep), "transcripts after filter\n")

cat("\nRunning DESeq2 (Wald test, condition CT20 vs CT8)...\n")
dds <- DESeq(dds, fitType = "local")
res <- DESeq2::results(dds, contrast = c("condition", "CT20", "CT8"),
                        independentFiltering = TRUE, alpha = FDR_CUTOFF, tidy = TRUE)

results_df <- res %>%
  dplyr::rename(transcript_id = row) %>%
  left_join(tx2gene, by = "transcript_id")

write.csv(results_df, file.path(OUT_DIR, "deseq2_R_tx_level_ct8_vs_ct20.csv"), row.names = FALSE)

results_valid <- results_df %>% filter(!is.na(padj))
sig <- results_valid %>% filter(abs(log2FoldChange) > LOG2FC_CUTOFF & padj < FDR_CUTOFF)

cat("\nTranscripts with valid padj:", nrow(results_valid), "\n")
cat(sprintf("Significant (|log2FC|>%s & padj<%s): %d transcripts, %d unique genes\n",
            LOG2FC_CUTOFF, FDR_CUTOFF, nrow(sig), n_distinct(sig$gene_id)))
cat("\nSaved:", file.path(OUT_DIR, "deseq2_R_tx_level_ct8_vs_ct20.csv"), "\n")
