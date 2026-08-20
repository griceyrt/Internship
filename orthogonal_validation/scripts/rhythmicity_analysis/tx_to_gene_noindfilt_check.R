#!/usr/bin/env Rscript
# Sanity check -- same whole-transcriptome transcript-level test as
# tx_to_gene_rhythmicity.R (906 sig transcripts, independentFiltering=TRUE,
# no stageR), rerun with independentFiltering=FALSE for comparison.
# Author: Gricey
# Run from orthogonal_validation/: Rscript scripts/tx_to_gene_noindfilt_check.R

library(tximport)
library(DESeq2)
library(stringr)
library(dplyr)

salmon_path <- "data/salmon_2024-09-25/"
gtf_path <- "data/transcriptome_productivity.gtf"
out_dir <- "results/rhythmicity_analysis/gene_rhythmicity/"

count_threshold <- 5
min_samples <- 2
fdr_cutoff <- 0.05

gtf_lines <- readLines(gtf_path)
tx_lines <- gtf_lines[grepl("\ttranscript\t", gtf_lines)]
gene_id <- str_match(tx_lines, 'gene_id "([^"]+)"')[, 2]
transcript_id <- str_match(tx_lines, 'transcript_id "([^"]+)"')[, 2]
tx2gene <- data.frame(transcript_id = transcript_id, gene_id = gene_id) %>% distinct()

files <- list.files(path = salmon_path, pattern = "\\.sf$", full.names = TRUE, recursive = TRUE)
names(files) <- str_extract(files, "([CKTO]+\\d+_\\d)(?=\\.sf$)")
wt_files <- files[grepl("^CT", names(files))]

txi_tx_wt <- tximport(wt_files, type = "salmon", txOut = TRUE, countsFromAbundance = "no",
                       ignoreTxVersion = TRUE, dropInfReps = TRUE, importer = function(x) read.delim(x))
counts_tx_wt <- txi_tx_wt$counts
mode(counts_tx_wt) <- "integer"
expressed_tx <- rownames(counts_tx_wt)[rowSums(counts_tx_wt > count_threshold) >= min_samples]
length(expressed_tx)

exp_design <- data.frame(sample_id = colnames(txi_tx_wt$counts),
                         batch = factor(str_extract(colnames(txi_tx_wt$counts), "(?<=_)\\d")),
                         time = as.numeric(str_extract(colnames(txi_tx_wt$counts), "\\d+"))) %>%
  mutate(inphase = cos(2 * pi * time / 24), outphase = sin(2 * pi * time / 24))
rownames(exp_design) <- exp_design$sample_id

dds_tx <- DESeqDataSetFromTximport(txi_tx_wt, colData = exp_design, design = ~ batch + inphase + outphase)
dds_tx <- dds_tx[rownames(dds_tx) %in% expressed_tx, ]
dds_tx <- DESeq(dds_tx, reduced = ~batch, test = "LRT", fitType = "local")

res_indfilt_on <- DESeq2::results(dds_tx, independentFiltering = TRUE, alpha = fdr_cutoff, tidy = TRUE)
res_indfilt_off <- DESeq2::results(dds_tx, independentFiltering = FALSE, alpha = fdr_cutoff, tidy = TRUE)

sig_on <- res_indfilt_on %>% filter(!is.na(padj) & padj < fdr_cutoff) %>%
  left_join(tx2gene, by = c("row" = "transcript_id"))
sig_off <- res_indfilt_off %>% filter(!is.na(padj) & padj < fdr_cutoff) %>%
  left_join(tx2gene, by = c("row" = "transcript_id"))

cat("independentFiltering=TRUE  -- significant transcripts:", nrow(sig_on), " unique genes:", n_distinct(sig_on$gene_id), "\n")
cat("independentFiltering=FALSE -- significant transcripts:", nrow(sig_off), " unique genes:", n_distinct(sig_off$gene_id), "\n")
