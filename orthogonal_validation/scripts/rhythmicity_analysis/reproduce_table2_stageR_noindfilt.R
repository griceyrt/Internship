#!/usr/bin/env Rscript
# Sensitivity check -- identical to reproduce_table2_stageR.R except
# independentFiltering=FALSE in the DESeq2::results() call, to see whether
# 501 is sensitive to that choice.
# Author: Gricey
# Run from orthogonal_validation/: Rscript scripts/reproduce_table2_stageR_noindfilt.R

library(tximport)
library(DESeq2)
library(stringr)
library(dplyr)
library(DEXSeq)
library(stageR)

salmon_path <- "data/salmon_2024-09-25/"
gtf_path <- "data/transcriptome_productivity.gtf"
out_dir <- "results/rhythmicity_analysis/gene_rhythmicity/"

count_threshold <- 5
min_samples <- 2
fdr_cutoff <- 0.05

perGeneQValueHack <- function(object, p = "pvalue", method = DEXSeq:::perGeneQValueExact) {
  wTest <- which(!is.na(object$padj))
  pvals <- object[[p]][wTest]
  geneID <- factor(object[["groupID"]][wTest])
  geneSplit <- split(seq(along = geneID), geneID)
  pGene <- sapply(geneSplit, function(i) min(pvals[i]))
  stopifnot(all(is.finite(pGene)))
  theta <- unique(sort(pGene))
  q <- method(pGene, theta, geneSplit)
  res <- rep(NA_real_, length(pGene))
  res <- q[match(pGene, theta)]
  res <- pmin(1, res)
  names(res) <- names(geneSplit)
  stopifnot(!any(is.na(res)))
  return(res)
}

gtf_lines <- readLines(gtf_path)
tx_lines <- gtf_lines[grepl("\ttranscript\t", gtf_lines)]
gene_id <- str_match(tx_lines, 'gene_id "([^"]+)"')[, 2]
transcript_id <- str_match(tx_lines, 'transcript_id "([^"]+)"')[, 2]
tx2gene <- data.frame(transcript_id = transcript_id, gene_id = gene_id) %>% distinct()

files <- list.files(path = salmon_path, pattern = "\\.sf$", full.names = TRUE, recursive = TRUE)
names(files) <- str_extract(files, "([CKTO]+\\d+_\\d)(?=\\.sf$)")
wt_files <- files[grepl("^CT", names(files))]

txi_all <- tximport(files, type = "salmon", txOut = TRUE, countsFromAbundance = "no",
                     ignoreTxVersion = TRUE, dropInfReps = TRUE, importer = function(x) read.delim(x))
counts_all <- txi_all$counts
mode(counts_all) <- "integer"
expressed_tx <- rownames(counts_all)[rowSums(counts_all > count_threshold) >= min_samples]
length(expressed_tx)

txi <- tximport(wt_files, type = "salmon", txOut = TRUE, countsFromAbundance = "no",
                 ignoreTxVersion = TRUE, dropInfReps = TRUE, importer = function(x) read.delim(x))

exp_design <- data.frame(sample_id = colnames(txi$counts),
                         batch = factor(str_extract(colnames(txi$counts), "(?<=_)\\d")),
                         time = as.numeric(str_extract(colnames(txi$counts), "\\d+"))) %>%
  mutate(inphase = cos(2 * pi * time / 24), outphase = sin(2 * pi * time / 24))
rownames(exp_design) <- exp_design$sample_id

counts <- txi$counts
mode(counts) <- "integer"

dds <- DESeqDataSetFromMatrix(countData = counts, colData = exp_design, design = ~ batch + inphase + outphase)
dds <- dds[rownames(dds) %in% expressed_tx, ]
dds <- DESeq(dds, reduced = ~batch, test = "LRT", fitType = "local")

res <- DESeq2::results(dds, independentFiltering = FALSE, alpha = fdr_cutoff, tidy = TRUE)
res <- res %>% dplyr::rename(transcript_id = row) %>% left_join(tx2gene, by = "transcript_id")
res$groupID <- res$gene_id

pScreen <- perGeneQValueHack(res)
pConfirmation <- matrix(res$pvalue, ncol = 1)
dimnames(pConfirmation) <- list(res$transcript_id, "transcript")

stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,
                       pScreenAdjusted = TRUE, tx2gene = tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method = "dte", alpha = fdr_cutoff, allowNA = TRUE)
circadian_analysis <- getAdjustedPValues(stageRObj, order = FALSE, onlySignificantGenes = FALSE)

write.csv(circadian_analysis, paste(out_dir, "csv/reproduce_table2_stageR_noindfilt.csv", sep = ""), row.names = FALSE)

sig_genes <- unique(circadian_analysis$geneID[circadian_analysis$gene < fdr_cutoff])
length(sig_genes)  # compare to the independentFiltering=TRUE version's 501
