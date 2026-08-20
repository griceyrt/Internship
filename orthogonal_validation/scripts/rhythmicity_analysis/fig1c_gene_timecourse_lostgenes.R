#!/usr/bin/env Rscript
# Fig 1c style, GENE level, restricted to the 43 "lost genes" -- same gene
# list/order as fig1c_isoform_timecourse_lostgenes.R, one line per gene.
# Author: Gricey
# Run from orthogonal_validation/: Rscript scripts/fig1c_gene_timecourse_lostgenes.R

library(tximport)
library(DESeq2)
library(stringr)
library(dplyr)
library(tidyr)
library(tibble)
library(limma)
library(ggplot2)

salmon_path <- "data/salmon_2024-09-25/"
gtf_path <- "data/transcriptome_productivity.gtf"
out_dir <- "results/rhythmicity_analysis/gene_rhythmicity/"

genes_per_page <- 12

gtf_lines <- readLines(gtf_path)
tx_lines <- gtf_lines[grepl("\ttranscript\t", gtf_lines)]
gene_id <- str_match(tx_lines, 'gene_id "([^"]+)"')[, 2]
transcript_id <- str_match(tx_lines, 'transcript_id "([^"]+)"')[, 2]
gene_name <- str_match(tx_lines, 'gene_name "([^"]+)"')[, 2]
tx2gene <- data.frame(transcript_id = transcript_id, gene_id = gene_id) %>% distinct()
gene2name <- data.frame(gene_id = gene_id, gene_name = gene_name) %>%
  filter(!is.na(gene_name)) %>% distinct(gene_id, .keep_all = TRUE)

gene_list <- read.csv(paste(out_dir, "csv/lost_genes_list.csv", sep = ""))$gene_id
length(gene_list)  # should be 43

files <- list.files(path = salmon_path, pattern = "\\.sf$", full.names = TRUE, recursive = TRUE)
names(files) <- str_extract(files, "([CKTO]+\\d+_\\d)(?=\\.sf$)")
wt_files <- files[grepl("^CT", names(files))]

txi <- tximport(wt_files, type = "salmon", tx2gene = tx2gene, txOut = FALSE, countsFromAbundance = "no",
                 ignoreTxVersion = TRUE, dropInfReps = TRUE, importer = function(x) read.delim(x))

exp_design <- data.frame(sample_id = colnames(txi$counts),
                         batch = factor(str_extract(colnames(txi$counts), "(?<=_)\\d")),
                         time = as.numeric(str_extract(colnames(txi$counts), "\\d+"))) %>%
  mutate(inphase = cos(2 * pi * time / 24), outphase = sin(2 * pi * time / 24))
rownames(exp_design) <- exp_design$sample_id

dds <- DESeqDataSetFromTximport(txi, colData = exp_design, design = ~ batch + inphase + outphase)
dds <- dds[rownames(dds) %in% gene_list, ]
dds <- DESeq(dds, reduced = ~batch, test = "LRT", fitType = "local")

abundance <- varianceStabilizingTransformation(dds, fitType = "local")
abundance_no_batch <- limma::removeBatchEffect(assay(abundance), batch = colData(dds)$batch,
  design = model.matrix(~inphase + outphase, data = colData(dds)))

df <- abundance_no_batch %>%
  as.data.frame() %>% rownames_to_column("gene_id") %>%
  pivot_longer(-gene_id, names_to = c("sample", "rep"), names_sep = "_", values_to = "value") %>%
  separate(sample, c("cond", "time"), convert = TRUE, sep = 2) %>%
  dplyr::select(-cond) %>%
  left_join(gene2name, by = "gene_id")

genes_ordered <- gene2name %>% filter(gene_id %in% gene_list) %>% arrange(gene_name) %>% pull(gene_name)
length(genes_ordered)

pdf(paste(out_dir, "figures/fig1c_gene_timecourse_lostgenes.pdf", sep = ""), width = 10, height = 7.5)
chunks <- split(genes_ordered, ceiling(seq_along(genes_ordered) / genes_per_page))
for (chunk in chunks) {
  p <- ggplot(data = filter(df, gene_name %in% chunk)) +
    geom_point(aes(x = time, y = value), size = 0.5, color = "black") +
    geom_smooth(aes(x = time, y = value, group = gene_id), linewidth = 0.7,
                formula = "y ~ x", method = "loess", se = FALSE, color = "black") +
    facet_wrap(~gene_name, scales = "free_y", ncol = 4, nrow = 3) +
    scale_x_continuous(breaks = seq(0, 20, 4)) +
    xlab("circadian time") + ylab(bquote("log"[2] ~ "gene expression")) +
    theme_classic(base_size = 10) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          strip.text = element_text(size = 9), strip.background = element_rect(fill = "white", color = "white"),
          aspect.ratio = 0.85)
  print(p)
}
dev.off()
