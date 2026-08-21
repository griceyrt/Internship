#!/usr/bin/env Rscript
# Fig 1c style, OVERLAP -- isoform-level (coding=black, non-coding=orange) and
# gene-level (blue) curves on the SAME panel per gene, for the 228 genes in
# Table5-rhythmic_AS_events.xlsx (rhythmic AS events, DEXSeq+DRIMSeq+stageR).
# Filtered version of ../fig1c_overlap_timecourse_dtugenes.R -- restricts the
# isoform side (tx_ids) to expressed transcripts (>5 reads in >=2 WT samples),
# same definition as tx_to_gene_rhythmicity.R (Khushi, 2026-08-20 email).
# Author: Gricey
# Run from orthogonal_validation/: Rscript scripts/rhythmicity_analysis/rhythmicity_filtered/fig1c_overlap_timecourse_dtugenes.R

library(tximport)
library(DESeq2)
library(stringr)
library(dplyr)
library(tidyr)
library(tibble)
library(limma)
library(ggplot2)
library(ggthemes)
library(ggnewscale)

salmon_path <- "data/salmon_2024-09-25/"
gtf_path <- "data/transcriptome_productivity.gtf"
productivity_path <- "data/files_from_khushi/productivity.txt"
in_dir <- "results/rhythmicity_analysis/gene_rhythmicity/"
out_dir <- "results/rhythmicity_analysis/rhythmicity_filtered/"

genes_per_page <- 12

gtf_lines <- readLines(gtf_path)
tx_lines <- gtf_lines[grepl("\ttranscript\t", gtf_lines)]
gene_id <- str_match(tx_lines, 'gene_id "([^"]+)"')[, 2]
transcript_id <- str_match(tx_lines, 'transcript_id "([^"]+)"')[, 2]
gene_name <- str_match(tx_lines, 'gene_name "([^"]+)"')[, 2]
tx2gene <- data.frame(transcript_id = transcript_id, gene_id = gene_id) %>% distinct()
gene2name <- data.frame(gene_id = gene_id, gene_name = gene_name) %>%
  filter(!is.na(gene_name)) %>% distinct(gene_id, .keep_all = TRUE)

productivity <- read.table(productivity_path, sep = "\t") %>%
  transmute(transcript_id = V4,
            productivity = ifelse(grepl("PRO", V13), "coding", "non-coding"))

gene_list <- read.csv(paste(in_dir, "csv/dtu_genes_list.csv", sep = ""))$gene_id
length(gene_list)  # 228

files <- list.files(path = salmon_path, pattern = "\\.sf$", full.names = TRUE, recursive = TRUE)
names(files) <- str_extract(files, "([CKTO]+\\d+_\\d)(?=\\.sf$)")
wt_files <- files[grepl("^CT", names(files))]

## ---- isoform-level abundance ----
txi_tx <- tximport(wt_files, type = "salmon", txOut = TRUE, countsFromAbundance = "no",
                    ignoreTxVersion = TRUE, dropInfReps = TRUE, importer = function(x) read.delim(x))

# expressed-transcript filter (>5 reads in >=2 WT samples), same definition as
# tx_to_gene_rhythmicity.R -- this was missing in the original script
counts_tx <- txi_tx$counts
mode(counts_tx) <- "integer"
expressed_tx <- rownames(counts_tx)[rowSums(counts_tx > 5) >= 2]

tx_ids <- tx2gene %>% filter(gene_id %in% gene_list, transcript_id %in% expressed_tx) %>% pull(transcript_id)

exp_design <- data.frame(sample_id = colnames(txi_tx$counts),
                         batch = factor(str_extract(colnames(txi_tx$counts), "(?<=_)\\d")),
                         time = as.numeric(str_extract(colnames(txi_tx$counts), "\\d+"))) %>%
  mutate(inphase = cos(2 * pi * time / 24), outphase = sin(2 * pi * time / 24))
rownames(exp_design) <- exp_design$sample_id

dds_tx <- DESeqDataSetFromTximport(txi_tx, colData = exp_design, design = ~ batch + inphase + outphase)
dds_tx <- dds_tx[rownames(dds_tx) %in% tx_ids, ]
dds_tx <- DESeq(dds_tx, reduced = ~batch, test = "LRT", fitType = "local")

abundance_tx <- varianceStabilizingTransformation(dds_tx, fitType = "local")
abundance_tx_no_batch <- limma::removeBatchEffect(assay(abundance_tx), batch = colData(dds_tx)$batch,
  design = model.matrix(~inphase + outphase, data = colData(dds_tx)))

## ---- gene-level abundance ----
txi_gene <- tximport(wt_files, type = "salmon", tx2gene = tx2gene, txOut = FALSE, countsFromAbundance = "no",
                      ignoreTxVersion = TRUE, dropInfReps = TRUE, importer = function(x) read.delim(x))

dds_gene <- DESeqDataSetFromTximport(txi_gene, colData = exp_design, design = ~ batch + inphase + outphase)
dds_gene <- dds_gene[rownames(dds_gene) %in% gene_list, ]
dds_gene <- DESeq(dds_gene, reduced = ~batch, test = "LRT", fitType = "local")

abundance_gene <- varianceStabilizingTransformation(dds_gene, fitType = "local")
abundance_gene_no_batch <- limma::removeBatchEffect(assay(abundance_gene), batch = colData(dds_gene)$batch,
  design = model.matrix(~inphase + outphase, data = colData(dds_gene)))

## ---- long-format data, both levels ----
df_tx <- abundance_tx_no_batch %>%
  as.data.frame() %>% rownames_to_column("transcript_id") %>%
  pivot_longer(-transcript_id, names_to = "sample_id", values_to = "value") %>%
  left_join(exp_design %>% dplyr::select(sample_id, time), by = "sample_id") %>%
  left_join(tx2gene, by = "transcript_id") %>%
  left_join(gene2name, by = "gene_id") %>%
  left_join(productivity, by = "transcript_id") %>%
  mutate(productivity = ifelse(is.na(productivity), "non-coding", productivity),
         series_id = transcript_id)

df_gene <- abundance_gene_no_batch %>%
  as.data.frame() %>% rownames_to_column("gene_id") %>%
  pivot_longer(-gene_id, names_to = "sample_id", values_to = "value") %>%
  left_join(exp_design %>% dplyr::select(sample_id, time), by = "sample_id") %>%
  left_join(gene2name, by = "gene_id") %>%
  mutate(series_id = paste0(gene_id, "_geneagg"))

genes_ordered <- gene2name %>% filter(gene_id %in% gene_list) %>% arrange(gene_name) %>% pull(gene_name)
length(genes_ordered)

## ---- overlaid isoform + gene time course, ~12 genes per page ----
pdf(paste(out_dir, "figures/fig1c_overlap_timecourse_dtugenes.pdf", sep = ""), width = 10, height = 7.5)
chunks <- split(genes_ordered, ceiling(seq_along(genes_ordered) / genes_per_page))
for (chunk in chunks) {
  df_tx_chunk <- filter(df_tx, gene_name %in% chunk)
  df_gene_chunk <- filter(df_gene, gene_name %in% chunk)

  p <- ggplot() +
    geom_point(data = df_tx_chunk, aes(x = time, y = value, color = productivity), size = 0.4) +
    geom_smooth(data = df_tx_chunk, aes(x = time, y = value, group = series_id, color = productivity),
                linewidth = 0.6, formula = "y ~ x", method = "loess", se = FALSE) +
    scale_color_colorblind() +
    new_scale_color() +
    geom_point(data = df_gene_chunk, aes(x = time, y = value, color = "gene (aggregate)"), size = 0.9) +
    geom_smooth(data = df_gene_chunk, aes(x = time, y = value, group = series_id, color = "gene (aggregate)"),
                linewidth = 1, formula = "y ~ x", method = "loess", se = FALSE) +
    scale_color_manual(values = c("gene (aggregate)" = "blue")) +
    facet_wrap(~gene_name, scales = "free_y", ncol = 4, nrow = 3) +
    scale_x_continuous(breaks = seq(0, 20, 4)) +
    xlab("circadian time") + ylab(bquote("log"[2] ~ "expression")) +
    theme_classic(base_size = 10) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          strip.text = element_text(size = 9), strip.background = element_rect(fill = "white", color = "white"),
          legend.title = element_blank(), legend.position = "top", aspect.ratio = 0.85)
  print(p)
}
dev.off()
