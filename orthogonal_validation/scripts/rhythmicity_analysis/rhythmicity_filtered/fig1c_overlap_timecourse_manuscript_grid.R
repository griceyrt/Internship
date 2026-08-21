#!/usr/bin/env Rscript
# Fig 1C manuscript grid (page 1, 16 genes) + supplementary grid (page 2, 16
# genes), OVERLAP style -- isoform-level (coding=black, non-coding=orange)
# and gene-level (blue) curves on the SAME panel per gene. Both grids are
# fixed 4-column layouts with genes placed in the exact requested order
# (filled top-to-bottom within each column, then left-to-right across columns).
# Filtered version of ../fig1c_overlap_timecourse_manuscript_grid.R --
# restricts the isoform side (tx_ids) to expressed transcripts (>5 reads in
# >=2 WT samples), same definition as tx_to_gene_rhythmicity.R. Khushi
# noticed (2026-08-20 email) that genes like Arntl/Insig2/Nampt/Noct showed
# more low-expression isoforms here than in her reference plot; the original
# script pulled every annotated isoform of a gene from the GTF with no
# expression check at all.
# Author: Gricey
# Run from orthogonal_validation/: Rscript scripts/rhythmicity_analysis/rhythmicity_filtered/fig1c_overlap_timecourse_manuscript_grid.R

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
out_dir <- "results/rhythmicity_analysis/rhythmicity_filtered/"

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

# column-major order: top-to-bottom within each column, then left-to-right across columns
page1_genes <- c("Arntl", "Per1", "Dbp", "Rorc",
                  "Clock", "Hlf", "Nampt", "Ciart",
                  "Insig2", "Wee1", "Tef", "Noct",
                  "Rpp21", "Spry1", "Rif1", "Fchsd2")

page2_genes <- c("Abhd16a", "BC002059", "Cdc42bpa", "Cyp2d13",
                  "Styx", "Elp6", "Gga1", "Mlh1",
                  "Ddo", "Mrpl15", "Ercc8", "Fmr1",
                  "Mul1", "Macroh2a1", "Arpc3", "Micos13")

all_requested <- unique(c(page1_genes, page2_genes))

# case-insensitive match against GTF gene_name -- a couple of the requested
# symbols were typed in non-canonical case (e.g. macroh2a1)
name_lookup <- gene2name %>% mutate(name_upper = toupper(gene_name))
resolved <- name_lookup %>% filter(name_upper %in% toupper(all_requested))
missing <- setdiff(toupper(all_requested), resolved$name_upper)
if (length(missing) > 0) {
  cat("WARNING -- not found in GTF gene_name (check for alias/renaming):\n")
  print(missing)
}

resolved_map <- setNames(resolved$gene_name, resolved$name_upper)
gene_list <- unname(resolved_map[toupper(all_requested)])
gene_list <- gene_list[!is.na(gene_list)]

gene_ids <- gene2name %>% filter(gene_name %in% gene_list) %>% pull(gene_id)

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

tx_ids <- tx2gene %>% filter(gene_id %in% gene_ids, transcript_id %in% expressed_tx) %>% pull(transcript_id)

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
dds_gene <- dds_gene[rownames(dds_gene) %in% gene_ids, ]
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

## ---- one grid page, genes placed column-major via facet_wrap(dir = "v") ----
make_grid_page <- function(gene_order) {
  df_tx_chunk <- df_tx %>% filter(gene_name %in% gene_order) %>%
    mutate(gene_name = factor(gene_name, levels = gene_order, ordered = TRUE))
  df_gene_chunk <- df_gene %>% filter(gene_name %in% gene_order) %>%
    mutate(gene_name = factor(gene_name, levels = gene_order, ordered = TRUE))

  ggplot() +
    geom_point(data = df_tx_chunk, aes(x = time, y = value, color = productivity), size = 0.4) +
    geom_smooth(data = df_tx_chunk, aes(x = time, y = value, group = series_id, color = productivity),
                linewidth = 0.6, formula = "y ~ x", method = "loess", se = FALSE) +
    scale_color_colorblind() +
    new_scale_color() +
    geom_point(data = df_gene_chunk, aes(x = time, y = value, color = "gene (aggregate)"), size = 0.9) +
    geom_smooth(data = df_gene_chunk, aes(x = time, y = value, group = series_id, color = "gene (aggregate)"),
                linewidth = 1, formula = "y ~ x", method = "loess", se = FALSE) +
    scale_color_manual(values = c("gene (aggregate)" = "blue")) +
    facet_wrap(~gene_name, scales = "free_y", ncol = 4, nrow = 4, dir = "v") +
    scale_x_continuous(breaks = seq(0, 20, 4)) +
    xlab("circadian time") + ylab(bquote("log"[2] ~ "expression")) +
    theme_classic(base_size = 10) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          strip.text = element_text(size = 9), strip.background = element_rect(fill = "white", color = "white"),
          legend.title = element_blank(), legend.position = "top", aspect.ratio = 0.85)
}

pdf(paste(out_dir, "figures/fig1c_overlap_timecourse_manuscript_grid.pdf", sep = ""), width = 10, height = 10)
print(make_grid_page(page1_genes))
print(make_grid_page(page2_genes))
dev.off()

## ---- minimal / manuscript-styled version for pptx: same curves/points
## (unchanged, so the tiny multiples still read at this size), thicker axis
## lines and ticks, bigger axis numbers, gene labels in italic Arial, bigger
## legend, no axis titles (Kiran adds "log2 expression"/"circadian time"
## manually) ----
make_grid_page_minimal <- function(gene_order) {
  df_tx_chunk <- df_tx %>% filter(gene_name %in% gene_order) %>%
    mutate(gene_name = factor(gene_name, levels = gene_order, ordered = TRUE))
  df_gene_chunk <- df_gene %>% filter(gene_name %in% gene_order) %>%
    mutate(gene_name = factor(gene_name, levels = gene_order, ordered = TRUE))

  ggplot() +
    geom_point(data = df_tx_chunk, aes(x = time, y = value, color = productivity), size = 0.4) +
    geom_smooth(data = df_tx_chunk, aes(x = time, y = value, group = series_id, color = productivity),
                linewidth = 0.6, formula = "y ~ x", method = "loess", se = FALSE) +
    scale_color_colorblind() +
    new_scale_color() +
    geom_point(data = df_gene_chunk, aes(x = time, y = value, color = "gene (aggregate)"), size = 0.9) +
    geom_smooth(data = df_gene_chunk, aes(x = time, y = value, group = series_id, color = "gene (aggregate)"),
                linewidth = 1, formula = "y ~ x", method = "loess", se = FALSE) +
    scale_color_manual(values = c("gene (aggregate)" = "blue")) +
    facet_wrap(~gene_name, scales = "free_y", ncol = 4, nrow = 4, dir = "v") +
    scale_x_continuous(breaks = seq(0, 20, 4)) +
    theme_classic(base_size = 10) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(linewidth = 0.9), axis.ticks = element_line(linewidth = 0.9),
          axis.text = element_text(size = 10), axis.title = element_blank(),
          strip.text = element_text(family = "Arial", face = "italic", size = 12),
          strip.background = element_rect(fill = "white", color = "white"),
          legend.title = element_blank(), legend.text = element_text(size = 12),
          legend.key.size = unit(0.5, "cm"), legend.position = "top", aspect.ratio = 0.85)
}

# base pdf() device only knows the 14 built-in Type1 fonts and errors on
# "Arial" -- cairo_pdf() can embed real system fonts instead
cairo_pdf(paste(out_dir, "figures/fig1c_overlap_timecourse_manuscript_grid_minimal.pdf", sep = ""), width = 10, height = 10)
print(make_grid_page_minimal(page1_genes))
print(make_grid_page_minimal(page2_genes))
dev.off()
