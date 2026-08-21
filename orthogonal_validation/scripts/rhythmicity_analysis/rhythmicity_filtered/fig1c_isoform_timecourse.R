#!/usr/bin/env Rscript
# Fig 1c redo -- individual isoform expression time courses, coding (black) vs
# non-coding (yellow), for genes with >1 isoform among the 773 gene-level
# significant list, plus a few named clock genes (Cry2, Arntl2) on request.
# Filtered version of ../fig1c_isoform_timecourse.R -- that script pulled
# every annotated isoform of a gene from the GTF with no expression check, so
# unexpressed isoforms got plotted alongside real ones (Khushi, 2026-08-20
# email). Here tx_ids is restricted to expressed transcripts (>5 reads in >=2
# WT samples), same definition as tx_to_gene_rhythmicity.R.
# Author: Gricey
# Run from orthogonal_validation/: Rscript scripts/rhythmicity_analysis/rhythmicity_filtered/fig1c_isoform_timecourse.R

library(tximport)
library(DESeq2)
library(stringr)
library(dplyr)
library(tidyr)
library(tibble)
library(limma)
library(ggplot2)
library(ggthemes)

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

## ---- gene list: 773 significant AND >1 isoform, plus named clock genes ----
sig773 <- read.csv(paste(in_dir, "csv/gene_rhythmicity_significant.csv", sep = ""))$gene_id
n_tx_per_gene <- tx2gene %>% count(gene_id, name = "n_tx")
multi_iso_genes <- n_tx_per_gene %>% filter(n_tx > 1) %>% pull(gene_id)

manuscript_clock_genes <- c("Arntl", "Per1", "Dbp", "Rorc", "Clock", "Hlf", "Bhlhe40",
                             "Nr1d2", "Cry1", "Cry2", "Nr1d1", "Per2", "Arntl2",
                             "Nampt", "Ciart", "Insig2", "Fasn")
manuscript_clock_ids <- gene2name %>% filter(gene_name %in% manuscript_clock_genes) %>%
  filter(gene_id %in% multi_iso_genes) %>% pull(gene_id)

gene_list <- union(intersect(sig773, multi_iso_genes), manuscript_clock_ids)
length(gene_list)  # candidate genes to plot

## ---- VST + batch-corrected abundance, same design as the rest of this project ----
files <- list.files(path = salmon_path, pattern = "\\.sf$", full.names = TRUE, recursive = TRUE)
names(files) <- str_extract(files, "([CKTO]+\\d+_\\d)(?=\\.sf$)")
wt_files <- files[grepl("^CT", names(files))]

txi <- tximport(wt_files, type = "salmon", txOut = TRUE, countsFromAbundance = "no",
                 ignoreTxVersion = TRUE, dropInfReps = TRUE, importer = function(x) read.delim(x))

# expressed-transcript filter (>5 reads in >=2 WT samples), same definition as
# tx_to_gene_rhythmicity.R -- this was missing in the original script
counts_tx <- txi$counts
mode(counts_tx) <- "integer"
expressed_tx <- rownames(counts_tx)[rowSums(counts_tx > 5) >= 2]

tx_ids <- tx2gene %>% filter(gene_id %in% gene_list, transcript_id %in% expressed_tx) %>% pull(transcript_id)

exp_design <- data.frame(sample_id = colnames(txi$counts),
                         batch = factor(str_extract(colnames(txi$counts), "(?<=_)\\d")),
                         time = as.numeric(str_extract(colnames(txi$counts), "\\d+"))) %>%
  mutate(inphase = cos(2 * pi * time / 24), outphase = sin(2 * pi * time / 24))
rownames(exp_design) <- exp_design$sample_id

dds <- DESeqDataSetFromTximport(txi, colData = exp_design, design = ~ batch + inphase + outphase)
dds <- dds[rownames(dds) %in% tx_ids, ]
dds <- DESeq(dds, reduced = ~batch, test = "LRT", fitType = "local")

abundance <- varianceStabilizingTransformation(dds, fitType = "local")
abundance_no_batch <- limma::removeBatchEffect(assay(abundance), batch = colData(dds)$batch,
  design = model.matrix(~inphase + outphase, data = colData(dds)))

## ---- long-format plotting data, matching Figure2's df_3 ----
df <- abundance_no_batch %>%
  as.data.frame() %>% rownames_to_column("transcript_id") %>%
  pivot_longer(-transcript_id, names_to = c("sample", "rep"), names_sep = "_", values_to = "value") %>%
  separate(sample, c("cond", "time"), convert = TRUE, sep = 2) %>%
  dplyr::select(-cond) %>%
  left_join(tx2gene, by = "transcript_id") %>%
  left_join(gene2name, by = "gene_id") %>%
  left_join(productivity, by = "transcript_id") %>%
  mutate(productivity = ifelse(is.na(productivity), "non-coding", productivity))

genes_ordered <- gene2name %>% filter(gene_id %in% unique(df$gene_id)) %>% arrange(gene_name) %>% pull(gene_name)
length(genes_ordered)

## ---- paginated time-course plot, ~12 genes per page ----
pdf(paste(out_dir, "figures/fig1c_isoform_timecourse.pdf", sep = ""), width = 10, height = 7.5)
chunks <- split(genes_ordered, ceiling(seq_along(genes_ordered) / genes_per_page))
for (chunk in chunks) {
  p <- ggplot(data = filter(df, gene_name %in% chunk)) +
    geom_point(aes(x = time, y = value, color = productivity), size = 0.4) +
    geom_smooth(aes(x = time, y = value, group = transcript_id, color = productivity),
                linewidth = 0.7, formula = "y ~ x", method = "loess", se = FALSE) +
    facet_wrap(~gene_name, scales = "free_y", ncol = 4, nrow = 3) +
    scale_x_continuous(breaks = seq(0, 20, 4)) +
    xlab("circadian time") + ylab(bquote("log"[2] ~ "isoform expression")) +
    theme_classic(base_size = 10) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          strip.text = element_text(size = 9), strip.background = element_rect(fill = "white", color = "white"),
          legend.title = element_blank(), legend.position = "top", aspect.ratio = 0.85) +
    scale_color_colorblind()
  print(p)
}
dev.off()
