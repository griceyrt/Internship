#!/usr/bin/env Rscript
# Reviewer 3, point 1 -- gene-level circadian rhythmicity, WT 6-timepoint series.
# Adapted from Figure2_2024_12_03.R's cosinor design, gene level, no stageR.
# Author: Gricey
# Run from orthogonal_validation/: Rscript scripts/gene_to_tx_rhythmicity.R

library(tximport)
library(DESeq2)
library(stringr)
library(dplyr)
library(tibble)

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
gene_name <- str_match(tx_lines, 'gene_name "([^"]+)"')[, 2]
tx2gene <- data.frame(transcript_id = transcript_id, gene_id = gene_id) %>% distinct()
gene2name <- data.frame(gene_id = gene_id, gene_name = gene_name) %>%
  filter(!is.na(gene_name)) %>% distinct(gene_id, .keep_all = TRUE)

files <- list.files(path = salmon_path, pattern = "\\.sf$", full.names = TRUE, recursive = TRUE)
names(files) <- str_extract(files, "([CKTO]+\\d+_\\d)(?=\\.sf$)")

txi_gene_all <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut = FALSE,
                          countsFromAbundance = "no", ignoreTxVersion = TRUE, dropInfReps = TRUE,
                          importer = function(x) read.delim(x))
counts_gene_all <- txi_gene_all$counts
mode(counts_gene_all) <- "integer"
expressed_genes <- rownames(counts_gene_all)[rowSums(counts_gene_all > count_threshold) >= min_samples]
write.csv(data.frame(gene_id = expressed_genes), paste(out_dir, "csv/expressed_gene_list.csv", sep=""), row.names = FALSE)
length(expressed_genes)

wt_files <- files[grepl("^CT", names(files))]

txi_gene_wt <- tximport(wt_files, type = "salmon", tx2gene = tx2gene, txOut = FALSE,
                         countsFromAbundance = "no", ignoreTxVersion = TRUE, dropInfReps = TRUE,
                         importer = function(x) read.delim(x))

exp_design <- data.frame(sample_id = colnames(txi_gene_wt$counts),
                         batch = factor(str_extract(colnames(txi_gene_wt$counts), "(?<=_)\\d")),
                         time = as.numeric(str_extract(colnames(txi_gene_wt$counts), "\\d+"))) %>%
  mutate(inphase = cos(2*pi*time/24), outphase = sin(2*pi*time/24))
rownames(exp_design) <- exp_design$sample_id

dds <- DESeqDataSetFromTximport(txi_gene_wt, colData = exp_design, design = ~ batch + inphase + outphase)
dds <- dds[rownames(dds) %in% expressed_genes, ]
dds <- DESeq(dds, reduced = ~batch, test = "LRT", fitType = "local")

res <- DESeq2::results(dds, independentFiltering = TRUE, alpha = fdr_cutoff, tidy = TRUE)

results_df <- res %>%
  dplyr::rename(gene_id = row) %>%
  left_join(gene2name, by = "gene_id") %>%
  cbind(coef(dds)[, c("inphase", "outphase")]) %>%
  mutate(phase = 12/pi * (atan2(outphase, inphase) %% (2*pi)),
         amp = sqrt(inphase^2 + outphase^2)) %>%
  relocate(gene_id, gene_name)
write.csv(results_df, paste(out_dir, "csv/gene_rhythmicity_results.csv", sep=""), row.names = FALSE)

results_valid <- results_df %>% filter(!is.na(padj))
sig <- results_valid %>% filter(padj < fdr_cutoff)
write.csv(sig, paste(out_dir, "csv/gene_rhythmicity_significant.csv", sep=""), row.names = FALSE)

nrow(sig)
100 * nrow(sig) / nrow(results_valid)

##---------------Heat map of cycling genes (Fig 1B replacement)---------------##
library(ggplot2)
library(ggthemes)
library(ggnewscale)
library(tidyr)
library(limma)
library(grid)

abundance <- varianceStabilizingTransformation(dds, fitType = "local")
abundance_no_batch_effect <- limma::removeBatchEffect(assay(abundance), batch = colData(dds)$batch,
  design = model.matrix(~inphase + outphase, data = colData(dds)))

phase_sorted_genes <- sig %>% arrange(desc(phase)) %>% pull(gene_id)
n_genes_rhy <- length(phase_sorted_genes)

df_1 <- abundance_no_batch_effect %>%
  as.data.frame() %>% rownames_to_column("gene_id") %>%
  filter(gene_id %in% phase_sorted_genes) %>%
  pivot_longer(-gene_id, names_to = c("sample", "rep"), names_sep = "_", values_to = "gene_exp") %>%
  separate(sample, c("cond", "time"), convert = TRUE, sep = 2) %>%
  dplyr::select(-cond) %>%
  group_by(gene_id) %>%
  mutate(gene_exp = scale(gene_exp)) %>%
  ungroup() %>%
  mutate(gene_id = factor(gene_id, levels = phase_sorted_genes, ordered = TRUE),
         time = factor(time, levels = seq(0, 20, 4), ordered = TRUE)) %>%
  arrange(gene_id)

LD_bars <- data.frame(xmin = rep(c(0.5, 3.5), 2), xmax = rep(c(3.5, 6.5), 2),
                       col = c("grey50", "black"), rep = c(1, 1, 2, 2))

heatmap_gene <- ggplot(df_1) +
  geom_tile(aes(y = gene_id, x = time, fill = gene_exp), width = 1) +
  facet_wrap(~rep, nrow = 1, labeller = as_labeller(function(x) paste("replicate", x))) +
  scale_fill_gradient2_tableau(name = NULL, limits = c(-2.5, 2.5), palette = "Green-Blue Diverging",
                                oob = scales::squish,
                                guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                                        barwidth = unit(3, "cm"), barheight = unit(0.2, "cm"), reverse = TRUE)) +
  new_scale_fill() +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = n_genes_rhy, ymax = n_genes_rhy * 1.05, fill = col),
            data = LD_bars, inherit.aes = FALSE) +
  scale_fill_manual(values = c("black", "grey80"), guide = "none") +
  ylab(paste(n_genes_rhy, "rhythmic genes")) +
  xlab("Time (h)") + theme_minimal(base_size = 11, base_family = "Arial") +
  scale_x_discrete(expand = expansion()) + scale_y_discrete(expand = expansion()) +
  theme(text = element_text(family = "Arial"),
        axis.text.y = element_blank(), axis.line = element_blank(),
        axis.ticks.y = element_blank(), panel.spacing.x = unit(0.2, "lines"),
        legend.position = "top", panel.grid = element_blank(), panel.grid.major = element_blank(),
        strip.text = element_text(family = "Arial", size = 11),
        axis.text.x = element_text(family = "Arial", size = 9),
        axis.title.x = element_text(family = "Arial", size = 12),
        axis.title.y = element_text(family = "Arial", size = 12))

ggsave(paste(out_dir, "figures/gene_rhythmicity_heatmap.svg", sep=""), plot = heatmap_gene, width = 3, height = 4.5)

##---------------Heat map #2: genes + isoforms, Lies-style scale--------------##
library(patchwork)

tx_ids <- tx2gene %>% filter(gene_id %in% phase_sorted_genes) %>% pull(transcript_id)

txi_tx_wt <- tximport(wt_files, type = "salmon", txOut = TRUE,
                       countsFromAbundance = "no", ignoreTxVersion = TRUE, dropInfReps = TRUE,
                       importer = function(x) read.delim(x))
counts_tx_wt <- txi_tx_wt$counts
mode(counts_tx_wt) <- "integer"
keep_tx <- rownames(counts_tx_wt) %in% tx_ids & rowSums(counts_tx_wt > count_threshold) >= min_samples

dds_tx <- DESeqDataSetFromTximport(txi_tx_wt, colData = exp_design, design = ~ batch + inphase + outphase)
dds_tx <- dds_tx[keep_tx, ]
dds_tx <- DESeq(dds_tx, reduced = ~batch, test = "LRT", fitType = "local")

abundance_tx <- varianceStabilizingTransformation(dds_tx, fitType = "local")
abundance_tx_no_batch <- limma::removeBatchEffect(assay(abundance_tx), batch = colData(dds_tx)$batch,
  design = model.matrix(~inphase + outphase, data = colData(dds_tx)))

tx_gene_order <- tx2gene %>% filter(gene_id %in% phase_sorted_genes, transcript_id %in% rownames(dds_tx)) %>%
  mutate(gene_id = factor(gene_id, levels = phase_sorted_genes, ordered = TRUE)) %>%
  arrange(gene_id, transcript_id) %>%
  pull(transcript_id)
n_tx_rhy <- length(tx_gene_order)

df_2 <- abundance_tx_no_batch %>%
  as.data.frame() %>% rownames_to_column("transcript_id") %>%
  filter(transcript_id %in% tx_gene_order) %>%
  pivot_longer(-transcript_id, names_to = c("sample", "rep"), names_sep = "_", values_to = "tx_exp") %>%
  separate(sample, c("cond", "time"), convert = TRUE, sep = 2) %>%
  dplyr::select(-cond) %>%
  group_by(transcript_id) %>%
  mutate(tx_exp = scale(tx_exp)) %>%
  ungroup() %>%
  mutate(transcript_id = factor(transcript_id, levels = tx_gene_order, ordered = TRUE),
         time = factor(time, levels = seq(0, 20, 4), ordered = TRUE)) %>%
  arrange(transcript_id)

LD_bars_tx <- data.frame(xmin = rep(c(0.5, 3.5), 2), xmax = rep(c(3.5, 6.5), 2),
                          col = c("grey50", "black"), rep = c(1, 1, 2, 2))

heatmap_isoforms <- ggplot(df_2) +
  geom_tile(aes(y = transcript_id, x = time, fill = tx_exp), width = 1) +
  facet_wrap(~rep, nrow = 1, labeller = as_labeller(function(x) paste("replicate", x))) +
  scale_fill_gradient2_tableau(limits = c(-2.5, 2.5), palette = "Green-Blue Diverging", guide = "none", oob = scales::squish) +
  new_scale_fill() +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = n_tx_rhy, ymax = n_tx_rhy * 1.05, fill = col),
            data = LD_bars_tx, inherit.aes = FALSE) +
  scale_fill_manual(values = c("black", "grey80"), guide = "none") +
  ylab(paste(n_tx_rhy, "isoforms")) +
  xlab("Time (h)") + theme_minimal(base_size = 11, base_family = "Arial") +
  scale_x_discrete(expand = expansion()) + scale_y_discrete(expand = expansion()) +
  theme(text = element_text(family = "Arial"),
        axis.text.y = element_blank(), axis.line = element_blank(),
        axis.ticks.y = element_blank(), panel.spacing.x = unit(0.2, "lines"),
        legend.position = "none", panel.grid = element_blank(), panel.grid.major = element_blank(),
        strip.text = element_text(family = "Arial", size = 11),
        axis.text.x = element_text(family = "Arial", size = 9),
        axis.title.x = element_text(family = "Arial", size = 12),
        axis.title.y = element_text(family = "Arial", size = 12))

# shrink the gene panel's aspect ratio by n_small/n_large so both panels keep the same row height
heatmap_gene_scaled <- heatmap_gene + theme(aspect.ratio = 4 * n_genes_rhy / n_tx_rhy)

heatmap_genes_isoforms <- heatmap_gene_scaled | heatmap_isoforms
ggsave(paste(out_dir, "figures/gene_isoform_heatmap.svg", sep=""), plot = heatmap_genes_isoforms, width = 7, height = 5)

# composite-ready version: x-axis stripped from the top block (the masked-genes
# block stacked underneath keeps its own x-axis, so one label set serves both)
heatmap_genes_isoforms_composite <- (heatmap_gene_scaled + theme(axis.text.x = element_blank(), axis.title.x = element_blank())) |
  (heatmap_isoforms + theme(axis.text.x = element_blank(), axis.title.x = element_blank()))
ggsave(paste(out_dir, "figures/gene_isoform_heatmap_composite.svg", sep=""), plot = heatmap_genes_isoforms_composite, width = 7, height = 5)
