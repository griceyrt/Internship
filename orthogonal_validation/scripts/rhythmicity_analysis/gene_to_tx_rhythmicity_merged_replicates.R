#!/usr/bin/env Rscript
# Reviewer 3, point 1 -- gene-level circadian rhythmicity, WT 6-timepoint series.
# Merged-replicate variant of ../gene_to_tx_rhythmicity.R -- Bharath (2026-08-21
# email) didn't like panel (c) splitting the 2 replicates into two separate
# heatmap panels, since it visually implies two different/separate runs when
# there's nothing special about a sample being "replicate 1" vs "replicate 2".
# Here all 12 samples (6 timepoints x 2 reps) are columns in ONE heatmap,
# ordered by time then replicate, instead of facet_wrap(~rep). Original script
# left untouched for reproducibility -- this writes to the same figures/ dir
# with a "_merged" suffix so nothing gets overwritten.
# Author: Gricey
# Run from orthogonal_validation/: Rscript scripts/rhythmicity_analysis/gene_to_tx_rhythmicity_merged_replicates.R

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

results_valid <- results_df %>% filter(!is.na(padj))
sig <- results_valid %>% filter(padj < fdr_cutoff)

nrow(sig)
100 * nrow(sig) / nrow(results_valid)

##---------------Heat map of cycling genes (Fig 1B replacement), MERGED-------##
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

# 12 columns (6 timepoints x 2 reps), ordered time-then-rep -- no facet, so
# there's no visual break implying two separate runs
sample_levels <- paste(rep(seq(0, 20, 4), each = 2), rep(c(1, 2), times = 6), sep = "_")

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
         sample_col = factor(paste(time, rep, sep = "_"), levels = sample_levels, ordered = TRUE)) %>%
  arrange(gene_id)

# light/dark blocks: columns 1-6 = CT0/CT4/CT8 pairs (subjective day), 7-12 =
# CT12/CT16/CT20 pairs (subjective night) -- one continuous block each, no
# per-replicate duplication needed now that there's a single panel
LD_bars <- data.frame(xmin = c(0.5, 6.5), xmax = c(6.5, 12.5), col = c("grey50", "black"))

heatmap_gene <- ggplot(df_1) +
  geom_tile(aes(y = gene_id, x = sample_col, fill = gene_exp), width = 1) +
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
  scale_x_discrete(expand = expansion(), labels = rep(seq(0, 20, 4), each = 2)) +
  scale_y_discrete(expand = expansion()) +
  theme(text = element_text(family = "Arial"),
        axis.text.y = element_blank(), axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top", panel.grid = element_blank(), panel.grid.major = element_blank(),
        axis.text.x = element_text(family = "Arial", size = 9),
        axis.title.x = element_text(family = "Arial", size = 12),
        axis.title.y = element_text(family = "Arial", size = 12))

ggsave(paste(out_dir, "figures/gene_rhythmicity_heatmap_merged.svg", sep=""), plot = heatmap_gene, width = 3, height = 4.5)

##---------------Heat map #2: genes + isoforms, MERGED-------------------------##
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
         sample_col = factor(paste(time, rep, sep = "_"), levels = sample_levels, ordered = TRUE)) %>%
  arrange(transcript_id)

heatmap_isoforms <- ggplot(df_2) +
  geom_tile(aes(y = transcript_id, x = sample_col, fill = tx_exp), width = 1) +
  scale_fill_gradient2_tableau(limits = c(-2.5, 2.5), palette = "Green-Blue Diverging", guide = "none", oob = scales::squish) +
  new_scale_fill() +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = n_tx_rhy, ymax = n_tx_rhy * 1.05, fill = col),
            data = LD_bars, inherit.aes = FALSE) +
  scale_fill_manual(values = c("black", "grey80"), guide = "none") +
  ylab(paste(n_tx_rhy, "isoforms")) +
  xlab("Time (h)") + theme_minimal(base_size = 11, base_family = "Arial") +
  scale_x_discrete(expand = expansion(), labels = rep(seq(0, 20, 4), each = 2)) +
  scale_y_discrete(expand = expansion()) +
  theme(text = element_text(family = "Arial"),
        axis.text.y = element_blank(), axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none", panel.grid = element_blank(), panel.grid.major = element_blank(),
        axis.text.x = element_text(family = "Arial", size = 9),
        axis.title.x = element_text(family = "Arial", size = 12),
        axis.title.y = element_text(family = "Arial", size = 12))

# shrink the gene panel's aspect ratio by n_small/n_large so both panels keep the same row height
heatmap_gene_scaled <- heatmap_gene + theme(aspect.ratio = 4 * n_genes_rhy / n_tx_rhy)

heatmap_genes_isoforms <- heatmap_gene_scaled | heatmap_isoforms
ggsave(paste(out_dir, "figures/gene_isoform_heatmap_merged.svg", sep=""), plot = heatmap_genes_isoforms, width = 7, height = 5)

# composite-ready version: x-axis stripped from the top block (the masked-genes
# block stacked underneath keeps its own x-axis, so one label set serves both)
heatmap_genes_isoforms_composite <- (heatmap_gene_scaled + theme(axis.text.x = element_blank(), axis.title.x = element_blank())) |
  (heatmap_isoforms + theme(axis.text.x = element_blank(), axis.title.x = element_blank()))
ggsave(paste(out_dir, "figures/gene_isoform_heatmap_composite_merged.svg", sep=""), plot = heatmap_genes_isoforms_composite, width = 7, height = 5)

# minimal version for pptx -- y-axis titles ("N rhythmic genes"/"N isoforms")
# and x-axis title ("Time (h)") removed; Kiran adds those back as textboxes
heatmap_gene_minimal <- heatmap_gene_scaled +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
heatmap_isoforms_minimal <- heatmap_isoforms +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
heatmap_genes_isoforms_minimal <- heatmap_gene_minimal | heatmap_isoforms_minimal
ggsave(paste(out_dir, "figures/gene_isoform_heatmap_merged_minimal.svg", sep=""), plot = heatmap_genes_isoforms_minimal, width = 7, height = 5)
