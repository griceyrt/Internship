#!/usr/bin/env Rscript
# Reviewer 3, point 1 -- "backwards" heatmap: rhythmic transcripts, zoomed out
# to their parent genes (no gene-level filter). Companion to gene_to_tx_rhythmicity.R.
# Two row-ordering versions: independent (genes by own phase) and grouped
# (transcripts grouped by parent gene).
# Author: Gricey
# Run from orthogonal_validation/: Rscript scripts/tx_to_gene_rhythmicity.R

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
wt_files <- files[grepl("^CT", names(files))]

exp_design <- data.frame(sample_id = names(wt_files),
                         batch = factor(str_extract(names(wt_files), "(?<=_)\\d")),
                         time = as.numeric(str_extract(names(wt_files), "\\d+"))) %>%
  mutate(inphase = cos(2*pi*time/24), outphase = sin(2*pi*time/24))
rownames(exp_design) <- exp_design$sample_id

##---------------Transcript-level rhythmicity, whole WT transcriptome---------##
txi_tx_wt <- tximport(wt_files, type = "salmon", txOut = TRUE,
                       countsFromAbundance = "no", ignoreTxVersion = TRUE, dropInfReps = TRUE,
                       importer = function(x) read.delim(x))
counts_tx_wt <- txi_tx_wt$counts
mode(counts_tx_wt) <- "integer"
expressed_tx <- rownames(counts_tx_wt)[rowSums(counts_tx_wt > count_threshold) >= min_samples]
length(expressed_tx)

dds_tx <- DESeqDataSetFromTximport(txi_tx_wt, colData = exp_design, design = ~ batch + inphase + outphase)
dds_tx <- dds_tx[rownames(dds_tx) %in% expressed_tx, ]
dds_tx <- DESeq(dds_tx, reduced = ~batch, test = "LRT", fitType = "local")

res_tx <- DESeq2::results(dds_tx, independentFiltering = TRUE, alpha = fdr_cutoff, tidy = TRUE)

tx_results <- res_tx %>%
  dplyr::rename(transcript_id = row) %>%
  left_join(tx2gene, by = "transcript_id") %>%
  left_join(gene2name, by = "gene_id") %>%
  cbind(coef(dds_tx)[, c("inphase", "outphase")]) %>%
  mutate(phase = 12/pi * (atan2(outphase, inphase) %% (2*pi)),
         amp = sqrt(inphase^2 + outphase^2)) %>%
  relocate(transcript_id, gene_id, gene_name)
write.csv(tx_results, paste(out_dir, "csv/tx_rhythmicity_results.csv", sep=""), row.names = FALSE)

tx_valid <- tx_results %>% filter(!is.na(padj))
sig_tx <- tx_valid %>% filter(padj < fdr_cutoff)
write.csv(sig_tx, paste(out_dir, "csv/tx_rhythmicity_significant.csv", sep=""), row.names = FALSE)

nrow(sig_tx)
100 * nrow(sig_tx) / nrow(tx_valid)
n_distinct(sig_tx$gene_id)

##---------------Gene-level abundance, zoomed-out parent genes----------------##
gene_results <- read.csv(paste(out_dir, "csv/gene_rhythmicity_results.csv", sep=""))
parent_genes <- unique(sig_tx$gene_id)

txi_gene_wt <- tximport(wt_files, type = "salmon", tx2gene = tx2gene, txOut = FALSE,
                         countsFromAbundance = "no", ignoreTxVersion = TRUE, dropInfReps = TRUE,
                         importer = function(x) read.delim(x))
dds_gene <- DESeqDataSetFromTximport(txi_gene_wt, colData = exp_design, design = ~ batch + inphase + outphase)
dds_gene <- dds_gene[rownames(dds_gene) %in% gene_results$gene_id[!is.na(gene_results$padj)], ]
dds_gene <- DESeq(dds_gene, reduced = ~batch, test = "LRT", fitType = "local")

abundance_gene <- varianceStabilizingTransformation(dds_gene, fitType = "local")
abundance_gene_no_batch <- limma::removeBatchEffect(assay(abundance_gene), batch = colData(dds_gene)$batch,
  design = model.matrix(~inphase + outphase, data = colData(dds_gene)))

abundance_tx <- varianceStabilizingTransformation(dds_tx, fitType = "local")
abundance_tx_no_batch <- limma::removeBatchEffect(assay(abundance_tx), batch = colData(dds_tx)$batch,
  design = model.matrix(~inphase + outphase, data = colData(dds_tx)))

##---------------Heat maps: two row-ordering versions--------------------------##
library(ggplot2)
library(ggthemes)
library(ggnewscale)
library(tidyr)
library(limma)
library(grid)
library(patchwork)

LD_bars <- data.frame(xmin = rep(c(0.5, 3.5), 2), xmax = rep(c(3.5, 6.5), 2),
                       col = c("grey50", "black"), rep = c(1, 1, 2, 2))

# ---- version 1: independent -- genes ordered by their own phase -------------
tx_order_indep <- sig_tx %>% arrange(desc(phase)) %>% pull(transcript_id)
n_tx_indep <- length(tx_order_indep)
gene_order_indep <- gene_results %>% filter(gene_id %in% parent_genes) %>% arrange(desc(phase)) %>% pull(gene_id)
n_gene_indep <- length(gene_order_indep)

df_tx_indep <- abundance_tx_no_batch %>%
  as.data.frame() %>% rownames_to_column("transcript_id") %>%
  filter(transcript_id %in% tx_order_indep) %>%
  pivot_longer(-transcript_id, names_to = c("sample", "rep"), names_sep = "_", values_to = "tx_exp") %>%
  separate(sample, c("cond", "time"), convert = TRUE, sep = 2) %>%
  dplyr::select(-cond) %>%
  group_by(transcript_id) %>%
  mutate(tx_exp = scale(tx_exp)) %>%
  ungroup() %>%
  mutate(transcript_id = factor(transcript_id, levels = tx_order_indep, ordered = TRUE),
         time = factor(time, levels = seq(0, 20, 4), ordered = TRUE)) %>%
  arrange(transcript_id)

df_gene_indep <- abundance_gene_no_batch %>%
  as.data.frame() %>% rownames_to_column("gene_id") %>%
  filter(gene_id %in% gene_order_indep) %>%
  pivot_longer(-gene_id, names_to = c("sample", "rep"), names_sep = "_", values_to = "gene_exp") %>%
  separate(sample, c("cond", "time"), convert = TRUE, sep = 2) %>%
  dplyr::select(-cond) %>%
  group_by(gene_id) %>%
  mutate(gene_exp = scale(gene_exp)) %>%
  ungroup() %>%
  mutate(gene_id = factor(gene_id, levels = gene_order_indep, ordered = TRUE),
         time = factor(time, levels = seq(0, 20, 4), ordered = TRUE)) %>%
  arrange(gene_id)

panel_tx_indep <- ggplot(df_tx_indep) +
  geom_tile(aes(y = transcript_id, x = time, fill = tx_exp), width = 1) +
  facet_wrap(~rep, nrow = 1, labeller = as_labeller(function(x) paste("replicate", x))) +
  scale_fill_gradient2_tableau(limits = c(-2.5, 2.5), palette = "Green-Blue Diverging", guide = "none", oob = scales::squish) +
  new_scale_fill() +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = n_tx_indep, ymax = n_tx_indep * 1.05, fill = col),
            data = LD_bars, inherit.aes = FALSE) +
  scale_fill_manual(values = c("black", "grey80"), guide = "none") +
  ylab(paste(n_tx_indep, "rhythmic transcripts")) +
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

panel_gene_indep <- ggplot(df_gene_indep) +
  geom_tile(aes(y = gene_id, x = time, fill = gene_exp), width = 1) +
  facet_wrap(~rep, nrow = 1, labeller = as_labeller(function(x) paste("replicate", x))) +
  scale_fill_gradient2_tableau(name = NULL, limits = c(-2.5, 2.5), palette = "Green-Blue Diverging",
                                oob = scales::squish,
                                guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                                        barwidth = unit(3, "cm"), barheight = unit(0.2, "cm"), reverse = TRUE)) +
  new_scale_fill() +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = n_gene_indep, ymax = n_gene_indep * 1.05, fill = col),
            data = LD_bars, inherit.aes = FALSE) +
  scale_fill_manual(values = c("black", "grey80"), guide = "none") +
  ylab(paste(n_gene_indep, "genes")) +
  xlab("Time (h)") + theme_minimal(base_size = 11, base_family = "Arial") +
  scale_x_discrete(expand = expansion()) + scale_y_discrete(expand = expansion()) +
  theme(text = element_text(family = "Arial"),
        axis.text.y = element_blank(), axis.line = element_blank(),
        axis.ticks.y = element_blank(), panel.spacing.x = unit(0.2, "lines"),
        legend.position = "top", panel.grid = element_blank(), panel.grid.major = element_blank(),
        strip.text = element_text(family = "Arial", size = 11),
        axis.text.x = element_text(family = "Arial", size = 9),
        axis.title.x = element_text(family = "Arial", size = 12),
        axis.title.y = element_text(family = "Arial", size = 12)) +
  theme(aspect.ratio = 4 * n_gene_indep / n_tx_indep)

heatmap_indep <- panel_tx_indep | panel_gene_indep
ggsave(paste(out_dir, "figures/tx_to_gene_heatmap_independent.svg", sep=""), plot = heatmap_indep, width = 7, height = 5)

# ---- version 2: grouped -- transcripts grouped by parent gene, gene order follows the group ----
gene_block_order <- gene_results %>% filter(gene_id %in% parent_genes) %>% arrange(desc(phase)) %>% pull(gene_id)
n_gene_grouped <- length(gene_block_order)
tx_order_grouped <- sig_tx %>%
  mutate(gene_id = factor(gene_id, levels = gene_block_order, ordered = TRUE)) %>%
  arrange(gene_id, transcript_id) %>%
  pull(transcript_id)
n_tx_grouped <- length(tx_order_grouped)

df_tx_grouped <- abundance_tx_no_batch %>%
  as.data.frame() %>% rownames_to_column("transcript_id") %>%
  filter(transcript_id %in% tx_order_grouped) %>%
  pivot_longer(-transcript_id, names_to = c("sample", "rep"), names_sep = "_", values_to = "tx_exp") %>%
  separate(sample, c("cond", "time"), convert = TRUE, sep = 2) %>%
  dplyr::select(-cond) %>%
  group_by(transcript_id) %>%
  mutate(tx_exp = scale(tx_exp)) %>%
  ungroup() %>%
  mutate(transcript_id = factor(transcript_id, levels = tx_order_grouped, ordered = TRUE),
         time = factor(time, levels = seq(0, 20, 4), ordered = TRUE)) %>%
  arrange(transcript_id)

df_gene_grouped <- abundance_gene_no_batch %>%
  as.data.frame() %>% rownames_to_column("gene_id") %>%
  filter(gene_id %in% gene_block_order) %>%
  pivot_longer(-gene_id, names_to = c("sample", "rep"), names_sep = "_", values_to = "gene_exp") %>%
  separate(sample, c("cond", "time"), convert = TRUE, sep = 2) %>%
  dplyr::select(-cond) %>%
  group_by(gene_id) %>%
  mutate(gene_exp = scale(gene_exp)) %>%
  ungroup() %>%
  mutate(gene_id = factor(gene_id, levels = gene_block_order, ordered = TRUE),
         time = factor(time, levels = seq(0, 20, 4), ordered = TRUE)) %>%
  arrange(gene_id)

panel_tx_grouped <- ggplot(df_tx_grouped) +
  geom_tile(aes(y = transcript_id, x = time, fill = tx_exp), width = 1) +
  facet_wrap(~rep, nrow = 1, labeller = as_labeller(function(x) paste("replicate", x))) +
  scale_fill_gradient2_tableau(limits = c(-2.5, 2.5), palette = "Green-Blue Diverging", guide = "none", oob = scales::squish) +
  new_scale_fill() +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = n_tx_grouped, ymax = n_tx_grouped * 1.05, fill = col),
            data = LD_bars, inherit.aes = FALSE) +
  scale_fill_manual(values = c("black", "grey80"), guide = "none") +
  ylab(paste(n_tx_grouped, "rhythmic transcripts")) +
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

panel_gene_grouped <- ggplot(df_gene_grouped) +
  geom_tile(aes(y = gene_id, x = time, fill = gene_exp), width = 1) +
  facet_wrap(~rep, nrow = 1, labeller = as_labeller(function(x) paste("replicate", x))) +
  scale_fill_gradient2_tableau(name = NULL, limits = c(-2.5, 2.5), palette = "Green-Blue Diverging",
                                oob = scales::squish,
                                guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                                        barwidth = unit(3, "cm"), barheight = unit(0.2, "cm"), reverse = TRUE)) +
  new_scale_fill() +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = n_gene_grouped, ymax = n_gene_grouped * 1.05, fill = col),
            data = LD_bars, inherit.aes = FALSE) +
  scale_fill_manual(values = c("black", "grey80"), guide = "none") +
  ylab(paste(n_gene_grouped, "genes")) +
  xlab("Time (h)") + theme_minimal(base_size = 11, base_family = "Arial") +
  scale_x_discrete(expand = expansion()) + scale_y_discrete(expand = expansion()) +
  theme(text = element_text(family = "Arial"),
        axis.text.y = element_blank(), axis.line = element_blank(),
        axis.ticks.y = element_blank(), panel.spacing.x = unit(0.2, "lines"),
        legend.position = "top", panel.grid = element_blank(), panel.grid.major = element_blank(),
        strip.text = element_text(family = "Arial", size = 11),
        axis.text.x = element_text(family = "Arial", size = 9),
        axis.title.x = element_text(family = "Arial", size = 12),
        axis.title.y = element_text(family = "Arial", size = 12)) +
  theme(aspect.ratio = 4 * n_gene_grouped / n_tx_grouped)

heatmap_grouped <- panel_tx_grouped | panel_gene_grouped
ggsave(paste(out_dir, "figures/tx_to_gene_heatmap_grouped.svg", sep=""), plot = heatmap_grouped, width = 7, height = 5)

# ---- version 3: masked genes -- gene-non-significant but >=1 significant isoform ----
gene_sig_ids <- gene_results$gene_id[!is.na(gene_results$padj) & gene_results$padj < fdr_cutoff]
masked_genes <- setdiff(parent_genes, gene_sig_ids)
length(masked_genes)

tx_masked <- sig_tx %>% filter(gene_id %in% masked_genes)
nrow(tx_masked)

gene_block_masked <- gene_results %>% filter(gene_id %in% masked_genes) %>% arrange(desc(phase)) %>% pull(gene_id)
n_gene_masked <- length(gene_block_masked)
tx_order_masked <- tx_masked %>%
  mutate(gene_id = factor(gene_id, levels = gene_block_masked, ordered = TRUE)) %>%
  arrange(gene_id, transcript_id) %>%
  pull(transcript_id)
n_tx_masked <- length(tx_order_masked)

df_tx_masked <- abundance_tx_no_batch %>%
  as.data.frame() %>% rownames_to_column("transcript_id") %>%
  filter(transcript_id %in% tx_order_masked) %>%
  pivot_longer(-transcript_id, names_to = c("sample", "rep"), names_sep = "_", values_to = "tx_exp") %>%
  separate(sample, c("cond", "time"), convert = TRUE, sep = 2) %>%
  dplyr::select(-cond) %>%
  group_by(transcript_id) %>%
  mutate(tx_exp = scale(tx_exp)) %>%
  ungroup() %>%
  mutate(transcript_id = factor(transcript_id, levels = tx_order_masked, ordered = TRUE),
         time = factor(time, levels = seq(0, 20, 4), ordered = TRUE)) %>%
  arrange(transcript_id)

df_gene_masked <- abundance_gene_no_batch %>%
  as.data.frame() %>% rownames_to_column("gene_id") %>%
  filter(gene_id %in% gene_block_masked) %>%
  pivot_longer(-gene_id, names_to = c("sample", "rep"), names_sep = "_", values_to = "gene_exp") %>%
  separate(sample, c("cond", "time"), convert = TRUE, sep = 2) %>%
  dplyr::select(-cond) %>%
  group_by(gene_id) %>%
  mutate(gene_exp = scale(gene_exp)) %>%
  ungroup() %>%
  mutate(gene_id = factor(gene_id, levels = gene_block_masked, ordered = TRUE),
         time = factor(time, levels = seq(0, 20, 4), ordered = TRUE)) %>%
  arrange(gene_id)

# cells outside the color scale's limits=c(-2.5,2.5) -- find which gene/sample
df_gene_masked %>% filter(abs(gene_exp) > 2.5) %>% left_join(gene2name, by = "gene_id") %>%
  write.csv(paste(out_dir, "csv/masked_genes_clipped_cells.csv", sep=""), row.names = FALSE)

panel_tx_masked <- ggplot(df_tx_masked) +
  geom_tile(aes(y = transcript_id, x = time, fill = tx_exp), width = 1) +
  facet_wrap(~rep, nrow = 1, labeller = as_labeller(function(x) paste("replicate", x))) +
  scale_fill_gradient2_tableau(limits = c(-2.5, 2.5), palette = "Green-Blue Diverging", guide = "none", oob = scales::squish) +
  new_scale_fill() +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = n_tx_masked, ymax = n_tx_masked * 1.05, fill = col),
            data = LD_bars, inherit.aes = FALSE) +
  scale_fill_manual(values = c("black", "grey80"), guide = "none") +
  ylab(paste(n_tx_masked, "significant isoforms")) +
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

panel_gene_masked <- ggplot(df_gene_masked) +
  geom_tile(aes(y = gene_id, x = time, fill = gene_exp), width = 1) +
  facet_wrap(~rep, nrow = 1, labeller = as_labeller(function(x) paste("replicate", x))) +
  scale_fill_gradient2_tableau(name = NULL, limits = c(-2.5, 2.5), palette = "Green-Blue Diverging",
                                oob = scales::squish,
                                guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                                        barwidth = unit(3, "cm"), barheight = unit(0.2, "cm"), reverse = TRUE)) +
  new_scale_fill() +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = n_gene_masked, ymax = n_gene_masked * 1.05, fill = col),
            data = LD_bars, inherit.aes = FALSE) +
  scale_fill_manual(values = c("black", "grey80"), guide = "none") +
  ylab(paste(n_gene_masked, "genes (not gene-significant)")) +
  xlab("Time (h)") + theme_minimal(base_size = 11, base_family = "Arial") +
  scale_x_discrete(expand = expansion()) + scale_y_discrete(expand = expansion()) +
  theme(text = element_text(family = "Arial"),
        axis.text.y = element_blank(), axis.line = element_blank(),
        axis.ticks.y = element_blank(), panel.spacing.x = unit(0.2, "lines"),
        legend.position = "top", panel.grid = element_blank(), panel.grid.major = element_blank(),
        strip.text = element_text(family = "Arial", size = 11),
        axis.text.x = element_text(family = "Arial", size = 9),
        axis.title.x = element_text(family = "Arial", size = 12),
        axis.title.y = element_text(family = "Arial", size = 12)) +
  theme(aspect.ratio = 4 * n_gene_masked / n_tx_masked)

heatmap_masked <- panel_tx_masked | panel_gene_masked
ggsave(paste(out_dir, "figures/masked_genes_heatmap.svg", sep=""), plot = heatmap_masked, width = 7, height = 5)

write.csv(gene_results %>% filter(gene_id %in% masked_genes),
          paste(out_dir, "csv/masked_genes_list.csv", sep=""), row.names = FALSE)

# switched version for the composite figure: genes left, isoforms right
# (matches gene_isoform_heatmap.svg's convention), same fitted panels, reordered
heatmap_masked_switched <- panel_gene_masked | panel_tx_masked
ggsave(paste(out_dir, "figures/masked_genes_heatmap_switched.svg", sep=""), plot = heatmap_masked_switched, width = 7, height = 5)

# composite-ready: squished via aspect.ratio (keeps text legible, unlike
# post-hoc image resizing). Legend dropped, top block's legend covers it.
composite_aspect <- 0.18
panel_gene_masked_composite <- panel_gene_masked + theme(aspect.ratio = composite_aspect, legend.position = "none")
panel_tx_masked_composite <- panel_tx_masked + theme(aspect.ratio = composite_aspect)
heatmap_masked_composite <- panel_gene_masked_composite | panel_tx_masked_composite
ggsave(paste(out_dir, "figures/masked_genes_heatmap_composite.svg", sep=""), plot = heatmap_masked_composite, width = 7, height = 2.2)
