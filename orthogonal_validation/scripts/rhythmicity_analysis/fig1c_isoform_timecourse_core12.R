#!/usr/bin/env Rscript
# Fig 1c redo, curated version -- restricted to the 12 core circadian clock
# genes, one page, matching the original manuscript's single-figure scale.
# Author: Gricey
# Run from orthogonal_validation/: Rscript scripts/fig1c_isoform_timecourse_core12.R

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
out_dir <- "results/rhythmicity_analysis/gene_rhythmicity/"

core_clock_genes <- c("Arntl", "Clock", "Per1", "Per2", "Cry1", "Cry2",
                       "Nr1d1", "Nr1d2", "Dbp", "Rorc", "Ciart", "Bhlhe40")

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

gene_list <- gene2name %>% filter(gene_name %in% core_clock_genes) %>% pull(gene_id)
length(gene_list)  # should be 12

tx_ids <- tx2gene %>% filter(gene_id %in% gene_list) %>% pull(transcript_id)

files <- list.files(path = salmon_path, pattern = "\\.sf$", full.names = TRUE, recursive = TRUE)
names(files) <- str_extract(files, "([CKTO]+\\d+_\\d)(?=\\.sf$)")
wt_files <- files[grepl("^CT", names(files))]

txi <- tximport(wt_files, type = "salmon", txOut = TRUE, countsFromAbundance = "no",
                 ignoreTxVersion = TRUE, dropInfReps = TRUE, importer = function(x) read.delim(x))

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

df <- abundance_no_batch %>%
  as.data.frame() %>% rownames_to_column("transcript_id") %>%
  pivot_longer(-transcript_id, names_to = c("sample", "rep"), names_sep = "_", values_to = "value") %>%
  separate(sample, c("cond", "time"), convert = TRUE, sep = 2) %>%
  dplyr::select(-cond) %>%
  left_join(tx2gene, by = "transcript_id") %>%
  left_join(gene2name, by = "gene_id") %>%
  left_join(productivity, by = "transcript_id") %>%
  mutate(productivity = ifelse(is.na(productivity), "non-coding", productivity),
         gene_name = factor(gene_name, levels = core_clock_genes))

p <- ggplot(data = df) +
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

ggsave(paste(out_dir, "figures/fig1c_isoform_timecourse_core12.svg", sep = ""), plot = p, width = 10, height = 7.5)
