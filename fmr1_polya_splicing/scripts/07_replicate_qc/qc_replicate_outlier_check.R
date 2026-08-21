#!/usr/bin/env Rscript
# =============================================================================
# QC check: is one of the 3 WT replicates a technical outlier?
#
# Checks whether a WT replicate is an outlier in our own oarfish/mm39
# gene-level quantification, using three standard DESeq2 QC diagnostics
# (PCA, sample-to-sample distance, Cook's distance).
#
# Rebuilds the gene-level dds object fresh (the main pipeline never saved
# it, only the results table) and saves it as .rds so it doesn't need
# rebuilding again for future QC questions.
# =============================================================================

library(tximport)
library(DESeq2)
library(ggplot2)

# ggsave()'s default SVG writer requires the "svglite" package, not
# installed in `biotools` -- save SVG via base R's grDevices::svg()
# device instead, no extra package needed.
save_plot <- function(plot, path_no_ext, width = 6, height = 5, dpi = 200) {
  ggsave(paste0(path_no_ext, ".png"), plot, width = width, height = height, dpi = dpi)
  grDevices::svg(paste0(path_no_ext, ".svg"), width = width, height = height)
  print(plot)
  grDevices::dev.off()
}

BASE_DIR    <- "/home/grangel/fmr1_polya_splicing"
OARFISH_DIR <- file.path(BASE_DIR, "results", "gse188840_oarfish_2026-08-04")
GTF_PATH    <- file.path(BASE_DIR, "data", "reference", "Mus_musculus.GRCm39.116.gtf.gz")
OUTPUT_DIR  <- file.path(BASE_DIR, "results", "qc_replicate_check")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

sample_names <- c("WT_Rep1", "WT_Rep2", "WT_Rep3", "KO_Rep1", "KO_Rep2", "KO_Rep3")
conditions   <- c("WT", "WT", "WT", "KO", "KO", "KO")
col_data <- data.frame(sample = sample_names,
                        condition = factor(conditions, levels = c("WT", "KO")),
                        row.names = sample_names)

quant_files <- file.path(OARFISH_DIR, paste0(sample_names, ".quant"))
names(quant_files) <- sample_names
missing <- quant_files[!file.exists(quant_files)]
if (length(missing) > 0) stop("Missing .quant files:\n", paste(missing, collapse = "\n"))

# tx2gene -- same code as the main pipeline
cat("Building tx2gene from GTF...\n")
gtf <- read.table(gzfile(GTF_PATH), sep = "\t", quote = "", comment.char = "#",
                   col.names = c("chr", "source", "feature", "start", "end",
                                 "score", "strand", "frame", "attributes"))
gtf_tx <- gtf[gtf$feature == "transcript", ]
extract_attr <- function(attr_string, key) {
  pattern <- paste0(key, ' "([^"]+)"')
  m <- regmatches(attr_string, regexpr(pattern, attr_string))
  if (length(m) == 0) return(NA)
  sub(paste0(key, ' "([^"]+)"'), "\\1", m)
}
tx2gene <- data.frame(tx_id = sapply(gtf_tx$attributes, extract_attr, key = "transcript_id"),
                       gene_id = sapply(gtf_tx$attributes, extract_attr, key = "gene_id"),
                       stringsAsFactors = FALSE, row.names = NULL)
tx2gene <- tx2gene[!is.na(tx2gene$tx_id) & !is.na(tx2gene$gene_id), ]

cat("Running tximport + DESeq2 (gene-level)...\n")
txi_gene <- tximport(quant_files, type = "oarfish", tx2gene = tx2gene,
                      txOut = FALSE, ignoreTxVersion = TRUE)
dds <- DESeqDataSetFromTximport(txi_gene, colData = col_data, design = ~condition)
dds <- DESeq(dds)
saveRDS(dds, file.path(OUTPUT_DIR, "dds_gene_level.rds"))

# =============================================================================
# 1. PCA -- does one WT replicate sit apart from the other two?
# =============================================================================
vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
pca_data$sample <- rownames(pca_data)

p1 <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1.2, size = 3, show.legend = FALSE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal(base_size = 12) +
  ggtitle("PCA, all 6 samples (vst-transformed)")
save_plot(p1, file.path(OUTPUT_DIR, "qc_pca"))

# =============================================================================
# 2. Sample-to-sample distance heatmap
# =============================================================================
sample_dists <- dist(t(assay(vsd)))
dist_matrix <- as.matrix(sample_dists)
dist_df <- as.data.frame(as.table(dist_matrix))
colnames(dist_df) <- c("sample1", "sample2", "distance")

p2 <- ggplot(dist_df, aes(sample1, sample2, fill = distance)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#0072B2") +
  geom_text(aes(label = round(distance, 1)), size = 3) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Sample-to-sample Euclidean distance (vst)")
save_plot(p2, file.path(OUTPUT_DIR, "qc_sample_distance"))

# =============================================================================
# 3. Cook's distance per sample -- DESeq2's own per-gene outlier diagnostic.
# A sample with systematically higher/wider Cook's distances looks unusual
# across many genes relative to the rest of its group.
# =============================================================================
cooks <- assays(dds)[["cooks"]]
cooks_df <- as.data.frame(cooks)
cooks_long <- stack(cooks_df)
colnames(cooks_long) <- c("cooks_distance", "sample")

p3 <- ggplot(cooks_long, aes(x = sample, y = cooks_distance)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_y_log10() +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Cook's distance per sample (log scale)")
save_plot(p3, file.path(OUTPUT_DIR, "qc_cooks_distance"))

# =============================================================================
# 4. Text summary -- the concrete numbers to actually read
# =============================================================================
cat("\n=== Median Cook's distance per sample (higher = more outlier-like) ===\n")
print(sort(sapply(cooks_df, median, na.rm = TRUE), decreasing = TRUE))

cat("\n=== Mean distance to same-condition replicates (higher = more isolated) ===\n")
for (s in sample_names) {
  cond <- col_data[s, "condition"]
  same_cond_samples <- setdiff(sample_names[conditions == cond], s)
  mean_dist <- mean(dist_matrix[s, same_cond_samples])
  cat(s, "(", as.character(cond), "): mean distance to same-condition replicates =",
      round(mean_dist, 2), "\n")
}

cat("\n=== DONE ===\n")
cat("Figures + dds_gene_level.rds saved in", OUTPUT_DIR, "\n")
cat("Read: does one WT replicate stand out with (a) a distinct PCA position,\n")
cat("(b) larger distance to the other 2 WT reps than KO reps show to each\n")
cat("other, and (c) higher median Cook's distance? If yes on multiple counts,\n")
cat("that's the candidate to consider excluding for a 2+3 comparison.\n")
