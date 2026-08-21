#!/usr/bin/env Rscript
# =============================================================================
# Splicing entropy by GO category -- reproduces Figure 4E of Chikhaoui,
# Mamgain, Seki et al., "Circadian PERIOD proteins sculpt the mammalian
# alternative splicing landscape" (bioRxiv 2024.12.23.630108, IGFL Lyon).
#
# Method: for every GO:BP term, a paired Wilcoxon signed-rank test compares
# entropy_WT vs entropy_KO across that term's member genes (paired by
# gene), BH-corrected across all terms tested, kept at p<0.05 & q<0.05
# (matching the paper's figure legend; NOTE the paper's own main text
# states a looser "adj p<0.1" instead -- a real inconsistency in their
# paper, worth confirming with Kiran if results ever look sparse).
#
# Dot plot: x = -log10(p), size = |delta entropy| (median WT-KO shift
# within the category), color = number of genes in category -- matching
# Figure 4E's own legend (not clusterProfiler's usual count=size/p=color
# convention).
# =============================================================================

suppressMessages({
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(GO.db)
  library(ggplot2)
})

BASE_DIR    <- "/home/grangel/fmr1_polya_splicing"
ENTROPY_DIR <- file.path(BASE_DIR, "results", "isoform_proportions")
OUTPUT_DIR  <- file.path(BASE_DIR, "results", "entropy_go_enrichment")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

MIN_GENES_PER_TERM <- 5    # standard minGSSize-style floor, avoids testing
                            # terms with too few genes for a meaningful Wilcoxon
P_CUTOFF <- 0.05
Q_CUTOFF <- 0.05

run_go_entropy <- function(threshold) {
  cat("########################################\n")
  cat("=== min", threshold, "TPM ===\n", sep = "")
  cat("########################################\n")

  ranked_path <- file.path(ENTROPY_DIR, paste0("entropy_shift_ranked_min", threshold, "TPM.tsv"))
  ranked <- read.table(ranked_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  ranked <- ranked[!is.na(ranked$entropy_WT) & !is.na(ranked$entropy_KO), ]

  # Map genes (non-versioned Ensembl IDs, same convention as the rest of
  # this project) to GO Biological Process terms
  gene_go <- AnnotationDbi::select(org.Mm.eg.db,
                                    keys    = unique(ranked$gene_id),
                                    keytype = "ENSEMBL",
                                    columns = c("GO", "ONTOLOGY"))
  gene_go <- gene_go[!is.na(gene_go$ONTOLOGY) & gene_go$ONTOLOGY == "BP" & !is.na(gene_go$GO), ]

  go_terms <- unique(gene_go$GO)
  cat("Testing", length(go_terms), "GO:BP terms across", nrow(ranked), "genes...\n")

  results_list <- lapply(go_terms, function(term) {
    genes_in_term <- unique(gene_go$ENSEMBL[gene_go$GO == term])
    sub <- ranked[ranked$gene_id %in% genes_in_term, ]
    if (nrow(sub) < MIN_GENES_PER_TERM) return(NULL)
    test <- tryCatch(wilcox.test(sub$entropy_WT, sub$entropy_KO, paired = TRUE),
                      error = function(e) NULL)
    if (is.null(test) || is.na(test$p.value)) return(NULL)
    data.frame(GO_id         = term,
               n_genes       = nrow(sub),
               delta_entropy = median(sub$entropy_WT - sub$entropy_KO, na.rm = TRUE),
               p_value       = test$p.value)
  })
  results <- do.call(rbind, results_list)
  if (is.null(results) || nrow(results) == 0) {
    cat("No GO terms had enough genes to test at min", threshold, "TPM.\n\n")
    return(NULL)
  }
  results$q_value <- p.adjust(results$p_value, method = "BH")

  term_names <- tryCatch(
    AnnotationDbi::select(GO.db, keys = results$GO_id, keytype = "GOID", columns = "TERM"),
    error = function(e) NULL
  )
  if (!is.null(term_names)) {
    results <- merge(results, term_names, by.x = "GO_id", by.y = "GOID", all.x = TRUE)
  } else {
    results$TERM <- results$GO_id  # fallback: plot GO IDs if GO.db lookup fails
  }

  write.table(results,
              file = file.path(OUTPUT_DIR, paste0("go_entropy_all_terms_min", threshold, "TPM.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  sig <- results[results$p_value < P_CUTOFF & results$q_value < Q_CUTOFF, ]
  cat(nrow(results), "terms tested,", nrow(sig), "pass p<", P_CUTOFF, "& q<", Q_CUTOFF, "\n", sep = "")

  if (nrow(sig) == 0) {
    cat("No significant GO terms at min", threshold, "TPM -- skipping plot.\n")
    cat("If this happens at all 3 thresholds, consider the looser main-text\n")
    cat("cutoff (adj p<0.1) the paper also mentions, or confirm with Kiran\n")
    cat("which cutoff he actually wants.\n\n")
    return(NULL)
  }

  sig$neg_log10_p     <- -log10(sig$p_value)
  sig$abs_delta_entropy <- abs(sig$delta_entropy)

  plot_df <- sig[order(sig$p_value), ][1:min(20, nrow(sig)), ]
  plot_df$TERM <- factor(plot_df$TERM, levels = rev(plot_df$TERM))

  p <- ggplot(plot_df, aes(x = neg_log10_p, y = TERM, size = abs_delta_entropy, color = n_genes)) +
    geom_point() +
    scale_color_gradient(low = "#B2182B", high = "#2166AC", name = "no. of genes") +
    scale_size_continuous(name = "delta entropy", range = c(2, 10)) +
    labs(x = expression(-log[10](italic(p)*"-value")), y = NULL,
         title = paste0("Splicing entropy by GO category, KO vs WT (min", threshold, "TPM)")) +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_text(size = 9))

  ggsave(file.path(OUTPUT_DIR, paste0("go_entropy_dotplot_min", threshold, "TPM.png")),
         p, width = 9, height = 7, dpi = 300)
  ggsave(file.path(OUTPUT_DIR, paste0("go_entropy_dotplot_min", threshold, "TPM.svg")),
         p, width = 9, height = 7)

  write.table(sig,
              file = file.path(OUTPUT_DIR, paste0("go_entropy_significant_min", threshold, "TPM.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  cat("Dot plot + significant-terms table saved.\n\n")
}

for (t in c(10, 20, 50)) run_go_entropy(t)

cat("=== DONE ===\n")
cat("Outputs in", OUTPUT_DIR, "\n")
