#!/usr/bin/env Rscript
# Author: Gricey
#
# Runs clusterProfiler::enrichGO (Biological Process) on each significant-
# splicing-gene list (4 timepoints vs WT_SD + union). Parameters mirror
# Desmond's own DE-gene GO enrichment (rna_seq_pipeline.r go_enrichment())
# for comparability: ont="BP", minGSSize=1, pvalueCutoff=0.05, qvalueCutoff=1,
# readable=TRUE. One deviation: enrichGO is called directly with
# keyType="ENSEMBL" rather than converting symbol->ENTREZID first, since my
# gene lists are already Ensembl IDs.
#
# USAGE: Rscript 11_go_enrichment.R <scratch_dir>  (called by 11_go_enrichment.sh)

suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Mm.eg.db))

args <- commandArgs(trailingOnly = TRUE)
scratch_dir <- args[1]

sig_dir <- file.path(scratch_dir, "results/suppa/significant")
out_dir <- file.path(scratch_dir, "results/suppa/go_enrichment")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

comparisons <- c(
  "W8_WHFD_vs_WT_SD",
  "W18_WHFD_vs_WT_SD",
  "W24_WHFD_vs_WT_SD",
  "W42_WHFD_vs_WT_SD",
  "union_all_comparisons"
)

for (comp in comparisons) {
  gene_file <- file.path(sig_dir, paste0("significant_genes_", comp, ".txt"))
  if (!file.exists(gene_file)) {
    cat("Skipping", comp, "- file not found:", gene_file, "\n")
    next
  }

  genes <- readLines(gene_file)
  cat("=== GO enrichment:", comp, "(", length(genes), "genes) ===\n")

  enrich <- tryCatch({
    enrichGO(
      gene = genes,
      OrgDb = org.Mm.eg.db,
      keyType = "ENSEMBL",
      ont = "BP",
      minGSSize = 1,
      pvalueCutoff = 0.05,
      qvalueCutoff = 1,
      readable = TRUE
    )
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    NULL
  })

  if (is.null(enrich) || nrow(enrich@result) == 0) {
    cat("  No enrichment result for", comp, "\n")
    next
  }

  result_df <- enrich@result
  out_csv <- file.path(out_dir, paste0("go_enrichment_", comp, ".csv"))
  write.csv(result_df, out_csv, row.names = FALSE)

  n_sig <- sum(result_df$p.adjust < 0.05)
  cat("  Wrote", out_csv, "-", nrow(result_df), "terms tested,",
      n_sig, "significant (padj<0.05)\n")
}

cat("\nDone.\n")
