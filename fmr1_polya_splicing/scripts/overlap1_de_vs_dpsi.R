#!/usr/bin/env Rscript
# =============================================================================
# Overlap #1: DESeq2 (gene-level) x SUPPA dPSI, all 7 event types SEPARATELY
# Author: Gricey
#
# Kiran's core task #3: "define the overlap between dPSI and deseq". Per
# Lies's call on the follow-up email, tested per event type, not pooled
# into one combined "any event type" bucket -- each event type gets its
# own independent overlap report.
#
# Thresholds match this project's established convention throughout:
# DESeq2 padj < 0.05 (already multiple-testing corrected across genes),
# SUPPA raw p < 0.05 per event (same threshold used in every SUPPA summary
# step in this project so far, e.g. run_suppa.sh's own p<0.05 check).
#
# Biological framing (confirmed with Lies): does missplicing precede
# misexpression? A gene showing up in both lists is stronger evidence than
# either alone; not every missliced gene is expected to be misexpressed.
# =============================================================================

BASE_DIR   <- "/home/grangel/fmr1_polya_splicing"
DESEQ_PATH <- file.path(BASE_DIR, "results", "normalisation", "deseq2_KO_vs_WT.tsv")
SUPPA_DIR  <- file.path(BASE_DIR, "results", "gse188840_suppa_2026-08-05", "diff")
OUTPUT_DIR <- file.path(BASE_DIR, "results", "overlap1")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

EVENT_TYPES <- c("SE", "A3", "A5", "MX", "RI", "AF", "AL")

# =============================================================================
# LOAD DESEQ2 RESULTS, GET SIGNIFICANT GENE SET
# =============================================================================
cat("Loading DESeq2 results...\n")
deseq <- read.table(DESEQ_PATH, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Total genes tested:", nrow(deseq), "\n")
cat("Columns:", paste(colnames(deseq), collapse = ", "), "\n")

deseq_sig <- deseq[!is.na(deseq$padj) & deseq$padj < 0.05, ]
deseq_sig_genes <- unique(deseq_sig$gene_id)
cat("DESeq2 significant genes (padj<0.05):", length(deseq_sig_genes), "\n\n")

# =============================================================================
# FOR EACH EVENT TYPE: load dpsi file, extract gene_id from the event_id,
# get the set of genes with >=1 significant event, compute overlap with
# DESeq2's significant gene set
# =============================================================================
summary_rows <- list()
overlap_long <- list()

for (EVENT in EVENT_TYPES) {
  dpsi_path <- file.path(SUPPA_DIR, paste0("diff_", EVENT, ".dpsi.temp.0"))
  cat("=== Processing", EVENT, "===\n")

  if (!file.exists(dpsi_path)) {
    cat("  WARNING: file not found, skipping:", dpsi_path, "\n\n")
    next
  }

  dpsi <- read.table(dpsi_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # Diagnostic print, first event type only, so we can confirm the column
  # layout actually matches what we assume below (event_id, dpsi, p-value)
  if (EVENT == EVENT_TYPES[1]) {
    cat("  Columns in dpsi file:", paste(colnames(dpsi), collapse = ", "), "\n")
    cat("  First row:", paste(as.character(dpsi[1, ]), collapse = " | "), "\n")
  }

  # event_id format: "ENSMUSG00000014353;SE:2:128676134-...:+" -- gene_id
  # is everything before the first semicolon
  event_id_col <- colnames(dpsi)[1]
  dpsi_col     <- colnames(dpsi)[2]
  pval_col     <- colnames(dpsi)[3]

  dpsi$gene_id <- sub(";.*", "", dpsi[[event_id_col]])
  dpsi$pvalue  <- suppressWarnings(as.numeric(dpsi[[pval_col]]))
  dpsi$dpsi_value <- suppressWarnings(as.numeric(dpsi[[dpsi_col]]))

  n_events_total <- nrow(dpsi)
  sig_events <- dpsi[!is.na(dpsi$pvalue) & dpsi$pvalue < 0.05, ]
  sig_genes_this_type <- unique(sig_events$gene_id)

  overlap_genes <- intersect(deseq_sig_genes, sig_genes_this_type)

  cat("  Total events:", n_events_total, "\n")
  cat("  Significant events (p<0.05):", nrow(sig_events), "\n")
  cat("  Unique genes with a significant", EVENT, "event:", length(sig_genes_this_type), "\n")
  cat("  Overlap with DESeq2 significant genes:", length(overlap_genes), "\n\n")

  summary_rows[[EVENT]] <- data.frame(
    event_type              = EVENT,
    n_events_total           = n_events_total,
    n_events_significant     = nrow(sig_events),
    n_genes_significant      = length(sig_genes_this_type),
    n_genes_overlap_with_DE  = length(overlap_genes)
  )

  if (length(overlap_genes) > 0) {
    overlap_detail <- sig_events[sig_events$gene_id %in% overlap_genes, ]
    overlap_detail <- merge(overlap_detail,
                             deseq_sig[, c("gene_id", "log2FoldChange", "padj")],
                             by = "gene_id")
    overlap_detail$event_type <- EVENT
    overlap_long[[EVENT]] <- overlap_detail[, c("gene_id", "event_type",
                                                  "log2FoldChange", "padj",
                                                  "dpsi_value", "pvalue")]
  }
}

# =============================================================================
# SAVE SUMMARY TABLE (one row per event type)
# =============================================================================
summary_table <- do.call(rbind, summary_rows)
write.table(summary_table,
            file      = file.path(OUTPUT_DIR, "overlap1_summary_by_event_type.tsv"),
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)
cat("\n=== SUMMARY ACROSS ALL EVENT TYPES ===\n")
print(summary_table)

# =============================================================================
# SAVE DETAILED OVERLAP GENE LIST (long format, all event types combined
# for convenience, but event_type column keeps them distinguishable --
# NOT a merged/deduplicated "any event type" list)
# =============================================================================
if (length(overlap_long) > 0) {
  overlap_table <- do.call(rbind, overlap_long)
  overlap_table <- overlap_table[order(overlap_table$padj, overlap_table$pvalue), ]
  write.table(overlap_table,
              file      = file.path(OUTPUT_DIR, "overlap1_gene_list.tsv"),
              sep       = "\t",
              quote     = FALSE,
              row.names = FALSE)
  cat("\nSaved:", file.path(OUTPUT_DIR, "overlap1_gene_list.tsv"), "\n")
  cat("Total (gene, event_type) overlap rows:", nrow(overlap_table), "\n")
  cat("Unique genes appearing in at least one event type's overlap:",
      length(unique(overlap_table$gene_id)), "\n")
} else {
  cat("\nNo overlaps found in any event type.\n")
}

cat("\n=== DONE ===\n")
cat("Summary:", file.path(OUTPUT_DIR, "overlap1_summary_by_event_type.tsv"), "\n")
cat("Gene list:", file.path(OUTPUT_DIR, "overlap1_gene_list.tsv"), "\n")
