#!/usr/bin/env Rscript
# =============================================================================
# Entropy-by-gene-category comparison
#
# Replicates the method from PMC10099729 (Ferrari et al., "Splicing
# complexity as a pivotal feature of alternative exons in mammalian
# species", BMC Genomics 2023): split genes into pre-defined categories
# and run a two-sided Wilcoxon rank-sum test on splicing entropy between
# the two groups. Applied here to entropy_diff (KO-WT) instead of their
# Whippet splice-path entropy.
#
# Categories implemented: housekeeping vs other, and expression level
# (high vs low). Gene age and tissue-specificity would need extra external
# reference data (not currently pulled in) -- add as new blocks following
# the same run_wilcoxon_summary() pattern if needed later.
#
# Housekeeping gene list: HRT Atlas mouse gene list (Hounkpe et al. 2020),
# pulled as the underlying .RData object from the atlas's GitHub repo
# (the site's own CSV download link is not a stable static URL).
# =============================================================================

BASE_DIR    <- "/home/grangel/fmr1_polya_splicing"
ENTROPY_DIR <- file.path(BASE_DIR, "results", "isoform_proportions")
HK_PATH     <- file.path(BASE_DIR, "data", "reference", "Housekeeping_Genes_Mouse.RData")
OUTPUT_DIR  <- file.path(BASE_DIR, "results", "entropy_by_category")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------
# Load housekeeping gene list from the .RData object (object name
# Mouse_HK_genes, symbol column "Gene", per HRT Atlas's own source code).
# Fails loudly with the available column names if that assumption is
# wrong, rather than silently mis-joining.
# -----------------------------------------------------------------------
load_hk_genes <- function(path) {
  if (!file.exists(path)) {
    stop("Housekeeping gene file not found at: ", path,
         "\nDownload it first:\n",
         "  wget -O ", path,
         " https://raw.githubusercontent.com/Bidossessih/HRT_Atlas/master/www/Housekeeping_Genes_Mouse.RData")
  }
  env <- new.env()
  loaded_names <- load(path, envir = env)
  cat("Objects found in RData file:", paste(loaded_names, collapse = ", "), "\n")

  if (!"Mouse_HK_genes" %in% loaded_names) {
    stop("Expected an object named 'Mouse_HK_genes' in the RData file but found: ",
         paste(loaded_names, collapse = ", "),
         "\nInspect the actual object(s) above and update this script's object name if it's changed.")
  }
  hk_df <- env$Mouse_HK_genes
  cat("Mouse_HK_genes columns found:", paste(colnames(hk_df), collapse = ", "), "\n")

  candidate_cols <- c("Gene", "Gene.name", "Gene_name", "GeneSymbol", "Gene.Symbol",
                       "Symbol", "gene_name", "gene_symbol")
  match_col <- intersect(candidate_cols, colnames(hk_df))
  if (length(match_col) == 0) {
    stop("Could not auto-detect a gene-symbol column in Mouse_HK_genes.\n",
         "Columns present: ", paste(colnames(hk_df), collapse = ", "), "\n",
         "Edit candidate_cols in this script to add the correct column name, then re-run.")
  }
  cat("Using column '", match_col[1], "' as the gene symbol for the HK join.\n", sep = "")
  unique(toupper(trimws(hk_df[[match_col[1]]])))
}

hk_genes <- load_hk_genes(HK_PATH)
cat("Loaded", length(hk_genes), "mouse housekeeping gene symbols from HRT Atlas.\n\n")

# -----------------------------------------------------------------------
# Run the comparison for each TPM threshold table (10/20/50), matching the
# rest of this project's convention of testing all three, min50TPM primary
# -----------------------------------------------------------------------
run_wilcoxon_summary <- function(df, group_col, label) {
  groups <- unique(na.omit(df[[group_col]]))
  if (length(groups) != 2) {
    cat("  [", label, "] skipped: expected exactly 2 groups, found ",
        length(groups), "\n", sep = "")
    return(NULL)
  }
  g1 <- df$entropy_diff[df[[group_col]] == groups[1]]
  g2 <- df$entropy_diff[df[[group_col]] == groups[2]]
  test <- wilcox.test(g1, g2)
  data.frame(
    category   = label,
    group1     = groups[1], n1 = length(g1), median1 = median(g1, na.rm = TRUE),
    group2     = groups[2], n2 = length(g2), median2 = median(g2, na.rm = TRUE),
    wilcox_p   = test$p.value
  )
}

for (THRESHOLD in c(10, 20, 50)) {
  cat("########################################\n")
  cat("=== min", THRESHOLD, "TPM ===\n", sep = "")
  cat("########################################\n")

  ranked_path <- file.path(ENTROPY_DIR, paste0("entropy_shift_ranked_min", THRESHOLD, "TPM.tsv"))
  ranked <- read.table(ranked_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # --- Category 1: housekeeping vs other ---
  ranked$hk_status <- ifelse(toupper(trimws(ranked$gene_name)) %in% hk_genes,
                              "housekeeping", "other")

  # --- Category 2: expression level, median split on this gene's own
  # average total TPM across conditions (adapted from the paper's fixed
  # TPM-50 cutoff, since our own distribution/depth differs) ---
  ranked$avg_total_TPM <- (ranked$total_TPM_WT + ranked$total_TPM_KO) / 2
  med_tpm <- median(ranked$avg_total_TPM, na.rm = TRUE)
  ranked$expr_level <- ifelse(ranked$avg_total_TPM >= med_tpm, "high_expr", "low_expr")

  results <- rbind(
    run_wilcoxon_summary(ranked, "hk_status",  "housekeeping_vs_other"),
    run_wilcoxon_summary(ranked, "expr_level", "high_vs_low_expression")
  )
  print(results)

  write.table(ranked,
              file      = file.path(OUTPUT_DIR, paste0("entropy_with_categories_min", THRESHOLD, "TPM.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(results,
              file      = file.path(OUTPUT_DIR, paste0("category_wilcoxon_min", THRESHOLD, "TPM.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n")
}

cat("=== DONE ===\n")
cat("Per-gene category labels + Wilcoxon summaries saved to", OUTPUT_DIR, "\n")
