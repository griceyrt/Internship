#!/usr/bin/env Rscript
# =============================================================================
# Answers Lies's question (email 2026-08-10): how many total transcripts are
# "expressed" using the paper's own filter of >10 reads per condition?
#
# Reuses the exact same oarfish/tximport import as deseq2_isoform_level_test.R
# (same 6 samples, same txOut=TRUE), just stops after building the raw counts
# matrix instead of running DESeq2 -- this is a filtering-universe question,
# not a significance-testing one.
# =============================================================================

library(tximport)

BASE_DIR    <- "/home/grangel/fmr1_polya_splicing"
OARFISH_DIR <- file.path(BASE_DIR, "results", "gse188840_oarfish_2026-08-04")

sample_names <- c("WT_Rep1", "WT_Rep2", "WT_Rep3", "KO_Rep1", "KO_Rep2", "KO_Rep3")
conditions   <- c("WT", "WT", "WT", "KO", "KO", "KO")

quant_files <- file.path(OARFISH_DIR, paste0(sample_names, ".quant"))
names(quant_files) <- sample_names
missing <- quant_files[!file.exists(quant_files)]
if (length(missing) > 0) stop("Missing .quant files:\n", paste(missing, collapse = "\n"))

txi_tx <- tximport(quant_files, type = "oarfish", txOut = TRUE, ignoreTxVersion = TRUE)
counts <- txi_tx$counts
rownames(counts) <- sub("\\.[0-9]+$", "", rownames(counts))

wt_counts <- counts[, conditions == "WT"]
ko_counts <- counts[, conditions == "KO"]

# Two reasonable readings of ">10 reads per condition", report both since
# Lies's phrasing doesn't specify which -- pick whichever matches how the
# paper's own methods describe the filter once you check it.
wt_mean <- rowMeans(wt_counts)
ko_mean <- rowMeans(ko_counts)
wt_sum  <- rowSums(wt_counts)
ko_sum  <- rowSums(ko_counts)

n_mean_both  <- sum(wt_mean > 10 & ko_mean > 10)
n_mean_either <- sum(wt_mean > 10 | ko_mean > 10)
n_sum_both   <- sum(wt_sum > 10 & ko_sum > 10)
n_sum_either <- sum(wt_sum > 10 | ko_sum > 10)

cat("Total transcripts in the count matrix:", nrow(counts), "\n\n")
cat("Reading A -- mean per-replicate count > 10 in EACH condition:  ", n_mean_both, "\n")
cat("Reading B -- mean per-replicate count > 10 in EITHER condition:", n_mean_either, "\n")
cat("Reading C -- total (summed) count > 10 in EACH condition:      ", n_sum_both, "\n")
cat("Reading D -- total (summed) count > 10 in EITHER condition:    ", n_sum_either, "\n")
