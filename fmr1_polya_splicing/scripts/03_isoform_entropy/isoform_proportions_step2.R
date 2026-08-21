#!/usr/bin/env Rscript
# =============================================================================
# Isoform proportion stacked-bar figure -- STEP 2 of 4
#
# Builds on step 1's avg_tpm_per_condition.tsv: joins in gene_id (from the
# GTF), then for each gene computes each isoform's share of that gene's
# total TPM, separately for WT and KO, keeping only isoforms present
# (TPM>0) in that condition.
# =============================================================================

BASE_DIR   <- "/home/grangel/fmr1_polya_splicing"
GTF_PATH   <- file.path(BASE_DIR, "data", "reference", "Mus_musculus.GRCm39.116.gtf")
INPUT_PATH <- file.path(BASE_DIR, "results", "isoform_proportions", "avg_tpm_per_condition.tsv")
OUTPUT_DIR <- file.path(BASE_DIR, "results", "isoform_proportions")

# =============================================================================
# LOAD STEP 1 OUTPUT
# =============================================================================
cat("Loading Step 1 output...\n")
avg_tpm <- read.table(INPUT_PATH, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Transcripts loaded:", nrow(avg_tpm), "\n")

# =============================================================================
# BUILD tx2gene FROM GTF (also pulls gene_name for readable figure labels)
# =============================================================================
cat("Building tx2gene from GTF...\n")
gtf <- read.table(GTF_PATH, sep = "\t", quote = "", comment.char = "#",
                  col.names = c("chr","source","feature","start","end",
                                "score","strand","frame","attributes"))
gtf_tx <- gtf[gtf$feature == "transcript", ]

extract_attr <- function(attr_string, key) {
  pattern <- paste0(key, ' "([^"]+)"')
  m <- regmatches(attr_string, regexpr(pattern, attr_string))
  if (length(m) == 0) return(NA)
  sub(paste0(key, ' "([^"]+)"'), "\\1", m)
}

tx2gene <- data.frame(
  transcript_id = sapply(gtf_tx$attributes, extract_attr, key = "transcript_id"),
  gene_id       = sapply(gtf_tx$attributes, extract_attr, key = "gene_id"),
  gene_name     = sapply(gtf_tx$attributes, extract_attr, key = "gene_name"),
  stringsAsFactors = FALSE,
  row.names = NULL
)
tx2gene <- tx2gene[!is.na(tx2gene$transcript_id) & !is.na(tx2gene$gene_id), ]
cat("tx2gene built:", nrow(tx2gene), "transcripts.\n")

# =============================================================================
# JOIN gene_id ONTO THE TPM TABLE
# =============================================================================
cat("Joining gene_id onto TPM table...\n")
merged <- merge(avg_tpm, tx2gene, by = "transcript_id", all.x = TRUE)
cat("Transcripts with a matched gene_id:", sum(!is.na(merged$gene_id)), "\n")
cat("Transcripts with NO matched gene_id (dropped):", sum(is.na(merged$gene_id)), "\n")
merged <- merged[!is.na(merged$gene_id), ]

# =============================================================================
# COMPUTE PER-GENE TOTALS AND PER-ISOFORM PROPORTIONS, PER CONDITION
# =============================================================================
cat("Computing per-gene totals...\n")
gene_totals_WT <- tapply(merged$WT, merged$gene_id, sum)
gene_totals_KO <- tapply(merged$KO, merged$gene_id, sum)

merged$gene_total_WT <- gene_totals_WT[merged$gene_id]
merged$gene_total_KO <- gene_totals_KO[merged$gene_id]

merged$proportion_WT <- ifelse(merged$gene_total_WT > 0, merged$WT / merged$gene_total_WT, NA)
merged$proportion_KO <- ifelse(merged$gene_total_KO > 0, merged$KO / merged$gene_total_KO, NA)

# =============================================================================
# RESHAPE TO LONG FORMAT: one row per (gene, transcript, condition) where
# that transcript is actually present (TPM>0) in that condition
# =============================================================================
cat("Reshaping to long format (one row per isoform per condition present)...\n")

long_WT <- merged[merged$present_in_WT, c("gene_id","gene_name","transcript_id","WT","proportion_WT")]
colnames(long_WT) <- c("gene_id","gene_name","transcript_id","TPM","proportion")
long_WT$condition <- "WT"

long_KO <- merged[merged$present_in_KO, c("gene_id","gene_name","transcript_id","KO","proportion_KO")]
colnames(long_KO) <- c("gene_id","gene_name","transcript_id","TPM","proportion")
long_KO$condition <- "KO"

long_table <- rbind(long_WT, long_KO)
long_table <- long_table[order(long_table$gene_id, long_table$condition, -long_table$proportion), ]

cat("Total rows in long-format table:", nrow(long_table), "\n")
cat("Unique genes represented:", length(unique(long_table$gene_id)), "\n")

write.table(long_table,
            file      = file.path(OUTPUT_DIR, "isoform_proportions_long.tsv"),
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

# =============================================================================
# ALSO SAVE: how many isoforms are present per gene, per condition
# =============================================================================
count_WT <- aggregate(transcript_id ~ gene_id, data = long_table[long_table$condition == "WT", ], FUN = length)
colnames(count_WT) <- c("gene_id", "n_isoforms_WT")
count_KO <- aggregate(transcript_id ~ gene_id, data = long_table[long_table$condition == "KO", ], FUN = length)
colnames(count_KO) <- c("gene_id", "n_isoforms_KO")

isoform_counts <- merge(count_WT, count_KO, by = "gene_id", all = TRUE)
isoform_counts$n_isoforms_WT[is.na(isoform_counts$n_isoforms_WT)] <- 0
isoform_counts$n_isoforms_KO[is.na(isoform_counts$n_isoforms_KO)] <- 0
isoform_counts$isoform_count_diff <- isoform_counts$n_isoforms_KO - isoform_counts$n_isoforms_WT

write.table(isoform_counts,
            file      = file.path(OUTPUT_DIR, "isoform_counts_per_gene.tsv"),
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

cat("\n=== STEP 2 DONE ===\n")
cat("Saved:", file.path(OUTPUT_DIR, "isoform_proportions_long.tsv"), "\n")
cat("Saved:", file.path(OUTPUT_DIR, "isoform_counts_per_gene.tsv"), "\n")
cat("Genes with a different isoform count between WT and KO:",
    sum(isoform_counts$isoform_count_diff != 0), "\n")
