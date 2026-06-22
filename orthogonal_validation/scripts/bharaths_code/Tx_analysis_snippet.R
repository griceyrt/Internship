library(magrittr)
library(tidyverse)
library(tximport)
library(DESeq2)
library(stageR)
library(RColorBrewer)
library(ggthemes)
library(ggh4x)
library(ggnewscale)
library(patchwork)
library(DEXSeq)
library(dplyr)
library(ggplot2)
library(ggforce)
library(sjmisc)
library(stringr)
library(viridis)
library(GenomicFeatures)
library(DRIMSeq)
library(scales)
library(reshape2)
library(patchwork)
library(rtracklayer)
library(lemon)
library(ggmosaic)

fdr_cutoff <- 0.05
base_size <- 9
retrieve_from_biomart <- TRUE

##-------------Functions-----------------------------------------------##
perGeneQValueHack <- function(object, p = "pvalue", method = perGeneQValueExact) {
  wTest <- which( !is.na( object$padj ) )
  ## use only those exons that were testable
  pvals     = object[[p]][wTest]
  ## 'factor' removes ununsed levels
  geneID    = factor(object[["groupID"]][wTest])
  geneSplit = split(seq(along=geneID), geneID)
  
  ## summarise p-values of exons for one gene: take the minimum
  pGene = sapply(geneSplit, function(i) min(pvals[i]))
  stopifnot(all(is.finite(pGene)))
  
  ## Determine the thetas to be used
  theta = unique(sort(pGene))
  
  ## compute q-values associated with each theta
  q = method(pGene, theta, geneSplit)
  
  ## return a named vector of q-values per gene
  res        = rep(NA_real_, length(pGene))
  res        = q[match(pGene, theta)]
  res = pmin(1, res)
  names(res) = names(geneSplit)
  stopifnot(!any(is.na(res)))
  return(res)
}

perGeneQValueExact = function(pGene, theta, geneSplit) {
  stopifnot(length(pGene)==length(geneSplit))
  
  ## Compute the numerator \sum_{i=1}^M 1-(1-theta)^{n_i}
  ## Below we first identify the summands which are the same
  ## (because they have the same n_i), then do the sum via the
  ## mapply
  numExons     = listLen(geneSplit)
  tab          = tabulate(numExons)
  notZero      = (tab>0)
  numerator    = mapply(function(m, n) m * (1 - (1-theta)^n),
                        m = tab[notZero],
                        n = which(notZero))
  numerator    = rowSums(numerator)
  
  ## Compute the denominator: for each value of theta, the number
  ## of genes with pGene <= theta[i].
  ## Note that in cut(..., right=TRUE), the intervals are
  ## right-closed (left open) intervals.
  bins   = cut(pGene, breaks=c(-Inf, as.vector(theta)), right = TRUE, include.lowest = TRUE)
  counts = tabulate(bins, nbins = nlevels(bins))
  denom  = cumsum(counts)
  stopifnot(denom[length(denom)]==length(pGene))
  
  return(numerator/denom)
}


txdb <- makeTxDbFromGFF("../results/flair_2022-01-05/transcriptome_ext.gtf", format = "gtf")

annotation <- transcriptLengths(txdb, with.cds_len = TRUE, 
                                with.utr3_len = TRUE, with.utr5_len = TRUE) %>%
  dplyr::select(-tx_id) %>%
  dplyr::rename(all_of(c(transcript_id="tx_name")))

if (retrieve_from_biomart) {
  ensembl <- useMart("ensembl", "mmusculus_gene_ensembl", host = "https://sep2019.archive.ensembl.org")
  filters <- listFilters(ensembl)
  attribs <- listAttributes(ensembl)
  
  mart <- getBM(filters = "ensembl_transcript_id",
                attributes = c("ensembl_transcript_id", "ensembl_gene_id", "mgi_symbol", "transcript_biotype"),
                values = annotation$tx_name,
                mart = ensembl)
  
  colnames(mart) <- sub("^ensembl_","", colnames(mart))
  mart %<>% distinct(transcript_id, .keep_all = TRUE)
  write_rds(mart, paste0("../results/biomart_annotation_", format(Sys.time(), "%Y-%m-%d.rds")), compress="bz2")
} else {
  mart <- read_rds("../results/biomart_annotation_2023-02-11.rds")
}

annotation %<>% left_join(mart) %>%
  mutate(class=case_when(grepl("^ENS", transcript_id) ~ "Annotated", 
                         grepl("ENS", gene_id) ~ "Novel-Annotated", 
                         .default = "Novel-Novel"))

rownames(annotation) <- annotation$transcript_id

##-------------Loading Salmon data-------------------------------------##
files <- list.files(path = "../results/salmon_2023-02-11", pattern = "quant_CT", 
                    full.names = TRUE, include.dirs = FALSE, no.. = TRUE, recursive = TRUE)
names(files) <- str_extract(files, "(CT\\d+_\\d)(\\.sf)$", group = 1)

txi.salmon <- tximport(files, type = "salmon", txOut = TRUE, countsFromAbundance = "no", ignoreTxVersion = TRUE)

exp_design <- data.frame(batch = factor(str_extract(colnames(txi.salmon$abundance), "(?<=_)\\d")),
                         time = as.numeric(str_extract(colnames(txi.salmon$abundance), "\\d+")),
                         sample_id= colnames(txi.salmon$abundance))
rownames(exp_design) <- exp_design$sample_id


counts <- magrittr::set_colnames(txi.salmon$counts, exp_design$sample_id)
rownames(counts) <- str_extract(rownames(counts), "([A-Za-z0-9-]+)")
mode(counts) <- "integer"

##---------Expression levels in  transcriptome--------------------##
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = exp_design,
                              design= ~ batch)

rowData(dds) <- annotation[rownames(dds),]

keep <- rowSums(counts(dds)) > 20 #edgeR::filterByExpr(dds[, colData(dds)$cond=="WT"], min.count=10, min.prop=0.5)

dds <- dds[keep, ]

abundance_no_batch_effect <- txi.salmon$abundance %>% data.frame %>%
  magrittr::set_rownames(str_extract(rownames(.), "([A-Za-z0-9-]+)")) %>% 
  {log2(. + 0.001)} %>%
  {limma::removeBatchEffect(.[rownames(dds), ], batch = colData(dds)$batch)}

##---------Circadian analysis of transcriptome--------------------##
##-- modify this with appropriate design matrix to do DE analysis--##

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = exp_design[exp_design$cond == "WT",],
                              design= ~ batch + inphase + outphase)

row.names(tx2gene) <- tx2gene$transcript_stable_id

rowData(dds) <- tx2gene[rownames(dds),]

keep <- rowSums(counts(dds)) > 20

dds <- dds[keep, ]

dds <- DESeq(dds, reduced = ~batch, test = "LRT", fitType = "local")

res <- DESeq2::results(dds, independentFiltering = TRUE, alpha = fdr_cutoff)

res$groupID <- rowData(dds)$gene_stable_id

res <- cbind(res, coef(dds))

pScreen <- perGeneQValueHack(res)
pConfirmation <- matrix(res$pvalue, ncol = 1)

dimnames(pConfirmation) <- list(rowData(dds)$transcript_stable_id, "transcript")

stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=TRUE, tx2gene=tx2gene)

stageRObj <- stageWiseAdjustment(stageRObj, method="dte", alpha=fdr_cutoff, allowNA=TRUE)

suppressWarnings({
  rhy_padj <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=FALSE)
})

