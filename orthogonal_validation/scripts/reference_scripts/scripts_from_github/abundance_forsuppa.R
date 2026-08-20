#This is to generate no batch files to use in suppa
setwd("C:/Users/kmamgain/prettyCoolPores_new/scripts") 

#load required packages
library(magrittr)
library(tximport)
library(DESeq2)
library(stageR)
library(RColorBrewer)
library(ggthemes)
library(ggh4x)
library(ggnewscale)
library(biomaRt)
library(sjmisc)
library(scales)
library(tidyverse)
library(stringr)
library(rlang)
library(GenomicFeatures)

count_threshold<- 5 #threshold for reads
retrieve_from_biomart<- FALSE #do you want to retrive annotation from biomart? not needed with biomart annotation in data folder

txdb <- makeTxDbFromGFF("../results/FLAIR_2024-09-11/transcriptome_productivity.gtf", format = "gtf") #make txdb object
annotation <- transcriptLengths(txdb, with.cds_len = TRUE, 
  with.utr3_len = TRUE, with.utr5_len = TRUE) %>% 
  dplyr::select(-tx_id) %>%
  dplyr::rename(transcript_id="tx_name")
if (retrieve_from_biomart) {
  ensembl <- useMart("ensembl", "mmusculus_gene_ensembl", host = "https://sep2019.archive.ensembl.org")
  filters <- listFilters(ensembl)
  attribs <- listAttributes(ensembl)
  mart <- getBM(filters = "ensembl_transcript_id",
                attributes = c("ensembl_transcript_id", "ensembl_gene_id", "mgi_symbol", "entrezgene_id", "transcript_biotype"),
                values = annotation$transcript_id,
                mart = ensembl)
  
  colnames(mart) <- sub("^ensembl_","", colnames(mart))
  mart %<>% distinct(transcript_id, .keep_all = TRUE)
  write_rds(mart, paste0("../data/biomart_annotation_", format(Sys.time(), "%Y-%m-%d.rds")), compress="bz2")
} else {
  mart <- read_rds("../data/biomart_annotation_2024-11-05.rds")
}
annotation %<>% left_join(mart) %>%
  mutate(class=case_when(grepl("^ENS", transcript_id) ~ "Annotated", 
                         grepl("ENS", gene_id) ~ "Novel-Annotated", 
                         .default = "Novel-Novel")) %>%
  group_by(gene_id) %>% 
  mutate(mgi_symbol = ifelse(grepl("^ENS", gene_id), 
                             unique(na.omit(mgi_symbol)), mgi_symbol)) %>%
  as.data.frame
rownames(annotation) <- annotation$transcript_id

#load salmon quantification files
files <- list.files(path = "../results/salmon_2024-09-25/", pattern = ".sf", 
                    full.names = TRUE, include.dirs = FALSE, no.. = TRUE, recursive = TRUE)
names(files) <- str_extract(files, "([CK][TO]\\d+_\\d)(\\.sf)$", group = 1)  #extract name of files to keep condition and CT
tx2gene=annotation[c("transcript_id", "gene_id")]

BPPARAM <-  BiocParallel::MulticoreParam(2)
txi.salmon <- tximport(files, type = "salmon", txOut = TRUE,
                       countsFromAbundance = "no",
                       ignoreTxVersion = TRUE)   #importing salmon files

exp_design <- data.frame(batch = factor(str_extract(colnames(txi.salmon$abundance), "(?<=_)\\d")),
                         time = as.numeric(str_extract(colnames(txi.salmon$abundance), "\\d+")),
                         cond = factor(str_extract(colnames(txi.salmon$abundance), "^[A-Z]+")),
                         sample_id= colnames(txi.salmon$abundance)) ## df for experimental design
exp_design[exp_design$cond == "KO",]$batch = 2 #All KO is batch 2
rownames(exp_design) <- colnames(counts)
exp_design$cond <- fct_recode(exp_design$cond, "WT" = "CT", "Per_KO" = "KO")

#exp_design$inphase = cos(2 * pi * exp_design$time/24) #remove
#exp_design$outphase = sin(2 * pi * exp_design$time/24) #remove

counts <- magrittr::set_colnames(txi.salmon$counts, exp_design$sample_id)
mode(counts) <- "integer"

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = exp_design,
                              design= ~ batch + cond )

row.names(tx2gene) <- tx2gene$transcript_id

rowData(dds) <- tx2gene[rownames(dds),]

keep <- rowSums(counts > count_threshold) >= 2

dds <- dds[keep, ]

#dds <- DESeq(dds,  reduced = ~batch + cond, "LRT", fitType = "local")

abundance = varianceStabilizingTransformation(dds, fitType = "local")
abundance_no_batch_effect <- limma::removeBatchEffect(assay(abundance),
                                                      batch = colData(dds)$batch,
                                                      design = model.matrix(~cond , data = colData(dds)))
abundance_no_batch_effect_no_log <- 2^abundance_no_batch_effect

##saving the abundance table for suppa ioe##
not_expressed<- data.frame(x=setdiff(annotation$transcript_id, rownames(abundance_no_batch_effect_no_log )))
rownames(not_expressed)<- not_expressed$x
not_expressed[1:16]<-0
colnames(not_expressed)<-colnames(abundance_no_batch_effect_no_log)
abundance_no_log_forsuppa<-rbind(abundance_no_batch_effect_no_log,not_expressed)
#for  run_suppa.sh#
write.table(abundance_no_log_forsuppa %>% as.matrix , "../data/abundance_no_batch_for_suppa.txt", sep="\t", row.names=TRUE, quote= FALSE)

#if also want gene name and symbol
#abundance_no_log_forsuppa$geneid<- annotation$gene_id[match(rownames(abundance_no_log_forsuppa), annotation$transcript_id)]
#abundance_no_log_forsuppa$mgi<- annotation$mgi_symbol[match(rownames(abundance_no_log_forsuppa), annotation$transcript_id)]
#write.table(abundance_no_log_forsuppa %>% as.matrix , "../data/abundance_no_batch_for_suppa.txt", sep="\t", row.names=TRUE, quote= FALSE)
#for differential_suppa.sh#
write.table(abundance_no_log_forsuppa %>% dplyr::select(c(7:8)) %>% as.matrix() , "../data/WT_abundance.txt", sep="\t", row.names=TRUE, quote= FALSE)
write.table(abundance_no_log_forsuppa %>% dplyr::select(c(13:14)) %>% as.matrix() , "../data/KO_abundance.txt", sep="\t", row.names=TRUE, quote= FALSE)