library(tximport)
library(GenomicFeatures)
library(GenomicRanges)
library(DESeq2)
library(DRIMSeq)
library(DEXSeq)
library(tidyverse)
library(stageR)
library(biomaRt)
library(limma)
library(DEXSeq)
library(ggplot2)
library(ggbeeswarm)
library(dplyr)
library(cowplot)
library(dplyr)
library(patchwork)
library(VennDiagram)
library(ggrepel)
setwd("C:/Users/kmamgain/prettycoolPores_new/scripts")   #change here
expressed_t<-read.csv("../results/expr_tx_list.csv")

mycolors_updown <- c("purple", "grey", "#FFC94C")
names(mycolors_updown) <- c("DOWN","NO","UP")
mycolors_class <- c("#66c2a5", "#fc8d62", "#8da0cb")
names(mycolors_class) <- c("Annotated","Novel-Annotated","Novel-Novel")

###
fig_dir <- paste(getwd(), "/../Manuscript_figures_",Sys.Date(), sep="") #Path of directory where we want to save all the extra figures
theme_legend<- theme(legend.key.size = unit(0.1, 'cm'), #change legend key size
                     legend.key.height = unit(0.5, 'cm'), #change legend key height
                     legend.key.width = unit(0.1, 'cm'), #change legend key width
                     legend.title = element_blank())
flair_path<- "../results/FLAIR_2024-09-11"
biomart_path<-"../data/biomart_annotation_2024-11-05.rds"
salmon_path<-"../results/salmon_2024-09-25/"
suppa_path<- "../results/suppa_2024-11-06"
diff_suppa_path<- "../results/Differential_suppa_2024-11-07"
list_dir<- "../results/lists"


fdr_cutoff <- 0.05
count_threshold <- 5
base_size<-9
retrieve_from_biomart <- FALSE
tx_names_with_version <- TRUE
mycolors_updown <- c("DOWN" = "purple", "NO"="grey", "UP" ="#FFC94C")
mycolors_class <- c("Novel-Novel"= "#8da0cb","Novel-Annotated"= "#fc8d62","Annotated"= "#66c2a5")
scale_fill <- scale_fill_manual(values = c("#BDD7E7", "#6BAED6", "#3182BD", "#08519C", "#EFF3FF"))
scale_color <- scale_color_manual(values = c("#BDD7E7", "#6BAED6", "#3182BD", "#08519C", "#EFF3FF"))

#thresholds
base_size <- 9
fdr_cutoff <- 0.05
count_threshold <- 5
retrieve_from_biomart<- FALSE
####
txdb <- makeTxDbFromGFF(paste(flair_path,"/transcriptome_productivity.gtf", sep=""), format = "gtf")
retrieve_from_biomart<- FALSE
circadian_transcripts<- read.csv(paste(list_dir, "/circadian_transcripts.csv", sep=""))

annotation <- transcripts(txdb, columns=c("gene_id", "tx_name")) %>% 
  mcols() %>%
  as.data.frame %>% 
  dplyr::rename(transcript_id="tx_name") %>%
  mutate(gene_id = unlist(gene_id)) %>%
  magrittr::set_rownames(.$transcript_id)

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
  write_rds(mart, paste0("../results/biomart_annotation_", format(Sys.time(), "%Y-%m-%d.rds")), compress="bz2")
} else {
  mart <- read_rds(biomart_path)
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


####---------------------------------srsf1KO------------------------#####
####---------------------------aligned to transcriptome---------------#####
files <- list.files(path = "../results/SRP252896_salmon_2024-12-06/", 
                            pattern = "*.sf", full.names = TRUE, include.dirs = TRUE, no.. = TRUE, recursive = TRUE)
names(files) <- str_extract(files, "SRR\\d+")
          
txi_SRP252896 <- tximport(files, type = "salmon", txOut = TRUE, countsFromAbundance = "no", ignoreTxVersion = TRUE)
          
meta_data_SRP252896 <- read_csv("../meta/SRP252896_SraRunTable.txt") %>% 
            dplyr::select(!!!c("Run", "Genotype")) %>% 
            dplyr::mutate(Genotype = ifelse(grepl("^srsf1", Genotype), "SRSF1KO", "WT"),
                          Genotype = fct(Genotype, levels=c("WT", "SRSF1KO"))) %>%
            filter(Run %in% colnames(txi_SRP252896$abundance))
counts <- txi_SRP252896$counts
mode(counts) <- "integer"
srsf1_dds <- DESeqDataSetFromMatrix(countData = counts,
                                              colData = meta_data_SRP252896,
                                              design= ~ Genotype)
rowData(srsf1_dds) <- annotation[rownames(srsf1_dds),]
keep <- rownames(counts) %in% expressed_t$transcript_id
srsf1_dds <- srsf1_dds[keep, ]
srsf1_dds <- DESeq(srsf1_dds, fitType = "local")
srsf1_res <- DESeq2::results(srsf1_dds, 
                             contrast= c("Genotype", "SRSF1KO", "WT"),
                             independentFiltering = TRUE, 
                             alpha = fdr_cutoff)
srsf1_res <- lfcShrink(srsf1_dds, 
                       coef="Genotype_SRSF1KO_vs_WT", 
                       res = srsf1_res,
                       type = "apeglm")
          
results_SRSF1K0 <- srsf1_res %>% as.data.frame() %>%
            filter(!is.na(padj)) %>%
            rownames_to_column(var="transcript_id") %>%
            dplyr::select(c("transcript_id", "log2FoldChange", "padj")) %>%
            dplyr::rename(padj_srsf1ko = "padj", log2FC_srsf1ko = "log2FoldChange") %>% left_join(., annotation) %>%
            mutate(diffexp= ifelse(.$log2FC_srsf1ko > 1.5 & .$padj_srsf1ko < 0.05, "UP", ifelse(.$log2FC_srsf1ko < - 1.5 & .$padj_srsf1ko < 0.05,"DOWN","No")),
                   cir=case_when(!is.na((match(.$transcript_id, circadian_transcripts$transcript_id))) ~ "Yes",
                                 .default = "No"))

srsf1_DE_summary <- results_SRSF1K0 %>% group_by(diffexp) %>% 
  summarize(gene = n_distinct(gene_id), 
            tx = n_distinct(transcript_id)) 

interesting_srsf1<- c("Npas2","Ciart", "Usp2")

srsf1p_common <- ggplot(data=results_SRSF1K0, aes(x=log2FC_srsf1ko, y=-log10(padj_srsf1ko), col=class)) + 
                      geom_point(size=0.4) + theme_minimal()+
                      geom_vline(xintercept=c(-1.5, 1.5), col="black", linetype = "dashed") +
                      geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed") + 
                      theme_classic() + labs(color="Annotation") +
                      theme(text=element_text(size= 10),strip.text.x = element_blank(), legend.position="bottom", 
                            legend.box.spacing = unit(0, "pt")) +
                      scale_colour_manual(values = mycolors_class) + 
                      xlab('log2(Expression Fold Change) WT vs SRSF1KO') +
                      ylab("-log10(padj)") 

srsf1p_cir <- ggplot(data=results_SRSF1K0, aes(x=log2FC_srsf1ko, y=-log10(padj_srsf1ko), col=cir)) + 
              geom_point(size=0.4) + theme_minimal()+
              geom_vline(xintercept=c(-1.5, 1.5), col="black", linetype = "dashed") +
              geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed") + 
              theme_classic() + labs(color="Circadian") +
              geom_text_repel(data=results_SRSF1K0 %>% filter(cir == "Yes")%>% filter(diffexp != "No") %>% filter(mgi_symbol %in% interesting_srsf1),
                  aes(y=-log10(padj_srsf1ko), x= log2FC_srsf1ko,label= mgi_symbol), col= "black", max.overlaps=20) +
              theme(text=element_text(size= 10),strip.text.x = element_blank(), legend.position="bottom",legend.box.spacing = unit(0, "pt")) +
              xlab('log2(Expression Fold Change) WT vs SRSF1KO') +
              ylab("-log10(padj)")  

ggsave(paste(fig_dir, "/srsf1KO_DE_circ.svg", sep=""), plot=srsf1p_cir, width=3.5, height = 3.5)
ggsave(paste(fig_dir, "/srsf1KO_DE_common.svg", sep=""), plot=srsf1p_common, width=3.5, height = 3.5)

write.csv(results_SRSF1K0,paste(list_dir,"/results_SRSF1K0_DE.csv", sep=""),row.names = FALSE)

#----------------------srsf3KO--------------------#
files <- list.files(path = "../results/SRP212166_salmon_2024-12-06/", 
                    pattern = "*.sf", full.names = TRUE, include.dirs = TRUE, no.. = TRUE, recursive = TRUE)
names(files) <- str_extract(files, "SRR\\d+")

txi_SRP212166 <- tximport(files, type = "salmon", txOut = TRUE, countsFromAbundance = "no", ignoreTxVersion = TRUE)

meta_data_SRP212166 <- read_csv("../meta/SRP212166_SraRunTable.txt") %>% 
  dplyr::select(!!!c("Run", "Genotype")) %>% 
  dplyr::mutate(Genotype = ifelse(grepl("KO", Genotype), "SRSF3KO", "WT"),
                Genotype = fct(Genotype, levels=c("WT", "SRSF3KO"))) %>%
  filter(Run %in% colnames(txi_SRP212166$abundance))

counts <- txi_SRP212166$counts
mode(counts) <- "integer"

srsf3_dds <- DESeqDataSetFromMatrix(countData = counts,
                                    colData = meta_data_SRP212166,
                                    design= ~ Genotype)
rowData(srsf3_dds) <- annotation[rownames(srsf3_dds),]
keep <- rownames(counts) %in% expressed_t$transcript_id
srsf3_dds <- srsf3_dds[keep, ]
srsf3_dds <- DESeq(srsf3_dds, fitType = "local")
srsf3_res <- DESeq2::results(srsf3_dds, c("Genotype", "SRSF3KO", "WT"), independentFiltering = TRUE, alpha = fdr_cutoff)
srsf3_res <- lfcShrink(srsf3_dds, coef="Genotype_SRSF3KO_vs_WT", res = srsf3_res,type = "apeglm")
results_SRSF3KO <- srsf3_res %>% as.data.frame() %>%
  filter(!is.na(padj)) %>%
  rownames_to_column(var="transcript_id") %>%
  dplyr::select(c("transcript_id", "log2FoldChange", "padj")) %>%
  dplyr::rename(padj_srsf3ko = "padj", log2FC_srsf3ko = "log2FoldChange") %>% left_join(., annotation) %>%
  mutate(diffexp= ifelse(.$log2FC_srsf3ko > 1.5 & .$padj_srsf3ko < 0.05, "UP", ifelse(.$log2FC_srsf3ko < -1.5 & .$padj_srsf3ko < 0.05,"DOWN","No")),
         cir=case_when(!is.na((match(.$transcript_id, circadian_transcripts$transcript_id))) ~ "Yes",
                       .default = "No"))

interesting<- c("Rorc","Fus","Cdk1a","Fasn")

srsf3p_common <- ggplot(data=results_SRSF3KO, aes(x=log2FC_srsf3ko, y=-log10(padj_srsf3ko), col=class)) + 
              geom_point(size=0.4) + 
              geom_vline(xintercept=c(-1.5, 1.5), col="black", linetype = "dashed") +
              geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed") + 
              theme_classic() + labs(color="Annotation") +
              theme(text=element_text(size= 10),strip.text.x = element_blank(), legend.position="bottom",
                    legend.box.spacing = unit(0, "pt")) +
              scale_colour_manual(values = mycolors_class) +
              xlab('log2(Expression Fold Change) WT vs SRSF3KO') +
              ylab("-log10(padj)") 

srsf3p_cir <- ggplot(data=results_SRSF3KO, aes(x=log2FC_srsf3ko, y=-log10(padj_srsf3ko), col=cir)) + 
              geom_point(size=0.4) + theme_minimal()+
              geom_vline(xintercept=c(-1.5, 1.5), col="black", linetype = "dashed") +
              geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed") + 
              theme_classic() + 
              geom_text_repel(data=results_SRSF3KO %>% filter(cir == "Yes")%>% filter(diffexp != "No") %>% filter(mgi_symbol %in% interesting),
                              aes(y=-log10(padj_srsf3ko), x= log2FC_srsf3ko,label= mgi_symbol), col= "black", max.overlaps=20) +
              theme(text=element_text(size= 10),strip.text.x = element_blank(),legend.position = "bottom", 
                    legend.box.spacing = unit(0, "pt")) +labs(color="Circadian") +
              xlab('log2(Expression Fold Change) WT vs SRSF3KO') +
              ylab("-log10(padj)") 

ggsave(paste(fig_dir, "/srsf3KO_DE_circ.svg", sep=""), plot=srsf3p_cir, width=3.5, height = 3.5)
ggsave(paste(fig_dir, "/srsf3KO_DE_common.svg", sep=""), plot=srsf3p_common, width=3.5, height = 3.5)

write.csv(results_SRSF3KO,paste(list_dir,"/results_SRSF3K0_DE.csv", sep=""), row.names = FALSE)

#srsf3_volcano<- srsf3p_A + srsf3p_NA + srsf3p_NA
#srsfs_volcano<-srsf1_volcano/srsf3_volcano
#srsfs_volcano<- plot_grid(srsf1_volcano, srsf3_volcano, rel_widths = c(1,1), rel_heights = c(1,1), nrow=2, ncol = 1)
#ggsave("../extra/srsfs_volcano.png", plot= srsfs_volcano, width=10, height = 4)


####
---------------------------------
  results_SRSF <- inner_join(results_SRSF1K0, results_SRSF3KO) %>%
  left_join(annotation) %>%
  left_join(mart)
ggplot(results_SRSF) + 
  geom_point(aes(x=log2FC_srsf1ko, y=log2FC_srsf3ko, alpha=abs(log2FC_srsf1ko)>1 | abs(log2FC_srsf3ko)>1), size=0.5) + 
  theme_classic(base_size = base_size)



#read per file here
results<- read.csv("../data/DE.csv")
Persig<-filter(results,diffexp=='UP'| diffexp=='DOWN')
srsf1sig<- filter(results_SRSF1K0,diffexp=='UP'| diffexp=='DOWN')
srsf3sig<- filter(results_SRSF3KO,diffexp=='UP'| diffexp=='DOWN')
persig<- filter(results,diffexp=='UP'| diffexp=='DOWN')
common_srsfs<- intersect(srsf3sig$transcript_id, srsf1sig$transcript_id)
bz<-srsf1sig %>% filter(transcript_id %in% common_srsfs)

common_srsf1<- intersect(Persig$transcript_id, srsf1sig$transcript_id)
common_srsf3<- intersect(Persig$transcript_id, srsf3sig$transcript_id)
common_per_srsf<- intersect(common_srsf1,common_srsf3)
venn.diagram(x = list(Persig$transcript_id,srsf1sig$transcript_id, srsf3sig$transcript_id),category.names = c("Per" , "srsf1", "srsf3"),filename = '#14_venn_diagramm_wtthtresh.png',
             output=TRUE,
             fill=c("#FFBEDE", "#BEE5FF", "green"),
             cex = 1.5,
)



#--------------------------------------#
#srsf1p_cir_A <- ggplot(data=results_SRSF1K0 %>% filter(cir=="Yes" & class == "Annotated"), aes(x=log2FC_srsf1ko, y=-log10(padj_srsf1ko), col=diffexp)) + geom_point() + theme_minimal()+ geom_vline(xintercept=c(-1, 1), col="red") +geom_hline(yintercept=-log10(0.05), col="red" )+ scale_colour_manual(values = mycolors_updown)+ theme(legend.position="none")+ylim(c(0,150))  
#srsf1p_cir_NA <- ggplot(data=results_SRSF1K0 %>% filter(cir=="Yes" & class == "Novel-Annotated"), aes(x=log2FC_srsf1ko, y=-log10(padj_srsf1ko), col=diffexp)) + geom_point() + theme_minimal()+ geom_vline(xintercept=c(-1, 1), col="red") +geom_hline(yintercept=-log10(0.05), col="red" )+ scale_colour_manual(values = mycolors_updown)+ theme(legend.position="none")+ylim(c(0,150))  
#srsf1p_cir_NN <- ggplot(data=results_SRSF1K0 %>% filter(cir=="Yes" & class == "Novel-Novel"), aes(x=log2FC_srsf1ko, y=-log10(padj_srsf1ko), col=diffexp)) + geom_point() + theme_minimal()+ geom_vline(xintercept=c(-1, 1), col="red") +geom_hline(yintercept=-log10(0.05), col="red" )+ scale_colour_manual(values = mycolors_updown)+ylim(c(0,150))  
#
#srsf1p_A <- ggplot(data=results_SRSF1K0 %>% filter(class=="Annotated"), 
#                   aes(x=log2FC_srsf1ko, y=-log10(padj_srsf1ko), col=diffexp)) +
#  geom_point(size=0.8) + theme_minimal()+ 
#  geom_vline(xintercept=c(-1, 1), col="black") + 
#  geom_hline(yintercept=-log10(0.05), col="black")+
#  scale_colour_manual(values = mycolors_updown)+ theme(legend.position="none") + 
#  geom_text(size=3,aes(label=ifelse(abs(log2FC_srsf1ko)> 5 & -log10(padj_srsf1ko) > 50, mgi_symbol, ""),hjust=0,vjust=0), 
#            color="black")
#srsf1p_NN <- ggplot(data=results_SRSF1K0 %>% filter(class=="Novel-Novel"), 
#                    aes(x=log2FC_srsf1ko, y=-log10(padj_srsf1ko), col=diffexp)) + geom_point(size=0.8) + theme_minimal()+ geom_vline(xintercept=c(-1, 1), col="black") +geom_hline(yintercept=-log10(0.05), col="black")+ scale_colour_manual(values = mycolors_updown)+ylim(c(0,150))  + geom_text(size=3,aes(label=ifelse(abs(log2FC_srsf1ko)> 3 & -log10(padj_srsf1ko) > 20, mgi_symbol, ""),hjust=0,vjust=0), color="black") + theme(axis.title.y = element_blank())
#srsf1p_NA <- ggplot(data=results_SRSF1K0 %>% filter(class=="Novel-Annotated"), aes(x=log2FC_srsf1ko, y=-log10(padj_srsf1ko), col=diffexp)) + geom_point(size=0.8) + theme_minimal()+ geom_vline(xintercept=c(-1, 1), col="black") +geom_hline(yintercept=-log10(0.05), col="black")+ scale_colour_manual(values = mycolors_updown)+ylim(c(0,150))+ theme(legend.position="none")+ geom_text(size=3,aes(label=ifelse(abs(log2FC_srsf1ko)> 4 & -log10(padj_srsf1ko) > 25, mgi_symbol, ""),hjust=0,vjust=0), color="black")+ theme(axis.title.y = element_blank())
##---------------------------------SRSF1 DTU--------------------------------------------------##

#cts <- data.frame(gene_id = rowData(srsf1_dds)$gene_id,
#                  feature_id = rowData(srsf1_dds)$transcript_id,
#                 counts(srsf1_dds), check.names = FALSE)
#exp_design_SRP252896 <- meta_data_SRP252896 %>% 
#  as.data.frame %>%
##  dplyr::rename(sample_id = "Run", group = "Genotype")
#d <- DRIMSeq::dmDSdata(counts = cts %>% filter(!is.na(feature_id)), samples = exp_design_SRP252896)
#d <- DRIMSeq::dmFilter(d,
#                       min_samps_feature_expr = 4,
#                       min_feature_expr=count_threshold,
#                       min_feature_prop=0.1,
#                       min_gene_expr=2*count_threshold)
#dxd <- DEXSeqDataSet(countData = as.matrix(counts(d)[,-c(1,2)]),
#                     sampleData = DRIMSeq::samples(d),
#                     design =~sample_id + exon + group:exon,
#                     featureID = counts(d)$feature_id,
#                     groupID = counts(d)$gene_id %>% sub(':','loc',.))
#
#BPPARAM<-  BiocParallel::MulticoreParam(6)
##dxd <- estimateSizeFactors(dxd)
#dxd <- estimateDispersions(dxd, quiet=TRUE, fitType = 'local', BPPARAM=BPPARAM)
#dxd <- testForDEU(dxd, reducedModel=~sample_id + exon, BPPARAM=BPPARAM)
#dxr <- DEXSeqResults(dxd, independentFiltering=TRUE)
#qval <- DEXSeq::perGeneQValue(dxr)
#pScreen <- qval
#annotation$gene_id<-  sub(':','loc',annotation$gene_id) #why not above?
#tx2gene<- annotation[c("transcript_id", "gene_id")]

#pConfirmation <- matrix(dxr$pvalue, ncol = 1)
#dimnames(pConfirmation) <- list(dxr$featureID, "transcript")
#stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
 #                     pScreenAdjusted=TRUE, tx2gene= tx2gene)
#stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.1, allowNA=TRUE)
#suppressWarnings({
#  dtu_padj <- getAdjustedPValues(stageRObj, order=FALSE,
 #                                onlySignificantGenes=TRUE)
#})
##drim.padj <-merge(tx2gene, dtu_padj, by.x = "transcript_id", by.y = "txID") %>%
#  arrange(gene)
#write_rds(dtu_padj, "../results/dtu_padj_srsf1.rds", compress = "bz2")
#dtu_padj %>% dplyr::select(txID, gene, geneID, transcript) %>%
#  filter(gene<fdr_cutoff) %>%
#  write_excel_csv(file = "DTU_genes_srsf1.csv")

#DTU_genes_srsf1Ko<-dtu_padj %>% dplyr::select(txID, gene, geneID, transcript) %>%
#  filter(gene<fdr_cutoff) %>% mutate(cir=case_when(!is.na((match(.$txID, circadian_transcripts$transcript_id))) ~ "Yes", 
#                                                   .default = "No")) 

#%>%  left_join(., annotation) not workingg even after putting by argument?????

#DTU_genes_srsf1Ko$genesymbol<- annotation$mgi_symbol[match(DTU_genes_srsf1Ko$txID, annotation$transcript_id)]
#DTU_genes_srsf1Ko$gene_id<- annotation$gene_id [match(DTU_genes_srsf1Ko$txID, annotation$transcript_id)]

#design_full <- model.matrix(~group, data=DRIMSeq::samples(d))
#set.seed(1)
#system.time({
#  d <- dmPrecision(d, design=design_full)
#  d <- dmFit(d, design=design_full)
#  d <- dmTest(d, coef="groupSRSF1KO")
#})

#srsf1_DTU_summary <- DTU_genes_srsf1Ko %>%  
#  group_by(cir) %>%  
#  summarize(gene = n_distinct(gene_id), tx = n_distinct(txID))

#dex.norm <- cbind(as.data.frame(stringr::str_split_fixed(rownames(counts(dxd)), ":", 2)), as.data.frame(counts(dxd, normalized = TRUE))[,1:4])
#colnames(dex.norm) <- c("groupID", "featureID", as.character(colData(dxd)$sample_id)[1:4])
#exp_dbp_srsf1<-plotExpression(dex.norm, "ENSMUSG00000059824",exp_design_SRP252896 , isProportion = FALSE)
#proportions_dbp_srsf1<- plotProportions(d,"ENSMUSG00000059824", "group")
#srsf1_transcriptome<- plot_grid(plot_grid(p_common, p_A,p_NA,p_NN,rel_widths = c(1.5,1,1,1), nrow = 1, ncol=5),plot_grid(p_cir,NULL, p_cir_A, p_cir_NA, p_cir_NN, nrow=1, ncol=4, rel_widths = c(1,0.5,1,1,1)), plot_grid(exp_dbp_srsf1, proportions_dbp_srsf1, rel_widths = c(1,1), ncol=2), nrow = 3)


###-------------------------------DTU srsf3Ko-------------------------------#
#cts <- data.frame(gene_id = rowData(srsf3_dds)$gene_id,
#                  feature_id = rowData(srsf3_dds)$transcript_id,
#                  counts(srsf3_dds), check.names = FALSE)
#exp_design_SRP212166 <- meta_data_SRP212166 %>% 
#  as.data.frame %>%
#  dplyr::rename(sample_id = "Run", group = "Genotype")
#d <- DRIMSeq::dmDSdata(counts = cts %>% filter(!is.na(feature_id)), samples = exp_design_SRP212166)
#d <- DRIMSeq::dmFilter(d,
#                       min_samps_feature_expr = 4,
#                       min_feature_expr=count_threshold,
#                       min_feature_prop=0.1,
#                       min_gene_expr=2*count_threshold)
#dxd <- DEXSeqDataSet(countData = as.matrix(counts(d)[,-c(1,2)]),
#                     sampleData = DRIMSeq::samples(d),
#                     design =~sample_id + exon + group:exon,
#                     featureID = counts(d)$feature_id,
#                     groupID = counts(d)$gene_id %>% sub(':','loc',.))

#BPPARAM<-  BiocParallel::MulticoreParam(6)
#dxd <- estimateSizeFactors(dxd)
#dxd <- estimateDispersions(dxd, quiet=TRUE, fitType = 'local', BPPARAM=BPPARAM)
#dxd <- testForDEU(dxd, reducedModel=~sample_id + exon, BPPARAM=BPPARAM)
#dxr <- DEXSeqResults(dxd, independentFiltering=TRUE)
#qval <- DEXSeq::perGeneQValue(dxr)
#pScreen <- qval
#annotation$gene_id<-  sub(':','loc',annotation$gene_id) #why not above?
#tx2gene<- annotation[c("transcript_id", "gene_id")]#

#pConfirmation <- matrix(dxr$pvalue, ncol = 1)
#dimnames(pConfirmation) <- list(dxr$featureID, "transcript")
#stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
#                      pScreenAdjusted=TRUE, tx2gene= tx2gene)
#stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.1, allowNA=TRUE)
#suppressWarnings({
#  dtu_padj <- getAdjustedPValues(stageRObj, order=FALSE,
#                                 onlySignificantGenes=TRUE)
#})
#drim.padj <-merge(tx2gene, dtu_padj, by.x = "transcript_id", by.y = "txID") %>%
#  arrange(gene)
#write_rds(dtu_padj, "../results/dtu_padj_srsf3.rds", compress = "bz2")
#dtu_padj %>% dplyr::select(txID, gene, geneID, transcript) %>%
#  filter(gene<fdr_cutoff) %>%
#  write_excel_csv(file = "DTU_genes_srsf3.csv")#
#
#DTU_genes_srsf3Ko<-dtu_padj %>% dplyr::select(txID, gene, geneID, transcript) %>%
#  filter(gene<fdr_cutoff) %>% mutate(cir=case_when(!is.na((match(.$txID, circadian_transcripts$transcript_id))) ~ "Yes", 
#                                                   .default = "No")) 
#%>%  left_join(., annotation,by.x="txID", by.y="transcript_id")


#design_full <- model.matrix(~group, data=DRIMSeq::samples(d))
#set.seed(1)
#system.time({
#  d <- dmPrecision(d, design=design_full)
#  d <- dmFit(d, design=design_full)
#  d <- dmTest(d, coef="group")
#})

##srsf3p_cir_A <- ggplot(data=results_SRSF3KO %>% filter(cir=="Yes" & class == "Annotated"), aes(x=log2FC_SRSF3KO, y=-log10(padj_SRSF3KO), col=diffexp)) + geom_point() + theme_minimal()+ geom_vline(xintercept=c(-1, 1), col="red") +geom_hline(yintercept=-log10(0.05), col="red" )+ scale_colour_manual(values = mycolors_updown)+ theme(legend.position="none")+ylim(c(0,150))  
#srsf3p_cir_NA <- ggplot(data=results_SRSF3KO %>% filter(cir=="Yes" & class == "Novel-Annotated"), aes(x=log2FC_SRSF3KO, y=-log10(padj_SRSF3KO), col=diffexp)) + geom_point() + theme_minimal()+ geom_vline(xintercept=c(-1, 1), col="red") +geom_hline(yintercept=-log10(0.05), col="red" )+ scale_colour_manual(values = mycolors_updown)+ theme(legend.position="none")+ylim(c(0,150))  
#srsf3p_cir_NN <- ggplot(data=results_SRSF3KO %>% filter(cir=="Yes" & class == "Novel-Novel"), aes(x=log2FC_SRSF3KO, y=-log10(padj_SRSF3KO), col=diffexp)) + geom_point() + theme_minimal()+ geom_vline(xintercept=c(-1, 1), col="red") +geom_hline(yintercept=-log10(0.05), col="red" )+ scale_colour_manual(values = mycolors_updown)+ylim(c(0,150))  



###
#--------------------------functions--------------------------------#
#https://ycl6.gitbook.io/guide-to-rna-seq-analysis/differential-expression-analysis/differential-transcript-usage/dtu-using-dexseq

#plotExpression <- function(expData = NULL, geneID = NULL, samps = NULL, isProportion = FALSE) {
#  colnames(expData)[1:2] = c("gid","tid")
#  sub = subset(expData, gid == geneID)
#  sub = reshape2::melt(sub, id = c("gid", "tid"))
#  sub = merge(samps, sub, by.x = "sample_id", by.y = "variable")
#  if(!isProportion) {
#    sub$value = log(sub$value)
#  }
#  clrs = c("dodgerblue3", "maroon2",  "forestgreen", "darkorange1", "blueviolet", "firebrick2",
#           "deepskyblue", "orchid2", "chartreuse3", "gold", "slateblue1", "tomato" , "blue", "magenta", "green3",
#           "yellow", "purple3", "red" ,"darkslategray1", "lightpink1", "lightgreen", "khaki1", "plum3", "salmon")
#  p = ggplot(sub, aes(tid, value, color = group, fill = group)) +
#    geom_boxplot(alpha = 0.4, outlier.shape = NA, width = 0.8, lwd = 0.5) +
#    stat_summary(fun = mean, geom = "point", color = "black", shape = 5, size = 3, position=position_dodge(width = 0.8)) +
#    scale_color_manual(values = clrs) + scale_fill_manual(values = clrs) +
#    geom_quasirandom(color = "black", size = 1, dodge.width = 0.8) + theme_bw() +
#    ggtitle(geneID) + xlab("Transcripts")
  
#  if(!isProportion) {
#    p = p + ylab("log(Expression)")
#  } else {
#    p = p + ylab("Proportions")
#  }
#  p
#}
#no.na <- function(x) ifelse(is.na(x), 1, x)