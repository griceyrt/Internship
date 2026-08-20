setwd("C:/Users/kmamgain/prettyCoolPores_new/scripts")

library(tximport)
library(stringr)
library(tidyr)
library(tidyverse)
library(DESeq2)
library(stageR)
library(DEXSeq)
library(ggplot2)
library(DRIMSeq)
library(dplyr)
library(cowplot)
library(stageR)
library(biomaRt)
library(hrbrthemes)
library(ggpubr)
library(clusterProfiler)
library(AnnotationHub)
library(grid)
library(org.Mm.eg.db)
library(GenomicFeatures)
library(clusterProfiler)
library(AnnotationHub)
library(ggrepel)
library(ggh4x)
library(ReactomePA)
library(cowplot)

fig_dir <- paste(getwd(), "/../Manuscript_figures_",Sys.Date(), sep="") #Path of directory where we want to save all the extra figures
theme_legend<- theme(legend.key.size = unit(0.1, 'cm'), #change legend key size
                     legend.key.height = unit(0.5, 'cm'), #change legend key height
                     legend.key.width = unit(0.1, 'cm'), #change legend key width
                     legend.title = element_blank())
flair_path<- "../results/FLAIR_2024-09-11"
biomart_path<-"../data/biomart_annotation_2024-11-05.rds"
salmon_path<-"../results/salmon_2024-09-25/"
suppa_path<- "../results/suppa_2024-11-06"
list_dir<- "../results/lists"

fdr_cutoff <- 0.05
retrieve_from_biomart <- FALSE
tx_names_with_version <- TRUE
mycolors_updown <- c("DOWN" = "purple", "NO"="grey", "UP" ="#FFC94C")
mycolors_productivtiy <- c("Productive" = "black", "Non-productive"="#ffbf00")
mycolors_class <- c("Novel-Annotated"= "#fc8d62","Novel-Novel"= "#8da0cb","Annotated"= "#66c2a5")
x_theme<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
expressed_list<- read.csv("../results/expr_tx_list.csv")

#load gtf and expression files
txdb <- makeTxDbFromGFF(paste(flair_path,"/transcriptome_productivity.gtf", sep=""), format = "gtf")

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

files <- list.files(path = salmon_path, pattern = ".sf", 
                    full.names = TRUE, include.dirs = FALSE, no.. = TRUE, recursive = TRUE)
names(files) <- str_extract(files, "([CKTO]\\d+_\\d)(\\.sf)$", group = 1)
txi.salmon <- tximport(files[c(7,8,13,14)], type = "salmon", txOut = TRUE, countsFromAbundance = "no", ignoreTxVersion = TRUE, 
                       importer = function(x) read.delim(x))
exp_design <- data.frame(batch = factor(str_extract(colnames(txi.salmon$abundance), "(?<=_)\\d")),
                         sample_id= colnames(txi.salmon$abundance),
                         time = as.numeric(str_extract(colnames(txi.salmon$abundance), "\\d+")),
                         genotype= ifelse(grepl("^T", colnames(txi.salmon$abundance)), "WT", "KO")) %>% mutate(genotype = fct(genotype, levels=c("WT", "KO")))
rownames(exp_design) <- exp_design$sample_id
counts <- magrittr::set_colnames(txi.salmon$counts, exp_design$sample_id)
if (tx_names_with_version){
  rownames(counts) <- str_extract(rownames(counts), "([A-Za-z0-9-]+)")
}
mode(counts) <- "integer"
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = exp_design,
                              design = ~ genotype)
keep <- rownames(counts) %in% expressed_list$transcript_id #>= 5 ncol(counts(dds))/2
dds <- dds[keep ]
dds <- DESeq(dds, fitType = "local")
res <- DESeq2::results(dds,c("genotype","KO","WT"),independentFiltering = TRUE, alpha = fdr_cutoff, tidy = TRUE)

#res <- lfcShrink(dds, coef="genotype_KO_vs_WT", res = res,type = "apeglm")    #what is the significance, is it okat to plot without?
res$groupID <- rowData(dds)$gene_id
res <- cbind(res, coef(dds))
results <- res %>% as.data.frame() %>%
  filter(!is.na(padj)) %>%
  rownames_to_column(var="transcript_id") %>%
  dplyr::select(c("transcript_id", "log2FoldChange", "padj")) %>%
  dplyr::rename(padj_per_ko = "padj", log2FC_perko = "log2FoldChange") %>% left_join(., annotation) %>%
  mutate(diffexp= ifelse(.$log2FC_perko > 1.5 & .$padj_per_ko < 0.05, "UP", ifelse(.$log2FC_perko < -1.5 & .$padj_per_ko < 0.05,"DOWN","No")))

productivity <- read.table(paste(flair_path,"/productivity.txt", sep="")) %>% mutate(transcript_id= .$V4, productivity_flair= .$V13) %>% 
    left_join(., annotation) %>% 
  mutate(productivity_flair = fct_relevel(productivity_flair, c("PRO","PTC","NGO","NST")),
         productivity_flair = fct_recode(productivity_flair, 
                                         !!!c("productive" = "PRO",
                                              "premature termination codon" = "PTC",
                                              "no start codon" = "NGO",
                                              "has start but no stop codon" = "NST")))%>% dplyr::select(c("transcript_id", "gene_id","productivity_flair", "class")) %>% 
  filter(transcript_id %in% expressed_list$transcript_id)


Perko_DE_summary <- results %>% filter(diffexp != "No") %>% left_join(.,productivity) %>% mutate(productivity= ifelse(productivity_flair =="productive", "Productive", "Non-productive"))%>% group_by(class,diffexp,productivity) %>% 
  summarize(gene = n_distinct(gene_id), 
            tx = n_distinct(transcript_id)) %>% mutate(percent = tx/sum(tx),
                                                       ymax=cumsum(percent), 
                                                       ymin=c(0, head(ymax, n=-1)),
                                                       labelPosition = (ymax + ymin) / 2,
                                                       label = paste0(tx)) %>% arrange(desc(productivity)) %>% 
                                                        mutate(text_y = 0.5*(ymax + ymin)) #added this because geom_label repel below shows wrong positions

volcanos_common_theme<-  theme(legend.position="none", axis.title.x = element_blank(), axis.text.x =  element_blank())
Persig<-filter(results,diffexp=='UP'| diffexp=='DOWN')

#save significant perKO an results
write.csv(results, "../results/DE_all.csv", row.names = FALSE)
write.csv(Persig, "../results/perKO_sig.csv")
write.csv(Perko_DE_summary, "../results/Perko_DE_summary.csv")


#not only clock but the ones that we want to label in the bigvolcano
#keep the one you want to have names in zoom as the first in this list or change it manually later
clock_genes <- c("Cyp2a5","Arntl", "Per1", "Dbp", "Rorc", "Clock", "Hlf", "Bhlhe40", "Nr1d2", "Cry1", "Cry2", "Nr1d1", "Per2", "Noct", "Bhlhe40", "Arntl2", "Nampt", "Ciart")

df<- results %>% 
  #filter(diffexp != "No") %>% 
  group_by(diffexp) %>% 
  summarise(n= n())
DE_clock<- c("Mcm6","Camk1d","Hnf4a", "Dbp", "Nr1d2", "Nr1d1", "Acly", "Hmga1b", "Comt", "Ppard", "Egr1")


results <- results %>%
  mutate(col_volcano = ifelse(diffexp == 'No', "grey80",
                              ifelse(diffexp != "No" & class == 'Annotated', '#66c2a5', 
                                     ifelse(diffexp != "No" & class == 'Novel-Novel', '#8da0cb', 
                                            ifelse(diffexp != "No" & class == 'Novel-Annotated', '#fc8d62', "black")))))

##
names<-c("Camk1d","Ciart","Acnat2","Mup20","Aldh2","Sult2a8","Cyp2a5","Nr1d2","Acy1","Nr1d1","Acly")


bigvolcano <- ggplot(data=results, aes(y=-log10(padj_per_ko), x= log2FC_perko, color=col_volcano)) + 
  geom_point(size=0.8) +
  scale_color_manual(values= c('#66c2a5', '#8da0cb', '#fc8d62', "grey80")) +
  facet_wrap('class', nrow=3, scale="free_y") + #, scale="free_y"
  geom_vline(xintercept=c(-1.5, 1.5), col="black", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed") + 
  theme_classic() + 
  theme(text=element_text(size= 10),strip.text.x = element_blank(), legend.title = element_blank(), legend.position = "none",legend.box.spacing = unit(0, "pt")) +
  xlab('log2(Expression Fold Change) WT vs PerKO') +
  ylab('-log10(padj)') + 
  geom_text_repel(data=results %>% filter(mgi_symbol %in% names) %>% filter(abs(log2FC_perko) > 1.5 & padj_per_ko < 0.05), aes(y=-log10(padj_per_ko), x= log2FC_perko, label= mgi_symbol), col= "black", max.overlaps=20)

ggsave(paste(fig_dir,"/4_volcanos_free_Y.svg",sep=""), plot= bigvolcano, width=5, height = 6)

col<- c("#66c2a5","#8da0cb","#fc8d62") #facet wrap doesnt allow separate colouring of the titles so used tis packages for the same.
strip <- strip_themed(background_y = elem_list_rect(fill = col), text_y = elem_list_text(color= col, size=c(0.5,0.5,0.5)))
p_upregulated<-ggplot(Perko_DE_summary[Perko_DE_summary$diffexp=="UP",], aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill= productivity)) + scale_fill_manual(values= mycolors_productivtiy)+
  geom_rect(color= "black") +
  coord_polar(theta="y") +
  geom_text(x=5, aes(y=text_y, label=label),fill="white", colour="black") + # x here controls label position (inner / outer)
  xlim(c(2, 6))+ # remove hash if to make it a donut
  theme_void() + theme(plot.title =  element_text(hjust = 0.5),text=element_text(size=5),legend.title = element_blank(), legend.position = "top", legend.box = "horizontal", legend.key.size = unit(0.3, "cm"),legend.box.spacing = unit(15,"pt"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  facet_grid2(class~., strip = strip, switch= "y") 
p_downregulated<-ggplot(Perko_DE_summary[Perko_DE_summary$diffexp=="DOWN",], aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill= productivity)) + scale_fill_manual(values= mycolors_productivtiy)+
  geom_rect(color= "black") +
  geom_label(x= 5, aes(y=labelPosition, label=label),fill="white", colour="black") + # x here controls label position (inner / outer)
  coord_polar(theta="y") +
  xlim(c(2, 6)) + #remove hash to make it donut
  theme_void() + theme(plot.title =  element_text(hjust = 0.5), legend.position = "none",legend.justification="centre", legend.title = element_blank(),
                       legend.box.spacing = unit(05, "pt"), legend.box = "horizontal", legend.key.size = unit(0.3, "cm"), text=element_text(size=5))+ 
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  facet_grid2(class~., strip = strip, switch="y") 

res %>% 
  as.data.frame() %>%
  rownames_to_column(var="transcript_id") %>%
  dplyr::select(c("transcript_id", "log2FoldChange", "padj")) %>%
  dplyr::rename(padj_per_ko = "padj", log2FC_perko = "log2FoldChange") %>% 
  left_join(., annotation) %>%
  mutate(diffexp= ifelse(.$log2FC_perko > 1.5 & .$padj_per_ko < 0.05, "UP", ifelse(.$log2FC_perko < -1.5 & .$padj_per_ko < 0.05,"DOWN","No"))) %>%
  dplyr::select(c("transcript_id", "log2FC_perko", "padj_per_ko", "gene_id", "mgi_symbol", "transcript_biotype", "class", "diffexp")) %>%
  write.table(., "../lists/Fig4A.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

#---------------------------------------#
#-------------- Lies -------------------# 
#---------------------------------------#
# Perko_DE_summary <- Perko_DE_summary %>%
#   mutate(class = str_replace(class,"-", "-\n"))

up<-ggplot(filter(Perko_DE_summary, diffexp == "UP"), aes(x=class, y=tx, group=productivity, fill=class, alpha=productivity))+ 
  geom_bar(stat='identity', width = 0.7, col="black") + 
  theme_classic()+ ylim(c(0,130)) + theme(axis.title.y = element_text(size=7), axis.title.x = element_blank(), axis.text.x = element_blank(), legend.key.size = unit(0.3, "cm"),
                                          legend.position="top", legend.direction="vertical", legend.box="vertical", legend.title = element_blank()) + 
  scale_fill_manual(values=c("#66c2a5","#fc8d62","#8da0cb")) + scale_alpha_manual(values = c(0.3,1)) + ylab("Up-regulated") + guides(fill="none", alpha = "none")
down<-ggplot(filter(Perko_DE_summary, diffexp == "DOWN"), aes(x=class, y=tx, group=productivity, fill=class, alpha=productivity)) + 
  geom_bar(stat='identity', width = 0.7, col="black") + theme_classic() + scale_x_discrete(position = "top") + 
  scale_y_reverse(limits=c(130,0)) + theme(axis.title.y = element_text(size=7), plot.margin=unit(c(0,0,0,0), 'cm'), axis.title.x = element_blank(), legend.position="bottom", 
  legend.direction="vertical", legend.box="vertical", axis.text.x = element_blank(), legend.title = element_blank(), legend.key.size = unit(0.3, "cm"))+ 
  scale_fill_manual(values=c("#66c2a5","#fc8d62","#8da0cb")) + scale_alpha_manual(values = c(0.3,1)) + ylab("Down-regulated") + guides(fill = "none")

nb_tx_barplot <- plot_grid(up,down,ncol=1, align = "v", axis = "lr", rel_heights = c(1,1))

####GO for persig### code below
hub <- AnnotationHub()
q <- query(hub, "org.Mm.eg.db")
id <- q$ah_id[length(q)]
org.Mm.eg.db <- hub[[id]]

###Reactome enrichment
#convert to entrezID
genes_forkegg <- bitr(Persig$gene_id,
                      fromType="ENSEMBL",
                      toType = "ENTREZID",
                      OrgDb = org.Mm.eg.db) 
genes_kegguniverse <- bitr(results$gene_id,
                           fromType="ENSEMBL",
                           toType = "ENTREZID",
                           OrgDb = org.Mm.eg.db)
reactome<-enrichPathway(
  genes_forkegg$ENTREZID,
  organism = "mouse",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  genes_kegguniverse$ENTREZID,
  qvalueCutoff = 0.2,
  readable = TRUE
)
reactome_result<- reactome@result %>% filter(p.adjust< 0.05)
#if you want in GO visualisation style
reactome_persig<-ggplot(reactome_result[1:10,]) + theme(panel.border = element_rect(color="black", fill="transparent"), 
                                                        panel.background = element_rect(fill="transparent"),  text = element_text(size=8),
                                                        legend.text = element_text(size=6),legend.key.height = unit(0.1, "cm") )+
                 geom_point(aes(x=-log(p.adjust), y=reorder(Description,Count ), size=Count, color= -log(p.adjust))) + ylab("")
reactome_cnet<- cnetplot(reactome, organism = "mouse", showCategory=20, category_node  = 1, #cnetplot
                         gene_node  = 0.5,  
                         layout="fr",    #different orientations/shapes of the cnet
                         category_label  = 0.8, #label size for category
                         gene_label  = 0.8, #label size for genes
                         hilight.params= list(category= "Biological oxidations",alpha_hilight = 1.5, alpha_no_hilight = 1),
                         cex_label_category = 0.7,
                         cex_label_gene = 0.5) +
  theme(legend.position = "none")

write.table(reactome_result, "../lists/Reactome_results.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

#---------------------------------------#
#---------- Distance to Koike ----------# 
#---------------------------------------#

# Import data
DE_res <-read.csv("../results/DE_all.csv")
productivity <- read.table(paste(flair_path,"/productivity.txt", sep="")) %>%
  setNames(c('chr', 'start', 'end', 'transcript_id', 'score', 'strand', 'CDS_start', 'CDS_end', 'itemRgb', 'nb_exons', 'size_exons', 'start_exons', 'productivity')) %>%
  mutate(chr = paste('chr', chr, sep="")) %>%
  filter(transcript_id %in% DE_res$transcript_id)

binding_factors <- c('Per1', 'Per2', 'Cry1', 'Cry2', 'Bmal1', 'Clock', 'Npas2')
Koike_df <- data.frame()
for (factor in binding_factors){
  table <- read.table(paste("../data/Koike_", factor, "_mm39.bed", sep="")) %>% 
    setNames(c('chr', 'start', 'end')) %>% 
    mutate(mid_position=(as.numeric(start)+as.numeric(end))/2, 
           binding_protein=factor)
  Koike_df <- rbind(Koike_df, table)
}
BS_df2 <- data.frame()

for (i in 1:nrow(productivity)){
  list_same_chr <- filter(Koike_df, chr == productivity[i,1]) # filter positions on the same chr
  df <- data.frame(data.frame(matrix(NA, nrow = 1)))
  for (factor in binding_factors){ # filter upstream positions for the given factor
    x <- filter(list_same_chr, binding_protein == factor) 
    if (productivity[i,6] == "+") {
      x <- filter(x, mid_position < productivity[i,2])
    } else {
      x <- filter(x, mid_position > productivity[i,2]) # get upstream for negative strand genes
    } 
    closest_BS2 <- x$mid_position[which.min(abs(x$mid_position - productivity[i,2]))] # get the closest binding site to the given TSS
    closest_BS2 <- ifelse(nrow(list_same_chr) == 0, "NA", closest_BS2)  
    df <- cbind(df, data.frame(factor=closest_BS2) %>% setNames(paste(factor, "up", sep="_")))
  }
  BS_df2 <- rbind(BS_df2, df)
}

#---- prepare the data.frame/ calculate the distance ----#
productivity <- cbind(productivity, BS_df2[,-1]) %>%
  mutate(across(starts_with(c("Per1_up", "Per2_up", "Cry1_up", "Cry2_up", "Bmal1_up", "Clock_up", "Npas2_up")), ~ as.numeric(.) - start, .names = "distance_upstream_{.col}_TSS"))

distance_promoter_lolipop <- productivity[,c(4, 13, 21:27)] %>%
  left_join(DE_res, ., by='transcript_id') %>% 
  dplyr::select(c('transcript_id', 'gene_id', 'class', 'diffexp', 'distance_upstream_Per1_up_TSS', 'distance_upstream_Per2_up_TSS', 'distance_upstream_Cry1_up_TSS', 'distance_upstream_Cry2_up_TSS', 'distance_upstream_Bmal1_up_TSS', 'distance_upstream_Clock_up_TSS', 'distance_upstream_Npas2_up_TSS')) %>%
  pivot_longer(!c('transcript_id', 'class', 'diffexp', 'gene_id'), names_to = 'factor', values_to = 'distance') %>%
  distinct(gene_id, diffexp, factor, .keep_all = TRUE) %>% # avoid plotting multiple isoform from the same gene/class/behavior
  mutate(promoter = ifelse(abs(distance) < 1000, "promoter", "Other")) %>% 
  group_by(across(c(diffexp, factor))) %>%
  mutate(promoter_ratio = sum(promoter == "promoter", na.rm = TRUE) / n()) %>%
  ungroup() %>%
  dplyr::select(c('diffexp', 'factor', 'promoter_ratio')) %>% 
  distinct(., .keep_all = TRUE) %>%
  pivot_wider(values_from = promoter_ratio, names_from = diffexp) %>%
  mutate(factor = str_remove(factor, "distance_upstream_") %>% str_remove(., "_up_TSS") %>% toupper()) %>%
  ggplot() +
  geom_segment(aes(x=factor, xend=factor, y=DOWN*100, yend=UP*100), color="grey") +
  geom_point(aes(x=factor, y=DOWN*100, color="down-regulated"), size=3 ) +
  geom_point(aes(x=factor, y=No*100, color="ns"), size=3) +
  geom_point(aes(x=factor, y=UP*100, color="up-regulated"), size=3) +
  coord_flip() +
  theme_ipsum(plot_margin = margin(0, 10, 0, 10)) +   #extrafont::loadfonts(device="win") to get rid of an error related to fonts in windows
  theme(legend.direction="vertical", legend.key.size = unit(0.01, "cm"), axis.title.x = element_text(size=10), axis.title.y = element_blank(), legend.position = "top", legend.key.width = unit(0.1, 'cm'), legend.title = element_blank()) +
  ylab("Fraction of binding sites in the promoter") +
  scale_color_manual(name='Expression',
                     breaks=c('up-regulated', 'ns', 'down-regulated'),
                     values=c('down-regulated'="#b50047", 'ns'=rgb(0.3,0.3,0.1,0.5), 'up-regulated'="#031972")) +
  guides(color = guide_legend(override.aes = list(size=1.5)))
#distance_promoter_lolipop

#--------------------#
#--- Supplemental ---#
#--------------------#
rhythmic_tx <- read.csv(paste(list_dir,"/circadian_transcripts.csv", sep=""), header=TRUE)
mis_reg <- results %>% filter(diffexp != "No")

library("ggvenn")
venn_KO_circa <- ggvenn(list(rhythmic_tx=rhythmic_tx$gene_id, mis_reg=mis_reg$gene_id), columns = c('rhythmic_tx', 'mis_reg'), show_percentage = FALSE, fill_color = c("#434299", "#822d30"), fill_alpha = 0.3, auto_scale = FALSE, text_size = 4, set_name_size=4) 

diff_exp<-read.csv("../results/perKO_sig.csv")  %>% pull(gene_id) 
venn_KO_circa <- ggvenn(list(rhythmic_tx=rhythmic_tx$gene_id, mis_reg=mis_reg$gene_id), columns = c('rhythmic_tx', 'mis_reg'), show_percentage = FALSE, fill_color = c("#434299", "#822d30"), fill_alpha = 0.3, auto_scale = FALSE, text_size = 4, set_name_size=4) 


### Work ###
ggsave(paste(fig_dir,"/4_distance_lollipop_plot.svg", sep=""), plot=distance_promoter_lolipop, width=3, height = 4)
ggsave(paste(fig_dir,"/4_nb_tx_barplot.svg", sep=""), plot=nb_tx_barplot, width=2, height = 4)
ggsave(paste(fig_dir,"/4_reactome_cnet.svg", sep=""), plot=reactome_cnet, width=7, height = 8)
ggsave(paste(fig_dir,"/4__distance_promoter_lolipop.svg", sep=""), plot=distance_promoter_lolipop, width=3, height = 4)
ggsave(paste(fig_dir,"/4_reactome_cnet.svg", sep=""), plot=reactome_cnet, width=7, height = 7)
ggsave(paste(fig_dir,"/4_venn_rhythmic_vs_misreg.svg", sep=""), plot=venn_KO_circa, width=7, height = 7)

#  volcano TE
#  TE_type <- read.table("../results/annotated_transcripts.tsv", header = TRUE) %>%
#  dplyr::filter(matching_transcript %in% filter(results, diffexp != "No")$transcript_id) %>%
#  dplyr::rename("transcript_id" = "matching_transcript")
#  TE_cov_misreg <- filter(results, diffexp != "No") %>%
#  left_join(., TE_type) %>%
#  mutate(tecov = tecov*100) %>%
#  ggplot(aes(y=-log10(padj_per_ko), x= log2FC_perko, color=tecov)) + 
#  geom_point(size=1) +
#  geom_vline(xintercept=c(-2, 2), col="black", linetype = "dashed") +
#  geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed") + 
#  theme_classic() + 
#  theme(text=element_text(size= 10), legend.position = "top",legend.box.spacing = unit(0, "pt")) +
#  xlab('log2(Expression Fold Change)') +
#  ylab('-log10(padj)') +
#  scale_color_continuous(na.value = "grey80", ) +
#  labs(color = "TE coverage")
#  TE_families <- read.table("../data/mm39_families.tsv") %>%
#  dplyr::select(V1, V2) %>%
#  setNames(c("te_name", "te_family"))
#
#  TE_type_misreg <- filter(results, diffexp != "No") %>%
#  left_join(., TE_type) %>%
#  left_join(., TE_families, by='te_name') %>%
#  filter(!is.na(tecov)) %>%
#  ggplot(aes(y=-log10(padj_per_ko), x= log2FC_perko, color=te_family)) + 
#  geom_point(size=1.5) +
#  geom_vline(xintercept=c(-2, 2), col="black", linetype = "dashed") +
#  geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed") + 
#  theme_classic() + 
#  theme(text=element_text(size= 10), legend.position = "bottom",legend.box.spacing = unit(0, "pt"), legend.title = element_blank()) +
#  xlab('log2(Expression Fold Change)') +
#  ylab('-log10(padj)') +
#  scale_color_discrete(na.value = "grey80", ) +
#  labs(color = "TE type")#