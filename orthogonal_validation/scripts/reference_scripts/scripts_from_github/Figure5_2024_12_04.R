
library(dplyr)
library(ggplot2)
library(cowplot)
library(EnhancedVolcano)
library(readxl)
library(biomaRt)
library(stringr)
library(tidyr)
library(ggrepel)
library(ggvenn)
library(org.Mm.eg.db)
library(clusterProfiler)
library(AnnotationHub)
library(scales)
library(ggpubr)

setwd("C:/Users/kmamgain/prettycoolPores_new/scripts")
#directory to save all the figures
fig_dir <- paste(getwd(), "/../Manuscript_figures_",Sys.Date(), sep="") #Path of directory where we want to save all the extra figures
theme_legend<- theme(legend.key.size = unit(0.1, 'cm'), #change legend key size
                     legend.key.height = unit(0.5, 'cm'), #change legend key height
                     legend.key.width = unit(0.1, 'cm'), #change legend key width
                     legend.title = element_blank())
flair_path<-"../results/FLAIR_2024-09-11"
biomart_path<-"../data/biomart_annotation_2024-11-05.rds"
salmon_path<-"../results/salmon_2024-09-25/"
suppa_path<- "../results/suppa_2024-11-06"
list_dir<- "../results/lists"
diff_suppa_path<- "../results/Differential_suppa_2024-11-07"

my_palette <- c()
my_theme <- theme_bw() + theme(legend.text = element_text(size=7), axis.title.y = element_text(size=9), axis.title.x = element_text(size=9), axis.text.y=element_text(size = 7), axis.text.x=element_text(size = 7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))

theme_transparent <- theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill='transparent'), legend.box.background = element_rect(fill='transparent'))

# _______ Import data _______ #

nuc_df <- read_excel("../data/220211-LC-0113-0121_Nuc.xlsx", sheet="Proteins") %>% 
  .[,c(3, 4, 17:19,23:28)] %>% 
  as.data.frame(.) %>%
  setNames(c("protein_id", "Description", "abundance_ratio 20/8", "abundance ratio 20/KO", "abundance ratio 8/KO", "pval 20/8", "pval 20/KO", "pval 8/KO", "padj 20/8", "padj 20/KO", "padj 8/KO")) %>%
  filter(!protein_id == "sp") %>% 
  mutate(gene_name = str_extract(Description, "(?<=GN=)[^ ]+"))

cyto_df <- read_excel("../data/220114-LC-0052-0060-(1)_Cyto.xlsx", sheet="Proteins") %>% 
  .[,c(3, 4, 17:19,23:28)] %>% 
  as.data.frame(.) %>%
  setNames(c("protein_id", "Description", "abundance_ratio 20/8", "abundance ratio 20/KO", "abundance ratio 8/KO", "pval 20/8", "pval 20/KO", "pval 8/KO", "padj 20/8", "padj 20/KO", "padj 8/KO")) %>%
  filter(!protein_id == "sp") %>% 
  mutate(gene_name = str_extract(Description, "(?<=GN=)[^ ]+"))

max_ln <- max(c(length(nuc_df$gene_name), length(cyto_df$gene_name))) 
data.frame(nuclei = c(nuc_df$gene_name,rep(NA, max_ln - length(nuc_df$gene_name))), 
           cytoplasm = c(cyto_df$gene_name,rep(NA, max_ln - length(cyto_df$gene_name)))) %>%
  write.table(., "../lists/Fig6B.txt", quote = FALSE, sep = "\t", row.names = FALSE)

proteomic_df <- rbind(nuc_df[,c(12, 3, 6, 9)] %>% mutate(condition="CT20/CT8_nuc") %>%  setNames(c('name', 'ratio', 'pval', 'padj', 'condition')),
                      nuc_df[,c(12, 4, 7, 10)] %>% mutate(condition="CT20/Mutant_nuc") %>%  setNames(c('name', 'ratio', 'pval', 'padj', 'condition')),
                      cyto_df[,c(12, 3, 6, 9)] %>% mutate(condition="CT20/CT8_cyto") %>%  setNames(c('name', 'ratio', 'pval', 'padj', 'condition')),
                      cyto_df[,c(12, 4, 7, 10)] %>% mutate(condition="CT20/Mutant_cyto") %>%  setNames(c('name', 'ratio', 'pval', 'padj', 'condition'))) %>%
  mutate(significance=ifelse(padj < 0.1 & ratio > 1, "down (padj)", ifelse(padj < 0.1 & ratio < 1, "up (padj)", ifelse(pval < 0.05 & ratio > 1, "down (pval)", ifelse(pval < 0.05 & ratio < 1, "up (pval)", "ns")))),
         label=ifelse(significance=="ns", "", name))

filter(proteomic_df, condition == "CT20/CT8_nuc") %>% 
  write.table(., "../lists/Fig6D.txt", quote = FALSE, sep = "\t", row.names = FALSE)
filter(proteomic_df, condition == "CT20/Mutant_nuc") %>% 
  write.table(., "../lists/Fig6I.txt", quote = FALSE, sep = "\t", row.names = FALSE)

filter(proteomic_df, condition == "CT20/CT8_cyto") %>% 
  write.table(., "../lists/Fig6F.txt", quote = FALSE, sep = "\t", row.names = FALSE)
filter(proteomic_df, condition == "CT20/Mutant_cyto") %>% 
  write.table(., "../lists/Fig6O.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# _______ Visualization _______ #
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

######################################
### Stats total detected proteins ###
######################################
venn_prot <- ggvenn(list(
  Nuclei=nuc_df$protein_id,
  Cytoplasm=cyto_df$protein_id),
  fill_alpha = 0.3, 
  stroke_size = 0.7,
  show_percentage = FALSE,
  auto_scale = TRUE,
  set_name_size = 3.5, 
  fill_color = c("white", "white"))

#####
hub <- AnnotationHub()
q <- query(hub, "org.Mm.eg.db")
id <- q$ah_id[length(q)]
org.Mm.eg.db <- hub[[id]]

ego_total_nuc <- enrichGO(gene = nuc_df$gene_name, keyType = 'SYMBOL', OrgDb = org.Mm.eg.db, ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff = 0.01) %>% as.data.frame(.) 
ego_total_cyto <- enrichGO(gene = cyto_df$gene_name, keyType = 'SYMBOL', OrgDb = org.Mm.eg.db, ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff = 0.01) %>% as.data.frame(.) 

GO_total_nuc <- ego_total_nuc %>%
  .[1:5,] %>%
  ggplot(., aes(x=-log10(p.adjust), y=factor(Description, levels=Description), col=rev(-log10(p.adjust)), size=-log10(p.adjust), label=Description)) +
  geom_text_repel(point.padding = 0, min.segment.length = 1, direction = "y", nudge_y = 0.3) +
  geom_point(alpha=0.8, size = 3) +
  scale_size(range = c(2, 4)) +
  my_theme +
  theme(legend.position="none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  ylab("-log10(pval)") +
  scale_color_gradientn(colors=c("black", "grey50")) +
  coord_cartesian(clip = "off")

GO_total_cyto <- ego_total_cyto %>%
  .[1:5,] %>%
  ggplot(., aes(x=-log10(p.adjust), y=factor(Description, levels=Description), col=rev(-log10(p.adjust)), size=-log10(p.adjust), label=Description)) +
  geom_text_repel(point.padding = 0, min.segment.length = 1, direction = "y", nudge_y = 0.3) +
  geom_point(alpha=0.8, size = 3) +
  scale_size(range = c(2, 4)) +
  my_theme +
  theme(legend.position="none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  ylab("-log10(pval)") +
  scale_color_gradientn(colors=c("black", "grey50")) +
  coord_cartesian(clip = "off")

Fig5b <- plot_grid(plot_grid(GO_total_nuc, GO_total_cyto, ncol=1), venn_prot, nrow=1)
Fig5b

rbind(as.data.frame(ego_total_nuc) %>% mutate(fraction = "nuclei", .before = 1),
      as.data.frame(ego_total_cyto) %>% mutate(fraction = "cytoplasm", .before = 1)) %>%
  write.table(., "../lists/Fig6C.txt", quote = FALSE, sep = "\t", row.names = FALSE)

######################################
### Circadian comparison in nuclei ###
######################################
## Volcano ##

candidates <- c("H1-3", "Arntl", "Clock", "Sap18", "Sf3b6", "Senp1", "Hnrnpd", "Hnrnpab", "Hnrnpdl", "Eif2s1", "Csnk2a1", "Nxt1", "Taf12", "Senp3", "Crebbp", "Rdh14", "Casq1")

volcano_time_nuc <- proteomic_df %>%
  filter(., condition == "CT20/CT8_nuc") %>%
  mutate(labels = ifelse(name %in% candidates, name, "")) %>%
  ggplot(., aes(x=-log2(ratio), y=-log10(pval), color=significance, label=labels)) +
  geom_point(size=3) +
  theme_classic() + 
  scale_color_manual(values=c("#aa9041", "#ffd862", "grey90", "#aa9041", "#ffd862", "grey90")) +
  geom_label_repel(size=3, force=2, max.overlaps = 100, min.segment.length = 0) +
  annotate(geom="text", x=-1.2, y=0, col="black", label = paste("n=", nrow(filter(proteomic_df, condition == "CT20/CT8_nuc" & significance == "down (padj)")))) +
  annotate(geom="text", x=1.2, y=0, col="black", label = paste("n=", nrow(filter(proteomic_df, condition == "CT20/CT8_nuc" & significance == "up (padj)")))) +
  theme(legend.position="none") +
  xlab("log2(ratio)") +
  theme_transparent
volcano_time_nuc

## GO ##

bckg <- nuc_df$gene_name

ego_circa_nuc <- enrichGO(gene = dplyr::filter(proteomic_df, condition == "CT20/CT8_nuc" & padj < 0.1)$name, universe = bckg, keyType = 'SYMBOL', OrgDb = org.Mm.eg.db, ont = "all", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff = 0.05) 

GO_circa_nuc <- as.data.frame(ego_circa_nuc) %>%
  .[1:5,] %>%
  ggplot(., aes(x=-log10(p.adjust), y=factor(Description, levels=Description))) +
  geom_bar(stat='identity', width = 0.6) +
  my_theme +
  theme(axis.title.y = element_blank())
GO_circa_nuc

## Example ##

list_examples <- as.data.frame(ego_circa_nuc) %>% .[2,9] %>% 
  gsub("\\[\\d+\\] \"|\"$", "", .) %>% 
  str_split("/") %>%
  unlist()

dplyr::filter(proteomic_df, name %in% list_examples) %>%
  ggplot(., aes(y=ratio, x = name)) +
  geom_bar(stat='identity', width=0.6) +
  theme_classic() +
  annotate(geom="line", x=0, y=c(-Inf, Inf)) 


######################################
### Phenotype comparison in nuclei ###
######################################

## Volcano ##

candidates <- c("Fus", "Srsf3", "Srsf4", "Srsf5", "Srsf6", "Srnpc", "Gsta4", "Suclg2")

volcano_pheno_nuc <- proteomic_df %>%
  filter(., condition == "CT20/Mutant_nuc") %>%
  mutate(labels = ifelse(name %in% candidates, name, "")) %>%
  ggplot(., aes(x=-log2(ratio), y=-log10(pval), color=significance, label=labels)) +
  geom_point(size=3) +
  theme_classic() + 
  scale_color_manual(values=c("#5680c9", "#c5d4ec", "grey90", "#5680c9", "#c5d4ec", "grey90")) +
  geom_label_repel(size=3, force=2, max.overlaps = 100, min.segment.length = 0) +
  annotate(geom="text", x=-2.5, y=0, col="black", label = paste("n=", nrow(filter(proteomic_df, condition == "CT20/Mutant_nuc" & significance == "down (padj)")))) +
  annotate(geom="text", x=2, y=0, col="black", label = paste("n=", nrow(filter(proteomic_df, condition == "CT20/Mutant_nuc" & significance == "up (padj)")))) +
  theme(legend.position="none") +
  xlab("log2(ratio)") +
  theme_transparent
volcano_pheno_nuc

## GO ##

bckg <- nuc_df$gene_name
ego_pheno_nuc_down <- enrichGO(gene = dplyr::filter(proteomic_df, condition == "CT20/Mutant_nuc" & padj < 0.1 & ratio < 1)$name, universe = bckg, keyType = 'SYMBOL', OrgDb = org.Mm.eg.db, ont = "all", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff = 0.05) 
ego_pheno_nuc_up <- enrichGO(gene = dplyr::filter(proteomic_df, condition == "CT20/Mutant_nuc" & padj < 0.1 & ratio > 1)$name, universe = bckg, keyType = 'SYMBOL', OrgDb = org.Mm.eg.db, ont = "all", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff = 0.05) 
GO_pheno_nuc_down <- as.data.frame(ego_pheno_nuc_down) %>%
  .[1:5,] %>%
  ggplot(., aes(x=-log10(p.adjust), y=factor(Description, levels=Description), col=rev(-log10(p.adjust)), size=-log10(p.adjust), label=Description)) +
  geom_text_repel(point.padding = 0, min.segment.length = 10, direction = "y", nudge_y = 0.3) +
  geom_point(alpha=0.8, size = 3) +
  scale_size(range = c(3, 6)) +
  my_theme +
  theme(legend.position="none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  ylab("-log10(pval)") +
  scale_color_gradientn(colors=c("black", "grey50")) +
  coord_cartesian(clip = "off") +
  theme_transparent
GO_pheno_nuc_down


## ratio plot ##

list_t <- dplyr::filter(proteomic_df, significance %in% c("up (padj)", "down (padj)"))

go_list <- getBM(attributes = c('mgi_symbol', 'go_id', 'name_1006', 'definition_1006', 'go_linkage_type', 'namespace_1003', 'reactome') , filters='mgi_symbol', values = list_t, mart = mart) 

plot_df <- dplyr::filter(proteomic_df, significance %in% c("up (padj)", "down (padj)")) %>%
  dplyr::filter(condition == "CT20/Mutant_nuc") %>%
  mutate(pos_lab = ifelse(ratio < 1, log2(ratio) - 0.5, log2(ratio) + 0.5),
         GO = ifelse(name %in% filter(go_list, name_1006 %in% c("RNA splicing", "mRNA processing", "mRNA splicing, via spliceosome", "negative regulation of mRNA splicing, via spliceosome"))$mgi_symbol, "RNA splicing", "Others")) %>%
  group_by(GO) %>%
  arrange(desc(log2(ratio)), .by_group = TRUE) %>%
  ungroup() %>%
  mutate(line=row.names(.))

ratio_nuc_sig <- ggplot(plot_df, aes(y=-log2(ratio), x = factor(name, levels=name), fill=GO)) +
  geom_bar(stat='identity', width=0.7) +
  theme_classic() + 
  geom_text(aes(y=-pos_lab, label=name), size=3) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.line.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none") + 
  annotate(geom="line", y=0, x=c(-Inf, Inf)) +
  ylim(-3.5, 3.5) +
  coord_flip() +
  scale_fill_manual(values = c("grey90", "#5680c9")) +
  xlab("log2(ratio)") +
  theme_transparent
 
files_splicing <- list.files(diff_suppa_path, pattern="0",full.names = TRUE) 
splicing_df <- data.frame()

for (file in files_splicing){
  event = gsub("../results/Differential_suppa_2024-11-07/res_", "", file) %>% gsub(".dpsi.temp.0", "", .)
  print(event)
  table_splicing <- read.table(file, header = TRUE) %>% 
    mutate(event=event,
           gene_id=gsub(";.*", "", Event_id)) %>% 
    setNames(c("Event_id", "dPSI", "pval", "event", "gene_id")) %>%
    dplyr::select(!Event_id)
  splicing_df <- rbind(splicing_df, table_splicing)
}

SRSF3_iCLIP <- read_excel("../data/public_datasets/supp_table2_Anko_2012.XLSX", sheet = "SRSF3") %>% as.data.frame(.) %>% .[,1]
SRSF4_iCLIP <- read_excel("../data/public_datasets/supp_table2_Anko_2012.XLSX", sheet = "SRSF4") %>% as.data.frame(.) %>% .[,1]
SRSF3_KO_liver <- read_excel("../data/public_datasets/supp_table2_Supriya_sen_2013.xls", sheet = "Combined Exon data") %>% as.data.frame(.) %>% .[-1:-20,4] %>% str_split(., "\\|") %>% unlist()
SRSF3_KO_cells_Kumar <- list(read_excel("../data/public_datasets/supp_table3_Kumar_2019.xlsx", sheet ="SRSF3 splicing events by Majiq") %>% as.data.frame(.) %>% .[-1:-2,1], read_excel("../data/public_datasets/supp_table3_Kumar_2019.xlsx", sheet ="SRSF3KO splicing events by JUM") %>% as.data.frame(.) %>% .[-1:-3,1]) %>% unlist() %>% unique()
SRSF3_KO_cells_Song <- read_excel("../data/public_datasets/supp_table3_Xiao_Song_2019.xlsx") %>% as.data.frame(.) %>% .[-1:-3,4]
SRSF3_iCLIP <- read_excel("../data/public_datasets/supp_table2_Anko_2012.XLSX", sheet = "SRSF3") %>% as.data.frame(.) %>% .[,1]
SRSF4_iCLIP <- read_excel("../data/public_datasets/supp_table2_Anko_2012.XLSX", sheet = "SRSF4") %>% as.data.frame(.) %>% .[,1]

splicing_df <- splicing_df %>%
  left_join(., getBM(attributes = c('mgi_symbol', 'ensembl_gene_id') , filters='ensembl_gene_id', values = .$gene_id, mart = mart) %>% setNames(c("name", "gene_id")), by='gene_id') %>%
  mutate(KO_liver = ifelse(tolower(name) %in% tolower(SRSF3_KO_liver), "yes", "ns"),
         KO_cells_Kumar = ifelse(tolower(name) %in% tolower(SRSF3_KO_cells_Kumar), "yes", "ns"),
         KO_cells_Song = ifelse(tolower(name) %in% tolower(SRSF3_KO_cells_Song), "yes", "ns"),
         any_KO = ifelse(KO_liver == "yes" | KO_cells_Kumar == "yes" | KO_cells_Song == "yes", "yes", "ns"),
         iCLIP_SRSF3 = ifelse(tolower(name) %in% tolower(SRSF3_iCLIP), "yes", "ns"),
         iCLIP_SRSF4 = ifelse(tolower(name) %in% tolower(SRSF4_iCLIP), "yes", "ns")) %>%
  unique()

splicing_df_sig <- filter(splicing_df, pval < 0.05) %>% filter(dPSI > 0.1 | dPSI < -0.1) 
splicing_df_others <- filter(splicing_df, !name %in% splicing_df_sig$name) 

### iCLIP analysis

events <- c("A3", "A5", "AF", "AL", "MX", "RI", "SE")

SRSF3_df <- data.frame(event = events, bckg = numeric(length(events)))

for (i in 1:length(events)) {
  bckg_value <- length(unique(filter(splicing_df_others, event == events[i] & iCLIP_SRSF3 == "yes")$name)) / length(unique(filter(splicing_df_others, event == events[i])$name))
  sig_value <- length(unique(filter(splicing_df_sig, event == events[i] & iCLIP_SRSF3 == "yes")$name)) / length(unique(filter(splicing_df_sig, event == events[i])$name))
  SRSF3_df$bckg[i] <- bckg_value
  SRSF3_df$sig[i] <- sig_value
}

SRSF3_targets_plot <- SRSF3_df %>% 
  pivot_longer(!event, values_to = "fraction", names_to = "group") %>%
  ggplot(aes(x=group, y=fraction*100, fill=group)) +
  geom_bar(stat='identity', width = 0.8, color="black") +
  facet_wrap("event", nrow=1, strip.position = "bottom") +
  scale_fill_manual(values=c("grey90", "#c5d4ec")) +
  my_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), legend.position="bottom") +
  theme_transparent +
  ylab("Fraction") +
  annotate(geom='line', y=0, x=c(-Inf, Inf))
SRSF3_targets_plot

## Fischer exact test

SRSF3_fischer_df <- data.frame()

for (i in 1:length(events)) {
  df <- data.frame(splicing=c("altered", 'others'),
                   SRSF3_binding = c(length(unique(filter(splicing_df_sig, event == events[i] & iCLIP_SRSF3 == "yes")$name)), length(unique(filter(splicing_df_others, event == events[i] & iCLIP_SRSF3 == "yes")$name))),
                   no_binding = c(length(unique(filter(splicing_df_sig, event == events[i])$name)) - length(unique(filter(splicing_df_sig, event == events[i] & iCLIP_SRSF3 == "yes")$name)), length(unique(filter(splicing_df_others, event == events[i])$name)) - length(unique(filter(splicing_df_others, event == events[i] & iCLIP_SRSF3 == "yes")$name))))
  fischer <- data.frame(fischer_pval = fisher.test(df[,2:3])$p.value,
                        event = events[i])
  SRSF3_fischer_df <- rbind(SRSF3_fischer_df, fischer)
}


#########################################
### Phenotype comparison in cytoplasm ###
#########################################

## GO ##

bckg <- cyto_df$gene_name

ego_pheno_cyto_up <- enrichGO(gene = dplyr::filter(proteomic_df, condition == "CT20/Mutant_cyto" & padj < 0.1 & ratio < 1)$name, universe = bckg, keyType = 'SYMBOL', OrgDb = org.Mm.eg.db, ont = "all", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff = 0.05) 
GO_df_up <- as.data.frame(ego_pheno_cyto_up) %>% dplyr::select(c("Description", "p.adjust", "Count", "geneID", "ONTOLOGY"))

df_up <- data.frame()
for (i in 1:nrow(GO_df_up)){
  Description <- GO_df_up[i,1]
  list_genes_GO <- str_split(GO_df_up[i,4], "\\/") %>% unlist()
  ratios <- dplyr::filter(proteomic_df, name %in% list_genes_GO & condition == "CT20/Mutant_cyto")$ratio 
  log2_mean_FC <- log2(sum(ratios)/length(ratios))
  df_up <- rbind(df_up, data.frame(Description, log2_mean_FC))
}

ego_pheno_cyto_down <- enrichGO(gene = dplyr::filter(proteomic_df, condition == "CT20/Mutant_cyto" & padj < 0.1 & ratio > 1)$name, universe = bckg, keyType = 'SYMBOL', OrgDb = org.Mm.eg.db, ont = "all", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff = 0.05) 
GO_df_down <- as.data.frame(ego_pheno_cyto_down) %>%dplyr::select(c("Description", "p.adjust", "Count", "geneID", "ONTOLOGY"))

df_down <- data.frame()
for (i in 1:nrow(GO_df_down)){
  Description <- GO_df_down[i,1]
  list_genes_GO <- str_split(GO_df_down[i,4], "\\/") %>% unlist()
  ratios <- dplyr::filter(proteomic_df, name %in% list_genes_GO & condition == "CT20/Mutant_cyto")$ratio 
  log2_mean_FC <- log2(sum(ratios)/length(ratios))
  df_down <- rbind(df_down, data.frame(Description, log2_mean_FC))
}


GO_df_up <- GO_df_up[!duplicated(GO_df_up$geneID), ] %>% # simplify the GO list by removing subclasses that are repeated
  left_join(., df_up, by="Description") 

GO_df_down <- GO_df_down[!duplicated(GO_df_down$geneID), ] %>% # simplify the GO list by removing subclasses that are repeated
  left_join(., df_down, by="Description") 


GO_pheno_cyto <- rbind(GO_df_down, GO_df_up) %>% 
  mutate(Description = str_replace(Description, ",.*", "")) %>%
  ggplot(., aes(x=-log2_mean_FC, y=-log10(p.adjust), size=Count, color=abs(log2_mean_FC), label=str_wrap(Description, width = 25))) +
  geom_point() +
  geom_text_repel(color="black", force = 4) +
  xlim(-1,1) +
  scale_size(range = c(1, 5), limits = c(5,25)) +
  scale_color_gradient(low = "grey80", high = "#5680c9") +
  my_theme +
  theme_transparent +
  theme(legend.position="none") +
  annotate(geom="line", x=0, y=c(-Inf,Inf), colour = "grey80", linewidth = 0.5, linetype=2) +
  xlab("log2(proteins mean FC)")
GO_pheno_cyto

#get the expression values of each transcript from the genes involved in the GO Terms

expression_df <- read.table("../results/lists/Fig1_abundance_no_batch_effect.txt", sep="\t", header=TRUE) %>%
  mutate(ensembl_transcript_id=.$transcript_id) %>%
  left_join(., getBM(attributes = c('mgi_symbol', 'ensembl_transcript_id') , filters='ensembl_transcript_id', values = .$ensembl_transcript_id, mart = mart), by='ensembl_transcript_id')
expression_df$mgi_symbol<- expression_df$mgi_symbol.y

dpsi_GO_up_df <- data.frame()
df_expr_up_genes_cyto <- data.frame()
for (i in 1:nrow(GO_df_up)){
  Description <- GO_df_up[i,1]
  list_expression_genes_GO <- str_split(GO_df_up[i,4], "\\/") %>% unlist()
  expr <- dplyr::filter(expression_df, mgi_symbol %in% list_expression_genes_GO) %>% dplyr::select(c("CT20_1", "CT20_2", "KO20_1", "KO20_2", "mgi_symbol")) %>% 
    rowwise() %>% mutate(WT = mean(c_across(CT20_1:CT20_2)), KO = mean(c_across(KO20_1:KO20_2)), Description=Description, )  %>% dplyr::select(c("mgi_symbol", "WT", "KO", "Description"))
  df_expr_up_genes_cyto <- rbind(df_expr_up_genes_cyto, expr)
  splicing_temp_df <- dplyr::filter(splicing_df, name %in% list_expression_genes_GO) %>%
    group_by(name) %>%
    dplyr::select(c(dPSI, pval, name)) %>%
    filter(abs(dPSI) == max(dPSI, na.rm = TRUE)) %>%
    mutate(Description = Description) %>%
    distinct(name, .keep_all = TRUE)
  dpsi_GO_up_df <- rbind(dpsi_GO_up_df, splicing_temp_df)
}

df_expr_up_genes_cyto <- df_expr_up_genes_cyto[rowSums(df_expr_up_genes_cyto[, c("WT", "KO")]) != 0, ] %>%
  mutate(log2_FC = log2(KO/WT)) 

expr_up_genes_cyto_plot <- ggplot(df_expr_up_genes_cyto, aes(x=str_wrap(Description, width = 35), y=log2_FC)) +
  geom_boxplot(outlier.shape=NA) +
  annotate(geom="line", y=0, x=c(-Inf,Inf), colour = "grey50", linewidth = 0.5) +
  ylim(-3,3) +
  my_theme +
  theme_transparent +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())

dpsi_GO_down_df <- data.frame()
df_expr_down_genes_cyto <- data.frame()
for (i in 1:nrow(GO_df_down)){
  Description <- GO_df_down[i,1]
  list_expression_genes_GO <- str_split(GO_df_down[i,4], "\\/") %>% unlist()
  expr <- dplyr::filter(expression_df, mgi_symbol %in% list_expression_genes_GO) %>% dplyr::select(c("CT20_1", "CT20_2", "KO20_1", "KO20_2", "mgi_symbol")) %>% rowwise() %>% mutate(WT = mean(c_across(CT20_1:CT20_2)), KO = mean(c_across(KO20_1:KO20_2)), Description=Description, ) %>% dplyr::select(c("mgi_symbol", "WT", "KO", "Description"))
  df_expr_down_genes_cyto <- rbind(df_expr_down_genes_cyto, expr)
  splicing_temp_df <- dplyr::filter(splicing_df, name %in% list_expression_genes_GO) %>%
    group_by(name) %>%
    dplyr::select(c(dPSI, pval, name)) %>%
    filter(abs(dPSI) == max(dPSI, na.rm = TRUE)) %>%
    mutate(Description = Description) %>%
    distinct(name, .keep_all = TRUE)
  dpsi_GO_down_df <- rbind(dpsi_GO_down_df, splicing_temp_df)
}

df_expr_down_genes_cyto <- df_expr_down_genes_cyto[rowSums(df_expr_down_genes_cyto[, c("WT", "KO")]) != 0, ] %>%
  mutate(log2_FC = log2(KO/WT)) 

expr_down_genes_cyto_plot <- ggplot(df_expr_down_genes_cyto, aes(x=str_wrap(Description, width = 35), y=log2_FC)) +
  geom_boxplot(outlier.shape=NA) +
  annotate(geom="line", y=0, x=c(-Inf,Inf), colour = "grey50", linewidth = 0.5) +
  ylim(-3,3) +
  my_theme +
  theme_transparent +
  theme(axis.title.x=element_blank() ,axis.text.x=element_blank())

plot_grid(GO_pheno_cyto, plot_grid(expr_down_genes_cyto_plot, expr_up_genes_cyto_plot, nrow=1, rel_widths = c(7,16)), ncol = 1, rel_heights = c(3,1))

ggsave(paste(fig_dir,"/5_grid1.svg", sep=""), plot_grid(GO_pheno_cyto, plot_grid(expr_down_genes_cyto_plot, expr_up_genes_cyto_plot, nrow=1, rel_widths = c(7,16)), ncol = 1, rel_heights = c(3,1)))


#########################################
### Circadian comparison in cytoplasm ###
#########################################

## GO ##

bckg <- cyto_df$gene_name

ego_circa_cyto_up <- enrichGO(gene = dplyr::filter(proteomic_df, condition == "CT20/CT8_cyto" & padj < 0.1 & ratio < 1)$name, universe = bckg, keyType = 'SYMBOL', OrgDb = org.Mm.eg.db, ont = "all", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff = 0.05) 
GO_df_up <- as.data.frame(ego_circa_cyto_up) %>% dplyr::select(c("Description", "p.adjust", "Count", "geneID", "ONTOLOGY"))

df_up <- data.frame()
for (i in 1:nrow(GO_df_up)){
  Description <- GO_df_up[i,1]
  list_genes_GO <- str_split(GO_df_up[i,4], "\\/") %>% unlist()
  ratios <- dplyr::filter(proteomic_df, name %in% list_genes_GO & condition == "CT20/CT8_cyto")$ratio 
  log2_mean_FC <- log2(sum(ratios)/length(ratios))
  df_up <- rbind(df_up, data.frame(Description, log2_mean_FC))
}

ego_circa_cyto_down <- enrichGO(gene = dplyr::filter(proteomic_df, condition == "CT20/CT8_cyto" & padj < 0.1 & ratio > 1)$name, universe = bckg, keyType = 'SYMBOL', OrgDb = org.Mm.eg.db, ont = "all", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff = 0.05) 
GO_df_down <- as.data.frame(ego_circa_cyto_down) %>%dplyr::select(c("Description", "p.adjust", "Count", "geneID", "ONTOLOGY"))

df_down <- data.frame()
for (i in 1:nrow(GO_df_down)){
  Description <- GO_df_down[i,1]
  list_genes_GO <- str_split(GO_df_down[i,4], "\\/") %>% unlist()
  ratios <- dplyr::filter(proteomic_df, name %in% list_genes_GO & condition == "CT20/CT8_cyto")$ratio 
  log2_mean_FC <- log2(sum(ratios)/length(ratios))
  df_down <- rbind(df_down, data.frame(Description, log2_mean_FC))
}


GO_df_up <- GO_df_up[!duplicated(GO_df_up$geneID), ] %>% # simplify the GO list by removing subclasses that are repeated
  left_join(., df_up, by="Description") 

GO_df_down <- GO_df_down[!duplicated(GO_df_down$geneID), ] %>% # simplify the GO list by removing subclasses that are repeated
  left_join(., df_down, by="Description") 


GO_circa_cyto <- rbind(GO_df_down, GO_df_up) %>% 
  mutate(Description = str_replace(Description, ",.*", "")) %>%
  ggplot(., aes(x=-log2_mean_FC, y=-log10(p.adjust), size=Count, color=abs(log2_mean_FC), label=str_wrap(Description, width = 25))) +
  geom_point() +
  geom_text_repel(color="black", force = 4) +
  xlim(-0.75,0.75) +
  scale_size(range = c(3, 8), limits = c(5,25)) +
  scale_color_gradient(low = "grey80", high = "#aa9041") +
  my_theme +
  theme_transparent +
  theme(legend.position="none") +
  annotate(geom="line", x=0, y=c(-Inf,Inf), colour = "grey80", linewidth = 0.5, linetype=2) +
  xlab("log2(proteins mean FC)")
GO_circa_cyto

#get the expression values of each transcript from the genes involved in the GO Terms

dpsi_GO_up_df <- data.frame()
df_expr_up_genes_cyto <- data.frame()
for (i in 1:nrow(GO_df_up)){
  Description <- GO_df_up[i,1]
  list_expression_genes_GO <- str_split(GO_df_up[i,4], "\\/") %>% unlist()
  expr <- dplyr::filter(expression_df, mgi_symbol %in% list_expression_genes_GO) %>% dplyr::select(c("CT20_1", "CT20_2", "CT8_1", "CT8_2", "mgi_symbol")) %>% rowwise() %>% mutate(night = mean(c_across(CT20_1:CT20_2)), day = mean(c_across(CT8_1:CT8_2)), Description=Description, ) %>% dplyr::select(c("mgi_symbol", "night", "day", "Description"))
  df_expr_up_genes_cyto <- rbind(df_expr_up_genes_cyto, expr)
  splicing_temp_df <- dplyr::filter(splicing_df, name %in% list_expression_genes_GO) %>%
    group_by(name) %>%
    dplyr::select(c(dPSI, pval, name)) %>%
    filter(abs(dPSI) == max(dPSI, na.rm = TRUE)) %>%
    mutate(Description = Description) %>%
    distinct(name, .keep_all = TRUE)
  dpsi_GO_up_df <- rbind(dpsi_GO_up_df, splicing_temp_df)
}

df_expr_up_genes_cyto <- df_expr_up_genes_cyto[rowSums(df_expr_up_genes_cyto[, c("night", "day")]) != 0, ] %>%
  mutate(log2_FC = log2(day/night)) 

expr_up_genes_cyto_plot <- ggplot(df_expr_up_genes_cyto, aes(x=str_wrap(Description, width = 35), y=log2_FC)) +
  geom_boxplot(outlier.shape=NA) +
  annotate(geom="line", y=0, x=c(-Inf,Inf), colour = "grey50", linewidth = 0.5) +
  ylim(-1.5,1.5) +
  my_theme +
  theme_transparent +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())

dpsi_GO_down_df <- data.frame()
df_expr_down_genes_cyto <- data.frame()
for (i in 1:nrow(GO_df_down)){
  Description <- GO_df_down[i,1]
  list_expression_genes_GO <- str_split(GO_df_down[i,4], "\\/") %>% unlist()
  expr <- dplyr::filter(expression_df, mgi_symbol %in% list_expression_genes_GO) %>% dplyr::select(c("CT20_1", "CT20_2", "CT8_1", "CT8_2", "mgi_symbol")) %>% rowwise() %>% mutate(night = mean(c_across(CT20_1:CT20_2)), day = mean(c_across(CT8_1:CT8_2)), Description=Description, ) %>% dplyr::select(c("mgi_symbol", "night", "day", "Description"))
  df_expr_down_genes_cyto <- rbind(df_expr_down_genes_cyto, expr)
  splicing_temp_df <- dplyr::filter(splicing_df, name %in% list_expression_genes_GO) %>%
    group_by(name) %>%
    dplyr::select(c(dPSI, pval, name)) %>%
    filter(abs(dPSI) == max(dPSI, na.rm = TRUE)) %>%
    mutate(Description = Description) %>%
    distinct(name, .keep_all = TRUE)
  dpsi_GO_down_df <- rbind(dpsi_GO_down_df, splicing_temp_df)
}

df_expr_down_genes_cyto <- df_expr_down_genes_cyto[rowSums(df_expr_down_genes_cyto[, c("night", "day")]) != 0, ] %>%
  mutate(log2_FC = log2(day/night)) 

expr_down_genes_cyto_plot <- ggplot(df_expr_down_genes_cyto, aes(x=str_wrap(Description, width = 35), y=log2_FC)) +
  geom_boxplot(outlier.shape=NA) +
  annotate(geom="line", y=0, x=c(-Inf,Inf), colour = "grey50", linewidth = 0.5) +
  ylim(-1.5,1.5) +
  #my_theme +  
  theme(axis.title.x = element_blank(), legend.position = "right",axis.text.x = element_text(angle = 25)) +
  theme_transparent #+
  #theme(axis.title.x=element_blank(), axis.text.x=element_blank())

plot_grid(GO_circa_cyto, plot_grid(expr_down_genes_cyto_plot, expr_up_genes_cyto_plot, nrow=1, rel_widths = c(5,11)), ncol = 1, rel_heights = c(3,1))

ggsave(paste(fig_dir,"/5_grid2.svg", sep=""), plot_grid(GO_circa_cyto, plot_grid(expr_down_genes_cyto_plot, expr_up_genes_cyto_plot, nrow=1, rel_widths = c(3,11)), ncol = 1, rel_heights = c(3,1)))

dPSI_GO_circa_cyto_down <- ggplot(dpsi_GO_down_df, aes(x=Description, y=dPSI)) +
  geom_boxplot(outlier.shape=NA) +
  annotate(geom="line", y=0.1, x=c(-Inf,Inf), colour = "grey50", linewidth = 0.5) +
  ylim(0,1) +
  my_theme +
  theme_transparent +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())

dPSI_GO_circa_cyto_up <- ggplot(dpsi_GO_up_df, aes(x=Description, y=dPSI)) +
  geom_boxplot(outlier.shape=NA) +
  annotate(geom="line", y=0.1, x=c(-Inf,Inf), colour = "grey50", linewidth = 0.5) +
  ylim(0,1) +
  my_theme +
  theme_transparent +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())

#########################
### work in progress ####
#########################

# correlation expression and protein

cor_cyto_nuc_circadian_df <- left_join(nuc_df, cyto_df, by='protein_id') %>%
  drop_na(gene_name.y) %>%
  dplyr::select(c('gene_name.y', 'abundance_ratio 20/8.x', 'abundance ratio 20/KO.x', 'abundance_ratio 20/8.y', 'abundance ratio 20/KO.y')) %>%
  setNames(c('gene_name', 'circa_nuc', 'pheno_nuc', 'circa_cyto', 'pheno_cyto'))

gene_expression_df <- expression_df %>%
  filter(rowSums(dplyr::select(., CT0_1:KO8_2)) != 0) %>%
  group_by(mgi_symbol) %>%
  summarise(across(starts_with("CT"), sum, na.rm = TRUE), across(starts_with("KO"), sum, na.rm = TRUE)) %>% 
  rowwise() %>%
  mutate(ratio_circa = mean(c(CT20_1, CT20_2))/mean(c(CT8_1, CT8_2)),
         ratio_pheno = mean(c(CT20_1, CT20_2))/mean(c(KO20_1, KO20_2))) %>%
  dplyr::select(c('mgi_symbol', 'ratio_circa', 'ratio_pheno')) %>%
  dplyr::filter(!mgi_symbol %in% cor_cyto_nuc_circadian_df$gene_name) %>%
  setNames(c('gene_name', 'ratio_circa_expr', 'ratio_pheno_expr')) %>% 
  left_join(rbind(nuc_df, cyto_df), by='gene_name') %>%
    dplyr::select(c('gene_name', 'ratio_circa_expr', 'ratio_pheno_expr', 'abundance_ratio 20/8', 'abundance ratio 20/KO')) %>%
    setNames(c('gene_name', 'ratio_circa_expr', 'ratio_pheno_expr', 'ratio_circa_prot', 'ratio_pheno_prot')) %>%
  drop_na(ratio_circa_prot)


cor_expr_prot_circa <- ggscatter(gene_expression_df, x="ratio_circa_expr", y="ratio_circa_prot", add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE) + 
  stat_cor() +
  xlab("Expression ratio") +
  ylab("Protein ratio")

cor_expr_prot_pheno <- ggscatter(gene_expression_df, x="ratio_pheno_expr", y="ratio_pheno_prot", add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE) + 
  stat_cor() +
  xlab("Expression ratio") +
  ylab("Protein ratio")

#####################
#Supplemental Figures
#####################

#SRSF4 iCLIP analysis
SRSF4_df <- data.frame(event = events, bckg = numeric(length(events)))

for (i in 1:length(events)) {
  bckg_value <- length(unique(filter(splicing_df_others, event == events[i] & iCLIP_SRSF4 == "yes")$name)) / length(unique(filter(splicing_df_others, event == events[i])$name))
  sig_value <- length(unique(filter(splicing_df_sig, event == events[i] & iCLIP_SRSF4 == "yes")$name)) / length(unique(filter(splicing_df_sig, event == events[i])$name))
  SRSF4_df$bckg[i] <- bckg_value
  SRSF4_df$sig[i] <- sig_value
}

SRSF4_targets_plot <- SRSF4_df %>% 
  pivot_longer(!event, values_to = "fraction", names_to = "group") %>%
  ggplot(aes(x=group, y=fraction*100, fill=group)) +
  geom_bar(stat='identity', width = 0.8, color="black") +
  facet_wrap("event", nrow=1, strip.position = "bottom") +
  scale_fill_manual(values=c("grey90", "#c5d4ec")) +
  my_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), legend.position="bottom") +
  theme_transparent +
  ylab("Fraction") +
  annotate(geom='line', y=0, x=c(-Inf, Inf))
SRSF4_targets_plot

## Fischer exact test

SRSF4_fischer_df <- data.frame()

for (i in 1:length(events)) {
  df <- data.frame(splicing=c("altered", 'others'),
                   SRSF4_binding = c(length(unique(filter(splicing_df_sig, event == events[i] & iCLIP_SRSF4 == "yes")$name)), length(unique(filter(splicing_df_others, event == events[i] & iCLIP_SRSF4 == "yes")$name))),
                   no_binding = c(length(unique(filter(splicing_df_sig, event == events[i])$name)) - length(unique(filter(splicing_df_sig, event == events[i] & iCLIP_SRSF4 == "yes")$name)), length(unique(filter(splicing_df_others, event == events[i])$name)) - length(unique(filter(splicing_df_others, event == events[i] & iCLIP_SRSF4 == "yes")$name))))
  fischer <- data.frame(fischer_pval = fisher.test(df[,2:3])$p.value,
                        event = events[i])
  SRSF4_fischer_df <- rbind(SRSF4_fischer_df, fischer)
}


## Expression of candidate factors

splicing_candidates <- c("Srsf3", "Srsf4", "Srsf5", "Srsf6", "Snrpc", "Fus")

expression_splicing_candidates <- dplyr::filter(expression_df, mgi_symbol %in% splicing_candidates) %>%
  dplyr::select(c("CT20_1", "CT20_2", "KO20_1", "KO20_2", "mgi_symbol")) %>%
  pivot_longer(!mgi_symbol, values_to = "expression", names_to = "condition") %>%
  mutate(condition=ifelse(condition %in% c("CT20_1", "CT20_2"), "WT", "PerKO")) %>%
  ggplot(aes(y=expression, x=factor(condition, levels=c("WT", "PerKO")), fill=factor(condition, levels=c("WT", "PerKO")))) +
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap("mgi_symbol", nrow=1, strip.position = "bottom") + 
  #ylim(0,100) +
  my_theme +
  theme_transparent + labs(fill="Condition") +
  scale_fill_manual(values=c("grey90", "#c5d4ec")) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

## Expression of up-reg proteins in KO nuc

upreg_MS <- dplyr::filter(proteomic_df, condition == "CT20/Mutant_nuc" & significance %in% c("up (padj)"))$name

expression_splicing_upreg_MS <- dplyr::filter(expression_df, mgi_symbol %in% upreg_MS) %>%
  dplyr::select(c("CT20_1", "CT20_2", "KO20_1", "KO20_2", "mgi_symbol")) %>%
  pivot_longer(!mgi_symbol, values_to = "expression", names_to = "condition") %>%
  mutate(condition=ifelse(condition %in% c("CT20_1", "CT20_2"), "WT", "PerKO")) %>%
  ggplot(aes(y=expression, x=factor(condition, levels=c("WT", "PerKO")), fill=factor(condition, levels=c("WT", "PerKO")))) +
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap("mgi_symbol", nrow=1, strip.position = "bottom") + 
  #ylim(0,10) +
  #my_theme +
  theme_transparent + labs(fill="Condition") +
  scale_fill_manual(values=c("grey90", "#c5d4ec")) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

## Expression of down-reg proteins in KO nuc

downreg_MS <- dplyr::filter(proteomic_df, condition == "CT20/Mutant_nuc" & significance %in% c("down (padj)"))$name

expression_splicing_downreg_MS <- dplyr::filter(expression_df, mgi_symbol %in% downreg_MS) %>%
  dplyr::select(c("CT20_1", "CT20_2", "KO20_1", "KO20_2", "mgi_symbol")) %>%
  pivot_longer(!mgi_symbol, values_to = "expression", names_to = "condition") %>%
  mutate(condition=ifelse(condition %in% c("CT20_1", "CT20_2"), "WT", "PerKO")) %>%
  ggplot(aes(y=expression, x=factor(condition, levels=c("WT", "PerKO")), fill=factor(condition, levels=c("WT", "PerKO")))) +
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap("mgi_symbol", nrow=1, strip.position = "bottom") + 
  #ylim(0,100) +
  #my_theme +
  theme_transparent +labs(fill="Condition") +
  scale_fill_manual(values=c("grey90", "#c5d4ec")) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

## Expression of up-reg proteins in nuc at CT20

upreg_MS_circ <- dplyr::filter(proteomic_df, condition == "CT20/CT8_nuc" & significance %in% c("up (padj)"))$name

expression_splicing_upreg_MS_circ <- dplyr::filter(expression_df, mgi_symbol.x %in% upreg_MS_circ) %>%
  dplyr::select(c("CT20_1", "CT20_2", "KO20_1", "KO20_2", "mgi_symbol.x")) %>%
  pivot_longer(!mgi_symbol.x, values_to = "expression", names_to = "condition") %>%
  mutate(condition=ifelse(condition %in% c("CT20_1", "CT20_2"), "WT", "PerKO")) %>%
  ggplot(aes(y=expression, x=factor(condition, levels=c("WT", "PerKO")), fill=factor(condition, levels=c("WT", "PerKO")))) +
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap("mgi_symbol.x", nrow=1, strip.position = "bottom") + 
  #ylim(0,100) +
  #my_theme +
  theme_transparent +labs(fill="Condition") +
  scale_fill_manual(values=c("grey90", "#c5d4ec")) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

## Expression of down-reg proteins in nuc at CT8

downreg_MS_circ <- dplyr::filter(proteomic_df, condition == "CT20/CT8_nuc" & significance %in% c("down (padj)"))$name

expression_splicing_downreg_MS_circ <- dplyr::filter(expression_df, mgi_symbol %in% downreg_MS_circ) %>%
  dplyr::select(c("CT20_1", "CT20_2", "KO20_1", "KO20_2", "mgi_symbol")) %>%
  pivot_longer(!mgi_symbol, values_to = "expression", names_to = "condition") %>%
  mutate(condition=ifelse(condition %in% c("CT20_1", "CT20_2"), "WT", "PerKO")) %>%
  ggplot(aes(y=expression, x=factor(condition, levels=c("WT", "PerKO")), fill=factor(condition, levels=c("WT", "PerKO")))) +
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap("mgi_symbol", nrow=1, strip.position = "bottom") + 
  #ylim(0,100) +
  #my_theme +
  theme_transparent +labs(fill="Condition") +
  scale_fill_manual(values=c("grey90", "#c5d4ec")) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), text= element_text(size=10))

## dPSI plot for the genes involved in the cytoplasmic GO terms

dPSI_GO_pheno_cyto_down <- ggplot(dpsi_GO_down_df, aes(x=Description, y=dPSI)) +
  geom_boxplot(outlier.shape=NA) +
  annotate(geom="line", y=0.1, x=c(-Inf,Inf), colour = "grey50", linewidth = 0.5) +
  ylim(0,1) +
  my_theme +
  theme_transparent +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())

dPSI_GO_pheno_cyto_up <- ggplot(dpsi_GO_up_df, aes(x=Description, y=dPSI)) +
  geom_boxplot(outlier.shape=NA) +
  annotate(geom="line", y=0.1, x=c(-Inf,Inf), colour = "grey50", linewidth = 0.5) +
  ylim(0,1) +
  my_theme +
  theme_transparent +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())

dPSI_GO_circa_cyto_down <- ggplot(dpsi_GO_down_df, aes(x=Description, y=dPSI)) +
  geom_boxplot(outlier.shape=NA) +
  annotate(geom="line", y=0.1, x=c(-Inf,Inf), colour = "grey50", linewidth = 0.5) +
  ylim(0,1) +
  my_theme +
  theme_transparent +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())

dPSI_GO_circa_cyto_up <- ggplot(dpsi_GO_up_df, aes(x=Description, y=dPSI)) +
  geom_boxplot(outlier.shape=NA) +
  annotate(geom="line", y=0.1, x=c(-Inf,Inf), colour = "grey50", linewidth = 0.5) +
  ylim(0,1) +
  #my_theme +
  theme_transparent #+
  #theme(axis.title.x=element_blank(), axis.text.x=element_blank())


# correlation nucleus cytoplasm

cor_cyto_nuc_pheno_plot <- ggscatter(cor_cyto_nuc_circadian_df, x="pheno_nuc", y="pheno_cyto", add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE) + 
  stat_cor() +
  xlab("Nucleus") +
  ylab("Cytoplasm")

cor_cyto_nuc_circa_plot <- ggscatter(cor_cyto_nuc_circadian_df, x="circa_nuc", y="circa_cyto", add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE) + 
  stat_cor() +
  xlab("Nucleus") +
  ylab("Cytoplasm")


###################################
### Alternative / not retained ####
###################################

GO_total_nuc <- ego_total_nuc %>%
  .[1:5,] %>%
  ggplot(., aes(x=-log10(p.adjust), y=rev(factor(Description, levels=Description)))) +
  geom_bar(stat='identity', width = 0.6) +
  my_theme +
  theme(axis.title.y = element_blank())

GO_total_cyto <- ego_total_cyto %>%
  .[1:5,] %>%
  ggplot(., aes(x=-log10(p.adjust), y=factor(Description, levels=Description))) +
  geom_bar(stat='identity', width = 0.6) +
  my_theme +
  theme(axis.title.y = element_blank())

### SRSF3 investigations ###

SRSF3_KO_liver <- read_excel("../data/public_datasets/supp_table2_Supriya_sen_2013.xls", sheet = "Combined Exon data") %>% as.data.frame(.) %>% .[-1:-20,4] %>% str_split(., "\\|") %>% unlist()
SRSF3_KO_cells_Kumar <- list(read_excel("../data/public_datasets/supp_table3_Kumar_2019.xlsx", sheet ="SRSF3 splicing events by Majiq") %>% as.data.frame(.) %>% .[-1:-2,1], read_excel("../data/public_datasets/supp_table3_Kumar_2019.xlsx", sheet ="SRSF3KO splicing events by JUM") %>% as.data.frame(.) %>% .[-1:-3,1]) %>% unlist() %>% unique()
SRSF3_KO_cells_Song <- read_excel("../data/public_datasets/supp_table3_Xiao_Song_2019.xlsx") %>% as.data.frame(.) %>% .[-1:-3,4]
SRSF3_iCLIP <- read_excel("../data/public_datasets/supp_table2_Anko_2012.XLSX", sheet = "SRSF3") %>% as.data.frame(.) %>% .[,1]
SRSF4_iCLIP <- read_excel("../data/public_datasets/supp_table2_Anko_2012.XLSX", sheet = "SRSF4") %>% as.data.frame(.) %>% .[,1]


plot_grid(ggplot(splicing_df_others, aes(x = event, y = ..count../sum(..count..), fill =KO_liver)) +
            geom_bar(position = 'fill', stat = 'count') +
            ylab("Fraction") +
            scale_y_continuous(labels = scales::percent_format(), limits = c(0,1)), 
          ggplot(splicing_df_sig, aes(x = event, y = ..count../sum(..count..), fill =KO_liver)) +
            geom_bar(position = 'fill', stat = 'count') +
            ylab("Fraction") +
            scale_y_continuous(labels = scales::percent_format(), limits = c(0,1)), nrow=1)


################
##### Work #####
################

downreg_MS_circ <- dplyr::filter(proteomic_df, condition == "CT20/CT8_nuc" & significance %in% c("down (padj)"))$name

expression_volcano_downreg_MS_circ <- dplyr::filter(expression_df, mgi_symbol %in% downreg_MS_circ) %>%
  dplyr::select(c("CT20_1", "CT20_2", "KO20_1", "KO20_2", "mgi_symbol", "ensembl_transcript_id")) %>% 
  rowwise() %>%
  mutate(WT = mean(c(CT20_1, CT20_2)),
         Mutant = mean(c(KO20_1, KO20_2)),
         ratio = log2(Mutant/WT)) %>% 
  dplyr::select(WT, Mutant, mgi_symbol, ensembl_transcript_id, ratio) %>% 
  pivot_longer(!c(mgi_symbol, ensembl_transcript_id, ratio), values_to = "expression", names_to = "condition") 


filter(expression_volcano_downreg_MS_circ, condition == "WT") %>% 
  ggplot(aes(y=expression, x=ratio, col=mgi_symbol)) +
  geom_point() 

ggsave(paste(fig_dir,"/5_Fig5A.svg", sep=""), plot=volcano_time_nuc, width=5, height = 5)
ggsave(paste(fig_dir,"/5_Fig5C.svg", sep=""), plot=expression_splicing_downreg_MS_circ, width=, height = 3)
ggsave(paste(fig_dir,"/5_Fig5E.svg", sep=""), plot=volcano_pheno_nuc, width=5, height = 5)
ggsave(paste(fig_dir,"/5_Fig5H.svg", sep=""), plot=GO_pheno_nuc_down, width=4, height = 3)
ggsave(paste(fig_dir,"/5_Fig5F.svg", sep=""), plot=expression_splicing_upreg_MS, width=5, height = 3)
ggsave(paste(fig_dir,"/5_Fig5Fbis.svg", sep=""), plot=expression_splicing_downreg_MS, width=8, height = 3)
ggsave(paste(fig_dir,"/5_expr_down_genes_cyto_plot.svg", sep=""), plot=expr_down_genes_cyto_plot, width=4, height = 4)
ggsave(paste(fig_dir,"/5_expression_splicing_candidates.svg", sep=""), plot=expression_splicing_candidates, width=4, height = 3)



Fig4D <- dplyr::select(nuc_df, c('gene_name', 'padj 20/KO', 'abundance ratio 20/KO')) %>%
  setNames(c('gene_name', 'padj', 'abundance_ratio')) %>% 
  ggplot(aes(x=reorder(gene_name, -abundance_ratio), y=log2(abundance_ratio), col=padj)) +
  geom_bar(stat='identity') +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank()) +
  scale_color_gradient(low = "darkblue", high = "grey80") +
  ylab("log2(WT/PerKO)")

Fig4D_inset <- dplyr::select(nuc_df, c('gene_name', 'padj 20/KO', 'abundance ratio 20/KO')) %>% 
  setNames(c('gene_name', 'padj', 'abundance_ratio')) %>% 
  filter(abundance_ratio < 0.5) %>%
  ggplot(aes(x=reorder(gene_name, -abundance_ratio), y=log2(abundance_ratio), col=padj, fill=padj)) +
  geom_bar(stat='identity', width=0.6) +
  theme_classic() +
  theme(axis.title.y = element_blank(), axis.text.x=element_text(angle=90), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(), legend.position = 'none') +
  scale_color_gradient(low = "darkblue", high = "grey80") +
  scale_fill_gradient(low = "darkblue", high = "grey80") +
  scale_x_discrete(position = "top") 
Fig4D_inset

Fig4D_inset2 <- dplyr::select(nuc_df, c('gene_name', 'padj 20/KO', 'abundance ratio 20/KO')) %>% 
  setNames(c('gene_name', 'padj', 'abundance_ratio')) %>% 
  filter(abundance_ratio > 3) %>%
  ggplot(aes(x=reorder(gene_name, -abundance_ratio), y=log2(abundance_ratio), col=padj, fill=padj)) +
  geom_bar(stat='identity', width=0.6) +
  theme_classic() +
  theme(axis.title.y = element_blank(), axis.text.x=element_text(angle=90), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(), legend.position = 'none') +
  scale_color_gradient(low = "darkblue", high = "grey80") +
  scale_fill_gradient(low = "darkblue", high = "grey80") +
  scale_x_discrete(position = "top") 

Fig4D <- ggdraw() + 
  draw_plot(Fig4D) +
  draw_plot(Fig4D_inset, x = 0.48, y = .6, width = .4, height = .35) +
  draw_plot(Fig4D_inset2, x = 0.07, y = 0.01, width = .4, height = .35)

ggsave(paste(fig_dir,"/Fig4D.pdf",sep=""), plot=Fig4D, width=10, height = 5)

dplyr::filter(expression_df, mgi_symbol %in% downreg_MS_circ) %>%
  dplyr::select(c("CT20_1", "CT20_2", "KO20_1", "KO20_2", "mgi_symbol")) %>%
  group_by(mgi_symbol) %>% mutate(iso=row_number()) %>%
  pivot_longer(!c(mgi_symbol, iso), values_to = "expression", names_to = "condition") %>%
  mutate(condition=ifelse(condition %in% c("CT20_1", "CT20_2"), "WT", "Mutant")) %>%
  ggplot(aes(y=expression, x=condition, col=as.character(iso))) +
  geom_point() + 
  facet_wrap("mgi_symbol", nrow=3, strip.position = "bottom") + 
  ylim(0,100) +
  my_theme +
  theme_transparent + theme(legend.position='none')

dplyr::filter(expression_df, mgi_symbol %in% upreg_MS_circ) %>%
  dplyr::select(c("CT20_1", "CT20_2", "KO20_1", "KO20_2", "mgi_symbol")) %>%
  group_by(mgi_symbol) %>% mutate(iso=row_number()) %>%
  pivot_longer(!c(mgi_symbol, iso), values_to = "expression", names_to = "condition") %>%
  mutate(condition=ifelse(condition %in% c("CT20_1", "CT20_2"), "WT", "Mutant")) %>%
  ggplot(aes(y=expression, x=condition, col=as.character(iso))) +
  geom_point() + 
  facet_wrap("mgi_symbol", nrow=3, strip.position = "bottom") + 
  ylim(0,100) +
  my_theme +
  theme_transparent + theme(legend.position='none')