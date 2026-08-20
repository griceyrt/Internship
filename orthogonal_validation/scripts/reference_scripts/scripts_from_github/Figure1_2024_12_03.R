#Loading required

library(UpSetR)
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
library(viridis)
library(GenomicFeatures)
library(scales)
library(patchwork)
library(ggridges)
library(clusterProfiler)
library(tidyverse)
library(stringr)
library(rlang)
library(cowplot)
library(extrafont)
library(ggpubr)
library(svglite)

#col_test <- c("#e2ddbc","#00263c", "#36857a") or try col_test <- c("#9e9784","#00263c", "#36857a")
setwd("C:/Users/kmamgain/prettycoolPores_new/scripts")   #change here
 
fdr_cutoff <- 0.05 #false discovery rate cutoff
base_size<-13  #size of the text for figures 
retrieve_from_biomart <- FALSE #do you want to retrive annotation from biomart?
tx_names_with_version <- TRUE #Transcripts are with version  or not 
class_names <- c("Annotated", "Novel-Annotated", "Novel-Novel") #list of three classes or 'categories'
count_threshold<-5 #threshold for number of reads
fig_dir <- paste(getwd(), "/../Manuscript_figures_",Sys.Date(), sep="") #Path of directory where we want to save all the extra figures
flair_path<- "../results/FLAIR_2024-09-11"
biomart_path<-"../data/biomart_annotation_2024-11-05.rds"
salmon_path<-"../results/salmon_2024-09-25/"
suppa_path<- "../results/suppa_2024-11-06"
list_dir<- "../results/lists"

txdb <- makeTxDbFromGFF(paste(flair_path,"/transcriptome_productivity.gtf", sep=""), format = "gtf") #Load the transcriptome gtf
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


#Uploading sf files from salmon quantification
files <- list.files(path = salmon_path, pattern = ".sf", 
                    full.names = TRUE, include.dirs = FALSE, no.. = TRUE, recursive = TRUE)
names(files) <- str_extract(files, "([CK][TO]\\d+_\\d)(\\.sf)$", group = 1) 
txi.salmon <- tximport(files, type = "salmon", txOut = TRUE, countsFromAbundance = "no", ignoreTxVersion = TRUE, 
                       importer = function(x) read.delim(x))
exp_design <- data.frame(batch = factor(str_extract(colnames(txi.salmon$abundance), "(?<=_)\\d")),
                         time = as.numeric(str_extract(colnames(txi.salmon$abundance), "\\d+")),
                         sample_id= colnames(txi.salmon$abundance))
rownames(exp_design) <- exp_design$sample_id
counts <- magrittr::set_colnames(txi.salmon$counts,exp_design$sample_id)

if (tx_names_with_version){
  rownames(counts) <- str_extract(rownames(counts), "([A-Za-z0-9-]+)")
  }
mode(counts) <- "integer"


##---------Expression levels in  transcriptome--------------------##
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = exp_design,
                              design= ~ batch)

rowData(dds) <- annotation[rownames(dds),]

keep <- rowSums(counts > count_threshold) >= 2   ##this is a universal threshold

dds <- dds[keep, ]

abundance_no_batch_effect <- txi.salmon$abundance %>% data.frame %>%
        magrittr::set_rownames(str_extract(rownames(.), "([A-Za-z0-9-]+)")) %>% 
        {log2(. + 0.001)} %>%
  {limma::removeBatchEffect(.[rownames(dds), ], batch = colData(dds)$batch)}

rowData(dds) %>% as_tibble() %>% 
  group_by(gene_id) %>% 
  mutate(mgi_symbol = ifelse(grepl("^ENS", gene_id), 
                             unique(na.omit(mgi_symbol)), mgi_symbol)) %>%
  write_excel_csv(file = "../results/expr_tx_list.csv")

#----------------------productivity by flair------------------------------##
productivity <- read.table(paste(flair_path,"/productivity.txt", sep=""), sep="\t") %>% mutate(transcript_id= .$V4, productivity_flair= .$V13) %>% left_join(., annotation) %>% 
      mutate(productivity_flair = fct_relevel(productivity_flair, c("PRO","PTC","NGO","NST")),
      productivity_flair = fct_recode(productivity_flair, 
                                !!!c("productive" = "PRO",
                                     "premature termination codon" = "PTC",
                                     "no start codon" = "NGO",
                                     "has start but no stop codon" = "NST")))%>% dplyr::select(c("transcript_id", "gene_id","productivity_flair", "class")) %>% 
      filter(transcript_id %in% row.names(dds))

df_plot <- productivity %>% 
   group_by(class, productivity_flair) %>% 
   summarize(number = dplyr::n()) %>%
   mutate(fraction = number/sum(number),
          ymax = cumsum(fraction),
          ymin = c(0, head(ymax, n=-1)),
          labelpos = 0.5*(ymax + ymin))
 
df_annot <- productivity %>% group_by(class) %>% 
               summarize(genes = n_distinct(gene_id), 
                     tx = n_distinct(transcript_id)) %>% 
              mutate(lab = paste(tx, "Tx\n", genes, "genes"))

productivity_plot <- ggplot(df_plot) +
   geom_rect(aes(xmin=3, xmax=4, ymin=ymin, ymax=ymax, alpha = fct_inorder(productivity_flair), fill=class), color="grey30") +
   geom_text(aes(x = 3.5, y = labelpos, label = number), nudge_x = 1.0, size=7/.pt) + 
   coord_polar(theta = "y") + xlim(c(1.5,4.5)) +
   geom_text(aes(x=1.5, y=0, label=lab), data=df_annot, size=8/.pt) + 
   facet_wrap(~class, nrow=1) + 
   theme_void(base_size = base_size)+ 
   scale_fill_brewer(palette="Set2") + #scale_fill_manual(values = c("#8da0cb", "#66c2a5", "#fc8d62")) +
   scale_alpha_manual(values = c(1.0, 0.5, 0.3,0.05), name="") +                                                        
   theme(legend.key.size = unit(5, 'pt'), legend.position = "right",legend.title = element_blank(), strip.placement = "inside", strip.text = element_text(face="bold",size = 10))

nb_exons_plot <- ggplot(annotation[rownames(dds), ], 
                        aes(group=class, y=nexon, fill=class)) + 
  geom_histogram(binwidth=1, center=0, position="dodge") + 
  scale_fill_brewer(palette="Set2", name="")  + theme_classic(base_size = base_size) + 
  ylab("number of exons") + xlab("number of isoforms") +
  theme(legend.position= "none", aspect.ratio = 1.0, legend.key.size = unit(15,"pt")) +  #legend position was c(0.7, 0.85) before
  scale_y_continuous(limits=c(1,24), expand=expansion(), breaks = c(1,seq(4,24,4))) +
  scale_x_continuous(expand = expansion(mult = 0.01))

tx_length_plot <- ggplot(annotation[rownames(dds), ]) + 
  geom_density_ridges(aes(y=class, x=tx_len/1000, fill=class), scale=1.6) + 
  scale_fill_brewer(palette="Set2") + theme_ridges() + xlim(c(0,8)) + 
  theme_classic(base_size = base_size) + scale_y_discrete(expand = expansion(mult = 0.05)) +
  xlab("transcript length (knt)") + ylab("density") + theme(aspect.ratio = 1.0, legend.position = "none", axis.text.y = element_blank())

##----------read distribution of the samples-------------------##
abundance_by_class <- abundance_no_batch_effect %>% 
                  rowMeans() %>% 
                  data.frame(value = .) %>% 
                  rownames_to_column(var="transcript_id") %>% 
                  left_join(annotation[c("transcript_id", "class")]) %>%
                  mutate(class = factor(class, levels = class_names)) %>%
                  filter(!is.na(class))     # 300 tx have a tx_id but no information on ensembl, so filtered out

abundance_by_class_plot <- ggplot(abundance_by_class, aes(y=value, x=class, fill=class)) +
  geom_violin() + scale_fill_brewer(palette = "Set2") +  ylim(c(0, 15)) + 
  theme_classic(base_size = base_size) + ylab(bquote(transcript~expression~(log[2]~TPM))) + 
  xlab("") + theme(legend.position = "none", axis.ticks.x = element_blank(), axis.text.x = element_blank(), aspect.ratio = 1.0)

abundance_by_class_plot_supp <- ggplot(abundance_by_class, aes(y=2^value, x= cond, fill=class), stat=unique, alpha = 0.6) + geom_boxplot(outlier.shape=NA, show.legend = FALSE, alpha = 0.6) + scale_fill_manual(values = c("#b9d5ff", "#b1e086", "#fd9997")) + ylim(0,50) + theme_classic() + ylab("Transcript Expression (CPM)") + scale_x_discrete(labels=c("CT"="WT")) + xlab("")

##----------read dsitribution per class-------------------##

number_reads <- counts(dds) %>% rowSums() %>% aggregate(by=list(annotation[rownames(dds), "class"]), sum)
names(number_reads) <- c("class", "number_reads")

nb_of_reads_per_class <- number_reads %>%
                         mutate(fraction=number_reads/sum(number_reads),
                                ymax = cumsum(fraction),
                                ymin = c(0, head(ymax, n=-1)),
                                labelpos=0.5*(ymax + ymin),
                                class = factor(class, levels = class_names))


nb_of_reads_per_class_plot <- ggplot(nb_of_reads_per_class, aes(xmin=3, xmax=4, ymin=ymin, ymax=ymax)) + 
  geom_rect(aes(fill=class), color="grey30") +
  geom_text(x=3.5, aes(y=labelpos, label=sprintf("%1.0f%%", fraction*100)), color="white", size=7/.pt, fontface="bold") + 
  scale_fill_brewer(palette = "Set2", name="") + 
  geom_text(x=1, y=0, label="Reads mapped\n to each class", size=8/.pt) + 
  coord_polar(theta="y") + xlim(c(1,4)) +
  theme_void(base_size = base_size) + xlab("") + ylab("") + 
  theme(legend.position = "none")


##suppa
suppa_files <- c(AS3 = paste(suppa_path, "/transcriptome_ext_A3_strict.ioe", sep=""),
                 AS5 = paste(suppa_path, "/transcriptome_ext_A5_strict.ioe", sep=""),
                 AF = paste(suppa_path, "/transcriptome_ext_AF_strict.ioe", sep=""),
                 AL =paste(suppa_path, "/transcriptome_ext_AL_strict.ioe", sep=""),
                 MX = paste(suppa_path, "/transcriptome_ext_MX_strict.ioe", sep=""),
                 RI = paste(suppa_path, "/transcriptome_ext_RI_strict.ioe", sep=""),
                 ES = paste(suppa_path, "/transcriptome_ext_SE_strict.ioe", sep=""))

expr_tx_list <- rownames(abundance_no_batch_effect)
suppa_data <- imap(suppa_files, function(f, x) {
                    read.table(f, header=TRUE) %>% rowwise() %>%
                             mutate(total_transcripts = strsplit(as.character(total_transcripts), ","),
                                    alternative_transcripts = strsplit(as.character(alternative_transcripts), ","),
                                    complement = list(setdiff(unlist(total_transcripts), unlist(alternative_transcripts))),
                                    in_tx_list = any(unlist(complement) %in% expr_tx_list) & any(unlist(alternative_transcripts) %in% expr_tx_list),
                                    novel_event = ((length(complement)==1) & all(!grepl("^ENSMUST", unlist(complement)))) | ((length(alternative_transcripts) == 1) & all(!grepl("^ENSMUST", unlist(alternative_transcripts))))) %>% 
                             filter(in_tx_list) %>%
                             mutate(event_type = x,
                                    list_transcripts = case_when(
                               event_type == "MX" ~ list(total_transcripts),
                               event_type == "ES" ~ list(complement),
                               TRUE ~ list(alternative_transcripts)
                             )) %>%
                             dplyr::select(-alternative_transcripts, -complement, -total_transcripts) %>% 
                             unnest(cols=c(list_transcripts)) %>%
                             dplyr::rename(tx_id = list_transcripts) %>% 
                             mutate(event_type = x)}) %>%
              bind_rows() %>%
              dplyr::select(-seqname) %>%
              filter(tx_id %in% expr_tx_list)

summarise(suppa_data, total_genes = n_distinct(gene_id), total_tx = n_distinct(tx_id), total_events = n_distinct(event_id))

suppa_data %>% group_by(event_type) %>% summarize(event_by_type = n_distinct(event_id))

upset_data <- suppa_data %>% dplyr::select(-event_id) %>%  
              mutate(value=1) %>% 
              pivot_wider(id_cols = tx_id, names_from = event_type, values_from = value, values_fn = function(x) as.integer(length(x)>0), values_fill = 0)
write.csv(upset_data, "../data/upset_all.csv" )

suppa_data$event_type <- factor(suppa_data$event_type, levels = colnames(upset_data[,-1])[order(colSums(upset_data[,-1]))])

upset_data %>% mutate(n_tx = rowSums(across(where(is.numeric))))

Upset_plot <- upset(as.data.frame(upset_data), main.bar.color = "grey30", sets.x.label = "total number of isoforms", 
                    mainbar.y.label = "number of isoforms", 
                    text.scale = 1.1, 
                    query.legend = "top", 
                    nsets = 7, 
                    sets.bar.color = "grey50", 
                    empty.intersections = FALSE,
                    show.numbers = "no", 
                    point.size = 2.5, 
                    line.size = 0.9, 
                    queries = list(list(query = function(row) grepl("^ENSMUST", row["tx_id"]), 
                                        color = "#66c2a5", active = T)))

                                        ####factor for the x axis
novel_event_plot <- ggplot(suppa_data %>% distinct(event_id, .keep_all = TRUE)) + 
                        stat_count(aes(x=factor(event_type,levels=c("AF", "AL", "AS3","AS5", "RI","MX", "ES")), fill=novel_event), width = 0.5) + 
                        scale_fill_manual("",labels=c("novel events"),values= c("TRUE"= "#3576aa"))+
                        #scale_fill_colorblind(name = "", breaks = c(TRUE), labels = "novel events") + 
                        theme_classic(base_size = 13) + 
                        theme(legend.position = c(0.5, 1.0), axis.text.x = element_text(margin = margin(t=10)), 
                              axis.line.x = element_blank(), axis.ticks.x = element_blank(), 
                              aspect.ratio = 1.15, axis.title.y.right = element_text(margin = margin(l=10)), 
                              legend.key.size = unit(12, "pt")) + xlab("") +
                        ylab("number of unique AS events") + scale_y_continuous(expand = expansion(), position = "right")
ggsave(paste(fig_dir,"/1_novel_event.svg",sep=""), plot=novel_event_plot, width=3, height = 3, create.dir = TRUE)

unique_event_plot <- ggplot(suppa_data %>% group_by(gene_id) %>% summarize(n_unique = n_distinct(event_id))) + 
                        stat_count(aes(x=n_unique), width = 0.5) + 
                        theme_classic(base_size = base_size) + scale_x_continuous(limits = c(0.2,20), expand = expansion()) +
                        xlab("number of unique AS events") + ylab("number of genes")

sizes_plot <- ggplot(data.frame(event_type = colnames(upset_data[,-1]), counts = colSums(upset_data[,-1])) %>% mutate(event_type = factor(event_type, levels = levels(suppa_data$event_type)))) +
              geom_col(aes(x=factor(event_type,levels=c("AF", "AL", "AS3","AS5", "RI","MX", "ES")), y = counts), fill = "gray50", width = 0.5) + 
              theme_classic(base_size = base_size) + 
              scale_y_reverse(position = "right", expand=expansion()) + xlab("") + 
              ylab("number of isoforms") + 
              theme(axis.line.x = element_blank(), axis.text.x=element_blank(),
                    axis.ticks.x = element_blank(), axis.title.y.right = element_text(margin = margin(l=10)),
                    plot.margin = margin(b=5), aspect.ratio = 1.15)
ggsave(paste(fig_dir,"/1_sizesplot.svg",sep=""), plot=sizes_plot, width=3, height = 3)

upset_plot_grob <- ((wrap_ggplot_grob(Upset_plot$Main_bar) + inset_element(unique_event_plot, 0.25, 0.4, 0.95, 0.95, ignore_tag = TRUE)) + wrap_ggplot_grob(Upset_plot$Matrix) + (novel_event_plot/sizes_plot) + 
                      plot_layout(heights = c(1,0.4), widths = c(1.0, 0.4), 
                                  design = c(area(1,1),area(2,1), area(1,2,2,2))))

##TE--plots##te_coverage.rds,onlyTEs.rds sent by F
colors <-c("Novel-novel" = "#8da0cb", "Annotated" = "#66c2a5", "Novel-annotated" = "#fc8d62")
te_coverage <- readRDS("../data/te_coverage.rds")
onlyTEs <-readRDS("../data/onlyTEs.rds")#
#
te_coverage$status = factor(te_coverage$status, sort(unique(te_coverage$status), decreasing = T))
plot1 <- ggplot(data = te_coverage, aes(x = trcov, y = status, fill = status)) + geom_violin(scale = "width", show.legend = F) + xlab("Repeat content/transcript") + ylab("Annotation status") + scale_fill_manual(values = colors) + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),                                                                                                                                                                                                                                            panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=8), axis.ticks = element_blank(), plot.margin = margin(0,0,0,20, "pt") )
plot2 <- as_grob(ggplot(data = te_coverage[te_coverage$status != "Novel-novel",], aes(x = trcov, y = status, fill = status)) + geom_violin(scale = "width", show.legend = F, bw= 0.002)  + coord_cartesian(xlim = c(0, 0.045)) + ylab(NULL) + xlab(NULL) + 
                   scale_x_continuous(breaks = c(0, 0.045))+scale_fill_manual(values = colors)+ theme(legend.position = "None", axis.text.y = element_blank(), text = element_text(size = 8) , axis.ticks = element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_rect(colour = "black", fill=NA, size=0.8)))
repeat_plot1 <- plot1 + draw_grob(plot2, x = 0.25, y = 1.6, height = 2.1, width = 0.70)
repeat_plot1
tmp <- aggregate(onlyTEs[onlyTEs$expressed, "x"], by = list("fam" = onlyTEs[onlyTEs$expressed, "fam"]), FUN = sum)
tmp <- tmp[order(tmp$x),]
tmp <- setNames(tmp$x, tmp$fam)
labels_plot <- names(tmp)
not_expressed <- unique(onlyTEs$fam)[!(unique(onlyTEs$fam) %in% labels_plot)]##

add_space <- function(label, max, table_fam) {
  spaces_to_add <- table_fam[match(label, names(table_fam))] * 1.9
  if(spaces_to_add < 0) { spaces_to_add = 0}
return(paste(strrep(" ", spaces_to_add), label, sep = ""))
}
max_axis = 31
counts = c(rep(0, length(not_expressed)), tmp)
names(counts) = c(not_expressed, labels_plot)
new_fam_names = sapply(unique(onlyTEs$fam), FUN = add_space, max_axis, counts)
names(new_fam_names) = unique(onlyTEs$fam)
onlyTEs$fam = new_fam_names[match(onlyTEs$fam, names(new_fam_names))]
plot1 <- ggplot(data = onlyTEs[onlyTEs$expressed,], aes(x = x, y = reorder(fam, -x, sum), fill = status)) + geom_col(show.legend = FALSE, width=0.7) + xlab("Number of tx that are 80% repeat") + ylab(NULL) + scale_fill_manual(values = colors) + coord_cartesian(xlim = c(0, max_axis)) + theme(plot.margin = margin(0,0,0,20, "pt"), axis.text.y = element_text(size = 6, hjust = 0), axis.title.x = element_text(size = 20))
yaxis <- get_y_axis(plot1) # , scale = 0.9
plot1 <- plot1 + theme(axis.text.y = element_blank()) +theme(legend.position = "None", axis.text.y = element_blank(), text = element_text(size = 8),axis.ticks = element_blank() ) + scale_fill_manual(values = colors)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                                                                                                               panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=8), plot.margin = margin(0,0,0,20, "pt"))
repeat_plot2<-ggdraw(plot1) + draw_grob(yaxis, x = 1.2, y = 1.08, hjust = 1, vjust = 1, scale = 0.75) #


#################################### 
############### Lies ############### 
#################################### 
data <- read.table("../data/data_for_final_plot.tsv", sep = "", header = TRUE)
te_families <- read.table("../data/mm39_families.tsv") %>% setNames(c("te_name", "family", "count", "mean size"))

theme_scale = scale_fill_gradient2(low = "#0053a1",  mid = "#daebf7", high = "#0053a1", midpoint = 0)
plotted_col = "mean_expr"
filter = TRUE 

data$famSub = te_families[match(data$te_name,
                                te_families$te_name), "family"]
tmp = str_split(data$famSub, "/", simplify = TRUE)
data$fam = tmp[, 1]
data$subfam = tmp[, 2]
data[is.na(data$famSub), "famSub"] = "merged"
data[is.na(data$subfam), "subfam"] = "merged"
data[is.na(data$fam), "fam"] = "merged"
data[data$subfam == "", "subfam"] = data[data$subfam == "", "famSub"]#

if (filter) { # shall we filter families ?
 keep_fam = c("SINE", "LINE", "LTR", "DNA")
  data = data[data$fam %in% keep_fam, ]
}

statuses = c("Annotated", "Novel-annotated", "Novel-novel")
data$status = statuses[1]
data[!grepl("ENS", data$matching_transcript) &
       !grepl("ENS", data$matching_gene), "status"] = statuses[3]
data[!grepl("ENS", data$matching_transcript) &
       grepl("ENS", data$matching_gene), "status"] = statuses[2]
statuses = statuses[match(unique(data$status), statuses)]
statuses = statuses[order(statuses)] #

data$matching_transcript = factor(x = data$matching_transcript,
         levels = unique(data[order(data$fam,
                      data[,plotted_col],
                      data$trcov),"matching_transcript"]))

widths = sapply(statuses, function(x) {
  nrow(data[data$status == x, ])
})  / nrow(data)

widths[3] <- 0.1#

discrete_factor = 1 / 2
data$disc_trcov = as.factor(round(data$trcov / discrete_factor)
                            * discrete_factor)
colors = c("grey90", "grey50", "black")
names(colors) = c(0, 0.5, 1)

labels = c("Low", "Half-half", "High")
names(labels) = c(0, 0.5, 1)
labels = as_labeller(labels)


relheights <- c(nrow(filter(data, status == 'Annotated')),nrow(filter(data, status == 'Novel-annotated')),nrow(filter(data, status == 'Novel-novel'))) %>% log2()
fam_TE_plots <- plot_grid(ggplot(filter(data, status == 'Annotated'), aes(y= matching_transcript, x = fam, fill = disc_trcov)) + 
                            geom_tile() + 
                            scale_fill_discrete(type = colors, labels = labels, drop = FALSE) +
                            xlab(NULL) + 
                            ylab(NULL) +
                            theme_classic() +
                            theme(axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  panel.grid = element_blank(),
                                  legend.position = "none",
                                  text = element_text(size = base_size), 
                                  plot.margin = margin(0,0,5,0, "pt")) + 
                            scale_y_discrete(position="left", limits = rev) +
                            scale_x_discrete(na.translate = FALSE, drop = FALSE, position = "top", limits = rev) +
                            guides(fill = guide_legend(title = "Transcript TE content")), 
                          ggplot(filter(data, status == 'Novel-annotated'), aes(y= matching_transcript, x = fam, fill = disc_trcov)) + 
                            geom_tile() + 
                            scale_fill_discrete(type = colors, labels = labels, drop = FALSE) +
                            xlab(NULL) + 
                            ylab(NULL) +
                            theme_classic() +
                            theme(axis.text.x = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  axis.ticks.y = element_blank(),
                                  panel.grid = element_blank(),
                                  legend.position = "none",
                                  text = element_text(size = base_size),
                                  plot.margin = margin(0,0,5,0, "pt")) + 
                            scale_y_discrete(position="left", limits = rev) +
                            scale_x_discrete(na.translate = FALSE, drop = FALSE, position = "top", limits = rev) +
                            guides(fill = guide_legend(title = "Transcript TE content")), 
                          ggplot(filter(data, status == 'Novel-novel'), aes(y= matching_transcript, x = fam, fill = disc_trcov)) + 
                            geom_tile() + 
                            scale_fill_discrete(type = colors, labels = labels, drop = FALSE) +
                            xlab(NULL) + 
                            ylab(NULL) +
                            theme_classic() +
                            theme(axis.text.x = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  axis.ticks.y = element_blank(),
                                  panel.grid = element_blank(),
                                  legend.position = "none",
                                  text = element_text(size = base_size),
                                  plot.margin = margin(0,0,5,0, "pt")) + 
                            scale_y_discrete(position="left", limits = rev) +
                            scale_x_discrete(na.translate = FALSE, drop = FALSE, position = "top", limits = rev) +
                            guides(fill = guide_legend(title = "Transcript TE content")), ncol=1, rel_heights = relheights)



expr_TE_plots <- plot_grid(ggplot(filter(data, status == 'Annotated'), aes(y = matching_transcript, x = 1, fill = mean_expr)) +
                             geom_tile(width=50) +
                             theme_classic() +
                             scale_fill_gradient2(low = "#66c2a4",  mid = "#f8fcfb", high = "#66c2a4", midpoint = 0) +
                             xlab("") + ylab(NULL)  +
                             theme(axis.text.y = element_blank(),
                                   axis.ticks = element_blank(),
                                   panel.grid = element_blank(),
                                   legend.position = "none",
                                   plot.margin = margin(0,0,5,20, "pt"),
                                   text = element_text(size = base_size),
                                   axis.text.x = element_blank(),
                                   axis.line.y = element_blank()) + scale_y_discrete(limits = rev) +
                             scale_x_discrete(position = "top"),
                           ggplot(filter(data, status == 'Novel-annotated'), aes(y = matching_transcript, x = 1, fill = mean_expr)) +
                             geom_tile(width=50) +
                             theme_classic() +
                             scale_fill_gradient2(low = "#fc8e62",  mid = "#fffaf8", high = "#fc8e62", midpoint = 0) +
                             xlab(NULL) + ylab(NULL)  +
                             theme(axis.text.y = element_blank(),
                                   axis.ticks = element_blank(),
                                   panel.grid = element_blank(),
                                   legend.position = "none",
                                   plot.margin = margin(0,0,5,20, "pt"),
                                   text = element_text(size = base_size),
                                   axis.text.x = element_blank(),
                                   axis.line.y = element_blank()) + scale_y_discrete(limits = rev) +
                            scale_x_discrete(position = "top"),
                           ggplot(filter(data, status == 'Novel-novel'), aes(y = matching_transcript, x = 1, fill = mean_expr)) +
                             geom_tile(width=50) +
                             theme_classic() +
                             scale_fill_gradient2(low = "#0053a1",  mid = "#daebf7", high = "#0053a1", midpoint = 0) +
                             xlab(NULL) + ylab(NULL)  +
                             theme(axis.text.y = element_blank(),
                                   axis.ticks = element_blank(),
                                   panel.grid = element_blank(),
                                   legend.position = "none",
                                   plot.margin = margin(0,0,5,20, "pt"),
                                   text = element_text(size = base_size),
                                   axis.text.x = element_blank(),
                                   axis.line.y = element_blank()) + scale_y_discrete(limits = rev) +
                             scale_x_discrete(position = "top") , ncol=1, rel_heights = relheights)


leg <- get_legend(ggplot(filter(data, status == 'Annotated'), aes(y = matching_transcript, x = 1, fill = mean_expr)) +
                    geom_tile(width=50) +
                    theme_classic() +
                    scale_fill_gradient2(low = "#fc8e62",  mid = "#fffaf8", high = "#fc8e62", midpoint = 0) +
                    xlab(NULL) + ylab(NULL)  +
                    theme(legend.key.size = unit(8, 'pt'),
                          axis.text.y = element_blank(),
                          axis.ticks = element_blank(),
                          panel.grid = element_blank(),
                          legend.position = "right",
                          text = element_text(size = base_size),
                          axis.text.x = element_blank(),
                          axis.line.y = element_blank()) + scale_y_discrete(limits = rev) +
                    scale_x_discrete(position = "top") +
                    labs(fill = "expression"))#

leg2 <- get_legend(ggplot(filter(data, status == 'Novel-novel'), aes(y= matching_transcript, x = fam, fill = disc_trcov)) + 
                     geom_tile() + 
                     scale_fill_discrete(type = colors, labels = labels, drop = FALSE) +
                     xlab(NULL) + 
                     ylab(NULL) +
                     theme_classic() +
                     theme(legend.key.size = unit(8, 'pt'),
                           axis.text.x = element_blank(),
                           axis.text.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           panel.grid = element_blank(),
                           legend.position = "top",
                           text = element_text(size = base_size)) + 
                     scale_y_discrete(position="right", limits = rev) +
                     scale_x_discrete(na.translate = FALSE, drop = FALSE, position = "top") +
                     guides(fill = guide_legend(title = "Transcript \nTE content")))


TE_heatmap <- plot_grid(plot_grid(expr_TE_plots, fam_TE_plots, as_ggplot(leg), nrow=1, rel_widths = c(1.5,2,1.5)), rel_heights = c(10,1), as_ggplot(leg2), ncol=1)


left_join(as.data.frame(abundance_no_batch_effect) %>%  mutate(transcript_id = row.names(.)), annotation, by='transcript_id') %>% 
  write.table(., paste(list_dir, "/Fig1_abundance_no_batch_effect.txt", sep=""), row.names = FALSE, quote = FALSE, sep = "\t")  #save the list as txt to read later as output of this figure



###

ggsave(paste(fig_dir,"/1_productivity.svg",sep=""), plot=productivity_plot, width=7, height = 2, create.dir = TRUE)
ggsave(paste(fig_dir,"/1_numofreads_perclass.svg",sep=""), plot=nb_of_reads_per_class_plot, width=2, height = 2)
ggsave(paste(fig_dir,"/1_no_of_exons.svg",sep=""), plot=nb_exons_plot, width=3, height = 3)
ggsave(paste(fig_dir,"/1_transctipt_lengthdensity.svg",sep=""), plot=tx_length_plot, width=3, height = 3)
ggsave(paste(fig_dir,"/1_abundance_by_classs.svg",sep=""), plot=abundance_by_class_plot, width=3, height = 3)
ggsave(paste(fig_dir,"/1_upset_plot_grob_mixed.svg",sep=""), plot=upset_plot_grob, width=7, height = 4)

###
ggsave(paste(fig_dir,"/1_repeatplot.svg",sep=""), plot=repeat_plot1, width=3, height = 3.5)
ggsave(paste(fig_dir,"/1_repeat_plot2.svg",sep=""), plot=repeat_plot2, width=2.5, height = 2)
ggsave(paste(fig_dir,"/1_heatmapTE.svg",sep=""), plot=TE_heatmap, width=2, height = 3)


suppa_data2 <- imap(suppa_files, function(f, x) {
  read.table(f, header=TRUE) %>% rowwise() %>%
    mutate(total_transcripts = strsplit(as.character(total_transcripts), ","),
           alternative_transcripts = strsplit(as.character(alternative_transcripts), ","),
           complement = list(setdiff(unlist(total_transcripts), unlist(alternative_transcripts))),
           in_tx_list = any(unlist(complement) %in% expr_tx_list)  |   any(unlist(alternative_transcripts) %in% expr_tx_list),
           novel_event = ((length(complement)==1) & all(!grepl("^ENSMUST", unlist(complement)))) | ((length(alternative_transcripts) == 1) & all(!grepl("^ENSMUST", unlist(alternative_transcripts))))) %>% 
    filter(in_tx_list) %>%
    mutate(event_type = x,
           list_transcripts = case_when(
             event_type == "MX" ~ list(total_transcripts),
             event_type == "ES" ~ list(complement),
             TRUE ~ list(alternative_transcripts)
           )) %>%
    dplyr::select(-alternative_transcripts, -complement, -total_transcripts) %>% 
    unnest(cols=c(list_transcripts)) %>%
    dplyr::rename(tx_id = list_transcripts) %>% 
    mutate(event_type = x)}) %>%
  bind_rows() 
#%>%
#dplyr::select(-seqname) %>%
#filter(tx_id %in% expr_tx_list)
summarise(suppa_data2, total_genes = n_distinct(gene_id), total_tx = n_distinct(tx_id), total_events = n_distinct(event_id))
suppa_data2 %>% group_by(event_type) %>% summarize(event_by_type = n_distinct(event_id))
upset_dataX <- suppa_data2 %>% dplyr::select(-event_id) %>%  
  mutate(value=1) %>% 
  pivot_wider(id_cols = tx_id, names_from = event_type, values_from = value, values_fn = function(x) as.integer(length(x)>0), values_fill = 0)

suppa_data2$event_type <- factor(suppa_data2$event_type, levels = colnames(upset_dataX[,-1])[order(colSums(upset_dataX[,-1]))])
upset_dataX %>% mutate(n_tx = rowSums(across(where(is.numeric))))
Upset_plot_larger <- upset(as.data.frame(upset_dataX), main.bar.color = "grey30", sets.x.label = "Total number of isoforms", 
                           mainbar.y.label = "Number of isoforms", text.scale = 1.1, 
                           query.legend = "top", nsets = 7, sets.bar.color = "grey50", empty.intersections = FALSE,
                           show.numbers = "yes",number.angles = 30, point.size = 2.5, line.size = 0.9, 
                           queries = list(list(query = function(row) grepl("^ENSMUST", row["tx_id"]), 
                                               color = "#66c2a5", active = T)))

ggsave(paste(fig_dir,"/1_upsetlarger_one_transcript_expressed.svg",sep=""), plot=wrap_ggplot_grob(Upset_plot_larger$Main_bar)/ wrap_ggplot_grob(Upset_plot_larger$Matrix), width=7, height = 4)
#----------------change the number of intersections with you want all or not--------###
#-------------------------saving only singlets and doublets------# not larger_normal ##
upSets<-c("AF","AL", "AS3","AS5","MX","RI","ES")
upSets <- factor(upSets, levels=upSets)
#suppa data by taking liver-spec here we just change ;like  instead of & we put OR
#to get if atleast one transcript is present for an AF even
#Upset_plot_single<- upset(as.data.frame(upset_data),order.by="degree",decreasing=F, sets = c("AF","AL" ,"AS3","AS5","RI","MX","ES"), main.bar.color = "black", sets.x.label = "Total number of isoforms", 
#                    mainbar.y.label = "No. of isoforms", text.scale = 1.1,nintersects=7,matrix.color = "black",
#                    query.legend = "top", nsets = 7, sets.bar.color = "black", empty.intersections = TRUE,
#                    show.numbers = "no", point.size = 2.5, line.size = 0.9,
#                    queries = list(list(query = function(row) grepl("^ENSMUST", row["tx_id"]), 
#                    color = "#66c2a5", active = T)), keep.order=T )
#ggsave(paste(fig_dir,"/1_upset_single_events.svg",sep=""), plot=wrap_ggplot_grob(Upset_plot_single$Main_bar)/ wrap_ggplot_grob(Upset_plot_single$Matrix), width=2, height = 3)#

write.table(upset_data, paste(list_dir, "/1_upset.tsv", sep=""), row.names = FALSE, quote = FALSE, sep = "\t")  #save the list as txt to read later as output of this figure


Upset_plot_single_double<- upset(as.data.frame(upset_data),order.by="degree",decreasing=F, sets = c("AF","AL" ,"AS3","AS5","RI","MX","ES"), main.bar.color = "black", sets.x.label = "Total number of isoforms", 
                    mainbar.y.label = "No. of isoforms", text.scale = 1.1,nintersects=28,matrix.color = "black",
                    query.legend = "top", nsets = 7, sets.bar.color = "black", empty.intersections = TRUE,
                    show.numbers = "yes", point.size = 2.5, line.size = 0.9,
                    queries = list(list(query = function(row) grepl("^ENSMUST", row["tx_id"]), 
                    color = "#66c2a5", active = T)), keep.order=T )#

ggsave(paste(fig_dir,"/1_upset_single_double.svg",sep=""), plot=wrap_ggplot_grob(Upset_plot_single_double$Main_bar)/ wrap_ggplot_grob(Upset_plot_single_double$Matrix), width=6, height = 3.5)

#----------------------------------saving third and fourth------------------------------#
#triplet intersections
upset_data_trip<- upset_data %>% filter(rowSums(upset_data[2:8]) == 3 )
Upset_plot_triple <- upset(as.data.frame(upset_data_trip),
                           order.by="degree",
                           decreasing=F, sets = c("AF","AL" ,"AS3","AS5","RI","MX","ES"), 
                           main.bar.color = "black", sets.x.label = "Total number of isoforms", 
                           mainbar.y.label = "Isoforms", 
                           text.scale = 1, matrix.color = "black",
                           query.legend = "top", nsets = 7, 
                           sets.bar.color = "black",
                          show.numbers = "yes", point.size = 0.8, line.size = 0.3, 
                           queries = list(list(query = function(row) grepl("^ENSMUST", row["tx_id"]), 
                                               color = "#66c2a5", active = T)), keep.order=T)
ggsave(paste(fig_dir,"/1_upset_triple.svg",sep=""), plot=wrap_ggplot_grob(Upset_plot_triple$Main_bar)/ wrap_ggplot_grob(Upset_plot_triple$Matrix), width=3.7, height = 2.5)

#doublet intersections
#upset_data_double<- upset_data %>% filter(rowSums(upset_data[2:8]) == 2 )
#Upset_plot_double <- upset(as.data.frame(upset_data_double),order.by="degree",decreasing=F, sets = c("AF","AL" ,"AS3","AS5","RI","MX","ES"), main.bar.color = "black", sets.x.label = "Total number of isoforms", 
#                           mainbar.y.label = "Isoforms", text.scale = 1, matrix.color = "black",
#                           query.legend = "top", nsets = 7, 
#                           sets.bar.color = "black",
#                           show.numbers = "no", point.size = 0.8, line.size = 0.3, 
#                           queries = list(list(query = function(row) grepl("^ENSMUST", row["tx_id"]), 
#                                               color = "#66c2a5", active = T)), keep.order=T)
#
#ggsave(paste(fig_dir,"/Fig1_upset_double.svg",sep=""), plot=wrap_ggplot_grob(Upset_plot_double$Main_bar)/ wrap_ggplot_grob(Upset_plot_double$Matrix), width=3, height = 1.5)#

# double and triple
#upset_data__triple<- upset_data %>% filter(rowSums(upset_data[2:8]) > 1 & rowSums(upset_data[2:8]) < 4)
#Upset_plot_triple <- upset(as.data.frame(upset_data_double_triple),
#                           order.by="degree",decreasing=F, 
#                           sets = c("AF","AL" ,"AS3","AS5","RI","MX","ES"), 
#                           main.bar.color = "black",
#                           sets.x.label = "Total number of isoforms", #
#                           mainbar.y.label = "Isoforms",
#                           text.scale = 1, 
 #                          matrix.color = "black",
#                           query.legend = "top",
#                           nsets = 7, 
#                           sets.bar.color = "black",  
#                           show.numbers = "no",
#                           point.size = 0.8, 
#                           line.size = 0.3, 
 #                          queries = list(list(query = function(row) grepl("^ENSMUST", row["tx_id"]), 
#                                               color = "#66c2a5", active = T)), keep.order=T)

#ggsave(paste(fig_dir,"/Fig1_upset_double_triple.svg",sep=""), plot=wrap_ggplot_grob(Upset_plot_double_triple$Main_bar)/ wrap_ggplot_grob(Upset_plot_double_triple$Matrix), width=3, height = 1.5)

#double triple without AF AL
#upset_data_double_triple_noAF<- upset_data %>% select(-c("AF", "AL")) %>% filter(rowSums(.[2:6]) > 1 & rowSums(.[2:6]) < 4) 
#upset_data_double_triple_noAF <- upset(as.data.frame(upset_data_double_triple_noAF),order.by="degree",decreasing=F, sets = c("AS3","AS5","RI","MX","ES"), main.bar.color = "black", sets.x.label = "Total number of isoforms", 
#                           mainbar.y.label = "Isoforms", text.scale = 1, matrix.color = "black",
#                           query.legend = "top", nsets = 7, 
#                           sets.bar.color = "black",
#                           show.numbers = "no", point.size = 0.8, line.size = 0.3, 
#                           queries = list(list(query = function(row) grepl("^ENSMUST", row["tx_id"]), 
#                                               color = "#66c2a5", active = T)), keep.order=T)

#ggsave(paste(fig_dir,"/Fig1_upset_double_triple_noAFAL.svg",sep=""), plot=wrap_ggplot_grob(upset_data_double_triple_noAF$Main_bar)/ wrap_ggplot_grob(upset_data_double_triple_noAF$Matrix), width=3, height = 1.5)

###################################################
#f1 <- plot_grid(NULL, plot_grid(nb_of_reads_per_class_plot, nb_exons_plot, tx_length_plot, abundance_by_class_plot, nrow=1, labels = c('B', 'C', 'D', 'E')), plot_grid(plot_grid(TE_heatmap, NULL, ncol = 1, rel_heights = c(10,1)), productivity_plot, nrow=1, labels=c('', 'G'), rel_widths = c(1,3)), plot_grid(plot_grid(repeat_plot1,repeat_plot2, nrow=2, rel_heights = c(1.2,1), labels = c('H', 'I')), upset_plot_grob, nrow=1, rel_widths = c(1,4), labels = c('', 'J')), rel_heights = c(0.75,0.6,1,1.75), ncol=1, labels = c('A', '', 'F'))
#f1
#ggsave("../results/Figure1.pdf", plot=f1, width=9, height = 11)

#extra_plot<- ggplot() + theme_void()
#layout <- c(area(1,1,2,2), area(1,2,1,3),area(1,5), area(2,3,2,5), area(3, 1, 3, 5), area(4,1,5,1), area(4,2,5,5))
#f1 <- extra_plot+ extra_plot+wrap_elements(panel=nb_of_reads_per_class_plot) + (nb_exons_plot | tx_length_plot | abundance_by_class_plot) +wrap_elements(plot = productivity_plot) + plot_grid(repeat_plot1,repeat_plot2, nrow=2, rel_heights = c(1.2,1)) + upset_plot_grob  + 
#  plot_layout(design = layout, widths = c(1,0.25,1,1,1), heights = c(1,0.8,1,0.8,1.0), ncol=5, nrow=5)+ plot_annotation(tag_levels = 'A' )
#ggsave("../results/Figure1.pdf", plot=f1, width=9.5, height = 12.5)
#ggsave("../results/Figure1.png", plot=f1, width=9.5, height = 12.5)
#ggsave("../extra/fig1.1.png", nb_of_reads_per_class_plot, width=9.5, height = 12.5)
#ggsave("../extra/fig1.2.png", nb_exons_plot,width=9.5, height = 12.5)
#ggsave("../extra/fig1.3.png", abundance_by_class_plot,width=9.5, height = 12.5)
#ggsave("../extra/fig1.4.png", tx_length_plot,width=9.5, height = 12.5)
#ggsave("../extra/fig1.5.png", productivity_plot,width=9.5, height = 12.5)
#ggsave("../extra/fig1.6.png", repeat_plot1,width=9.5, height = 12.5)
#ggsave("../extra/fig1.7.png", repeat_plot2,width=9.5, height = 12.5)
#ggsave("../extra/fig1.8.png", upset_plot_grob,width=9.5, height = 12.5)
###
#save all figures separately
#ggsave(paste(fig_dir,"/Fig1B.pdf",sep=""), plot=productivity_plot, width=7, height = 3)
#ggsave(paste(fig_dir,"/Fig1C.pdf",sep=""), plot=repeat_plot1, width=3, height = 4)
#ggsave(paste(fig_dir,"/FigS1-1.pdf",sep=""), plot=nb_of_reads_per_class_plot, width=2, height = 2)
#ggsave(paste(fig_dir,"/FigS1-2.pdf",sep=""), plot=nb_exons_plot, width=2.5, height = 2)
#ggsave(paste(fig_dir,"/FigS1-3.pdf",sep=""), plot=tx_length_plot, width=2, height = 2)
#ggsave(paste(fig_dir,"/FigS1-4.pdf",sep=""), plot=abundance_by_class_plot, width=2, height = 2)
#ggsave(paste(fig_dir,"/FigS1-5.pdf",sep=""), plot=repeat_plot2, width=2, height = 3)
#ggsave(paste(fig_dir,"/FigS1-6.pdf",sep=""), plot=TE_heatmap, width=2, height = 3)
#ggsave(paste(fig_dir,"/Fig2B-C.pdf",sep=""), plot=upset_plot_grob, width=7, height = 4)