#set working directory and load required packages
setwd("C:/Users/kmamgain/prettyCoolPores_new/scripts") 

library(patchwork)
library(ggplot2)
library(clusterProfiler)
library(tidyverse)
library(betareg)
library(multidplyr)
library(lmtest)
library(stageR)
library(msigdbr)
library(circular)
library(rtracklayer)
library(GenomicFeatures)
library(ggthemes)
library(ggtranscript)
library(biomaRt)
library(cowplot)
library(ggpubr)
library(dplyr)
library(tximport)
library(DESeq2)

event_order <- c("AF", "AL", "A3", "A5", "RI", "MX", "SE")

#define cutoffs
fdr_cutoff <- 0.05
base_size <- 9
retrieve_from_biomart<-FALSE

#path to the directory where we want to save our figures
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
#load expression list
expressed_t<-read.csv("../results/expr_tx_list.csv")

#-----------functions--------------------------------#

### round isoform structure rectangles

geom_rangeX <- function(mapping = NULL, data = NULL,
                        stat = "identity", position = "identity",
                        ...,
                        vjust = NULL,
                        linejoin = "mitre",
                        na.rm = FALSE,
                        radius = grid::unit(1, "pt"), #change here
                        show.legend = NA,
                        inherit.aes = TRUE) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomRange,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      vjust = vjust,
      linejoin = linejoin,
      radius = radius, # K added radius here
      na.rm = na.rm,
      
      ...
    )
  )
}
GeomRange <- ggplot2::ggproto("GeomRange", ggplot2::GeomTile,
                              required_aes = c("xstart", "xend", "y"),
                              default_aes = aes(
                                fill = "grey",
                                colour = "black",
                                size = 0.25,
                                linetype = 1,
                                alpha = NA,
                                height = NA
                              ),
                              setup_data = function(data, params) {
                                # modified from ggplot2::GeomTile
                                data$height <- data$height %||% params$height %||% 0.5
                                
                                transform(
                                  data,
                                  xmin = xstart,
                                  xmax = xend,
                                  ymin = y - height / 2,
                                  ymax = y + height / 2,
                                  height = NULL
                                )
                              },
                              draw_panel = function(self,
                                                    data,
                                                    panel_params,
                                                    coord,   radius = grid::unit(6, "pt"),
                                                    vjust = NULL,
                                                    lineend = "butt",
                                                    linejoin = "mitre") {
                                if (!coord$is_linear()) {
                                  # prefer to match geom_curve and warn
                                  # rather than copy the implementation from GeomRect for simplicity
                                  # also don'think geom_range would be used for non-linear coords
                                  warn("geom_ is not implemented for non-linear coordinates")
                                }
                                
                                
                                
                                coords <- coord$transform(data, panel_params)
                                # K added part here inspired from https://stackoverflow.com/questions/64355877/round-corners-in-ggplots-geom-tile-possible
                                lapply(1:length(coords$xmin), function(i) {      
                                  
                                  grid::roundrectGrob(
                                    coords$xmin[i], coords$ymax[i],
                                    width = (coords$xmax[i] - coords$xmin[i]),
                                    height = (coords$ymax[i] - coords$ymin)[i],
                                    r = radius,
                                    default.units = "native",
                                    just = c("left", "top"),
                                    gp = grid::gpar(
                                      col = coords$colour[i],
                                      fill = alpha(coords$fill[i], coords$alpha[i]),
                                      lwd = coords$size[i] * .pt,
                                      lty = coords$linetype[i],
                                      lineend = "butt"
                                    )
                                  )
                                  
                                }) -> gl
                                
                                grobs <- do.call(grid::gList, gl)
                                
                                ggname("geom_rrect", grid::grobTree(children = grobs))
                                
                              }
                              
)

ggname <- function(prefix, grob) {
  grob$name <- grid::grobName(grob, prefix)
  grob
}

library(ggplot2)

###
beta_omnibus_fit <- function(data, ...){
  if (any(data$psi == 0) || any(data$psi == 1)) {
    n <- nrow(data)
    data$psi <- (data$psi * (n - 1) + 0.5)/n
  }
  success <- try({
    if (length(unique(data$eventID))>1){
      fit1 <- betareg(psi~eventID*(inphase + eventID:outphase), data = data, ...)
      fit0 <- betareg(psi~eventID, data = data, ...)
    } else {
      fit1 <- betareg(psi~inphase + outphase, data = data, ...)
      fit0 <- betareg(psi~1, data = data, ...)
    }
    1
  })
  if (success == 1){
    p_val <- lrtest(fit1, fit0)[5][2,]
  }
  else {
    p_val <- NA
  }
  names(p_val) <- "p_value_omni"
  return(p_val)
}

quasibinom_fit <- function(data, ...){
  success <- try({
    fit1 <- glm(psi~inphase + outphase, data = data, family = quasibinomial(link = "logit"))
    fit0 <- glm(psi~1, data = data, family = quasibinomial(link = "logit"))
    1
  })
  if (success == 1){
    coefs <- coef(fit1)
    p_val <- anova(fit1, fit0, test = "LRT")[2,5]
  }
  else {
    coefs <- rep(NA, 4)
    p_val <- NA
  }
  names(p_val) <- "p_value"
  return(tibble::as_tibble_row(c(coefs, p_val)))
}

quasibinom_omnibus_fit <- function(data, ...){
  success <- try({
    if (length(unique(data$eventID))>1){
      fit1 <- glm(psi~eventID*(inphase + outphase), data = data, family = quasibinomial(link = "logit"))
      fit0 <- glm(psi~eventID, data = data, family = quasibinomial(link = "logit"))
    } else {
      fit1 <- glm(psi~inphase + outphase, data = data, family = quasibinomial(link = "logit"))
      fit0 <- glm(psi~1, data = data, family = quasibinomial(link = "logit"))
    }
    1
  })
  if (success == 1){
    p_val <- anova(fit1, fit0, test = "LRT")[2,5]
  }
  else {
    p_val <- NA
  }
  names(p_val) <- "p_value_omni"
  return(p_val)
}

########
txdb <- makeTxDbFromGFF(paste(flair_path,"/transcriptome_productivity.gtf", sep=""), format = "gtf")
annotation <- transcriptLengths(txdb, with.cds_len = TRUE, 
                                with.utr3_len = TRUE, with.utr5_len = TRUE) %>%
  dplyr::select(-tx_id) %>%
  dplyr::rename(transcript_id="tx_name")

if (retrieve_from_biomart) {
  ensembl <- useMart("ensembl", "mmusculus_gene_ensembl", host = "https://oct2022.archive.ensembl.org")
  filters <- listFilters(ensembl)
  attribs <- listAttributes(ensembl)
  mart <- getBM(filters = "ensembl_transcript_id",
                attributes = c("ensembl_transcript_id", "ensembl_gene_id", "mgi_symbol", "entrezgene_id" ,"transcript_biotype"),
                values = annotation$transcript_id,
                mart = ensembl)
  colnames(mart) <- sub("^ensembl_","", colnames(mart))
  mart %<>% distinct(transcript_id, .keep_all = TRUE)
  write_rds(mart, paste0("../data/biomart_annotation_", format(Sys.time(), "%Y-%m-%d.rds")), compress="bz2")
  } else {
  mart <- readRDS(biomart_path)
  }
  annotation %<>% left_join(mart) %>%
  mutate(class=case_when(grepl("^ENS", transcript_id) ~ "Annotated", 
                         grepl("ENS", gene_id) ~ "Novel-Annotated", 
                         .default = "Novel-Novel")) %>%
  group_by(gene_id) %>% 
  mutate(mgi_symbol = ifelse(grepl("^ENS", gene_id), 
                             unique(na.omit(mgi_symbol)), mgi_symbol)) %>%
  as.data.frame
  df<- data.frame(gene_id= annotation$gene_id, mgi_symbol=annotation$mgi_symbol) %>% filter(!is.na(mgi_symbol))
  annotation$mgi_symbol<- df$mgi_symbol[match(annotation$gene_id, df$gene_id)]

  gene_annot <- annotation %>%  dplyr::rename(geneID=gene_id) %>% dplyr::select(geneID, mgi_symbol)  %>%
    distinct(geneID, .keep_all = TRUE) 

  suppa_files <- c(AS3 = paste(suppa_path, "/transcriptome_ext_A3_strict.ioe", sep=""),
                   AS5 = paste(suppa_path, "/transcriptome_ext_A5_strict.ioe", sep=""),
                   AF = paste(suppa_path, "/transcriptome_ext_AF_strict.ioe", sep=""),
                   AL =paste(suppa_path, "/transcriptome_ext_AL_strict.ioe", sep=""),
                   MX = paste(suppa_path, "/transcriptome_ext_MX_strict.ioe", sep=""),
                   RI = paste(suppa_path, "/transcriptome_ext_RI_strict.ioe", sep=""),
                   ES = paste(suppa_path, "/transcriptome_ext_SE_strict.ioe", sep=""))
  
  
  suppa_annot <- map(suppa_files, function(f) read.table(f, header=TRUE)) %>% 
    bind_rows %>% dplyr::select(-seqname, -gene_id) %>% 
    separate_wider_delim(event_id, ";", names = c("geneID","eventID"))
  
#psi files
files <- list.files(path = paste(suppa_path, "/", sep=""), pattern = ".psi", full.names = TRUE, include.dirs = TRUE, no.. = TRUE)
psi_data <- map(files, function(f) read.delim(f, row.names = 1) %>% rownames_to_column(var = "eventID")) %>%
  bind_rows %>% 
  filter(if_all(starts_with("CT"), ~ !is.na(.))) %>%
  pivot_longer(starts_with("CT"), names_to = "time", values_to = "psi") %>%
  mutate(time = as.numeric(str_extract(time, "^CT(\\d+)_", group=1)),
         inphase = cos(2*pi*time/24),
         outphase = sin(2*pi*time/24)) %>%
  separate_wider_delim(eventID, names = c("geneID","eventID"), delim = ";")  %>% dplyr::select(-c(3:6))
#
cluster <- new_cluster(6)
cluster_library(cluster, c("tidyverse","betareg", "lmtest"))
cluster_copy(cluster, c("quasibinom_fit", "quasibinom_omnibus_fit"))
fits <- psi_data %>%
  nest_by(geneID, eventID) %>%
  filter(max(data$psi)>0.0001 & min(data$psi)<0.999) %>%
  partition(cluster) %>%
  mutate(fit = list(quasibinom_fit(data))) %>%
  collect() %>%
  unnest(c(data)) %>%
  ungroup %>%
  nest_by(geneID) %>%
  partition(cluster) %>%
  mutate(p_value_omni = quasibinom_omnibus_fit(data)) %>%
  collect() %>%
  unnest(c(data)) %>%
  group_by(geneID, eventID) %>% 
  slice_head(n=1) %>% 
  dplyr::select(-time, -inphase, -outphase, -psi) %>% 
  unnest(c(fit)) %>%
  mutate(type = str_extract(eventID, "\\w{2}"), 
         geneID = gsub(":","_",geneID), 
         eventID=gsub(":","_",eventID),
         phase = (atan2(outphase, inphase) %% (2*pi)) * 12/pi) %>%
  dplyr::select(-inphase, -outphase, -`(Intercept)`) %>%
  filter(!is.na(p_value) & !is.na(p_value_omni))

test <- psi_data %>% nest_by(geneID) %>% {.[1:500,]} %>% unnest(c(data)) %>% ungroup

pScreen <- fits$p_value_omni
names(pScreen) <- fits$geneID
pConfirmation <- matrix(fits$p_value, ncol = 1)
dimnames(pConfirmation) <- list(fits$eventID,c("transcript"))
event2gene <- data.frame(fits$eventID, fits$geneID)
stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted= FALSE, tx2gene=event2gene)
stageRObj <- stageWiseAdjustment(object=stageRObj, method="dte", alpha=fdr_cutoff)
rhy_event_analysis <- getAdjustedPValues(stageRObj, order=TRUE, onlySignificantGenes=FALSE) %>%
                  mutate(type = str_extract(txID, "\\w{2}")) %>%
                  dplyr::rename(eventID = txID, event=transcript) %>%
                  left_join(fits) %>%
                  dplyr::select(-p_value, -p_value_omni) %>%
                  left_join(gene_annot) %>%
                  mutate(eventID = gsub("_", ":", eventID)) %>%
                  left_join(suppa_annot)  %>%
                  rowwise() %>% 
                  mutate(alternative_transcripts = str_split(alternative_transcripts, ","), 
                        total_transcripts = str_split(total_transcripts, ","), 
                        primary_transcripts = list(setdiff(total_transcripts, alternative_transcripts)), 
                        novel = !any(grepl("^ENSMUST", alternative_transcripts)) & !any(grepl("^ENSMUST", primary_transcripts)),
                        n= ifelse(novel=="TRUE", "novel-event", "annotated-event"))

write_rds(rhy_event_analysis[,c(1:12)], "../data/rhy_event_analysis_all.rds")

anno_bed<-read.table(paste(flair_path,"/productivity.txt", sep=""), sep="\t") %>% 
  mutate(transcript_id= .$V4, productivity_flair= .$V13) %>% left_join(., annotation) %>% 
  mutate(productivity_flair = fct_relevel(productivity_flair, c("PRO","PTC","NGO","NST")),
         productivity_flair = fct_recode(productivity_flair, 
                                         !!!c("productive" = "PRO",
                                              "premature termination codon" = "PTC",
                                              "no start codon" = "NGO",
                                              "has start but no stop codon" = "NST"))) %>% filter(transcript_id %in% rownames(expressed_t)) %>% dplyr::select(c(1:13,16))
colnames(anno_bed) <- c("chrom", "chromStart", "chromEnd", "tx_name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "exons_nb", "exons_size", "exons_start",  "productivity","gene_name")

m_df <- msigdbr(species = "Mus musculus", category = "C3", subcategory = "MIR:MIRDB")
m_t2g_C2 <- m_df %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
msig_C2 <- clusterProfiler::enricher(rhy_event_analysis %>% filter(gene< fdr_cutoff) %>% pull(mgi_symbol) %>% unique(),
                     universe = rhy_event_analysis %>% pull(mgi_symbol) %>% unique(),
                     TERM2GENE = m_t2g_C2) #, qvalueCutoff = 1, pvalueCutoff = 1)

rhy_event_analysis %>% 
  filter(event<fdr_cutoff) %>% 
  group_by(geneID) %>% 
  summarize(s = var(circular(phase*pi/12, modulo = "2pi"))) %>% 
  pull(s) %>% qplot(bins=50)

df <-  rhy_event_analysis %>% data.frame() %>%
          filter(event<fdr_cutoff) %>% 
          mutate(phase_int = (ceiling(phase/2)*2) %% 24) %>% 
          group_by(type, phase_int) %>% summarize(counts = n()) %>% 
          mutate(ncounts = counts/sum(counts))

cycling_AS<-rhy_event_analysis %>% filter(event < fdr_cutoff)
write_rds(cycling_AS[,c(1:8, 12)], "../results/cycling_AS.rds")

df_X<- rhy_event_analysis %>% filter(event < fdr_cutoff) %>% group_by(geneID) %>% summarise(count= length(unique(type))) 
df_heh<- rhy_event_analysis %>% filter(event < fdr_cutoff) %>% group_by(geneID) %>% mutate(as= type)

##--------For PSEA-----------------------##
PSEA_Plot <- ggplot(df) +
  geom_tile(aes(x=phase_int, y=factor(type, levels = event_order), fill=ncounts), color="white", linewidth=0.8) + 
  #geom_text(aes(x=24, y=type, label= totals), data = df %>% summarize(totals = sum(counts)), hjust=0.75, size=3) +
  scale_fill_gradient(low = "white", high = "black", guide = "colourbar", aesthetics = "fill", label=scales::percent) + theme_classic(base_size = base_size) + 
  #theme(aspect.ratio = 0.5,legend.justification="centre",legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-10,-10,-10)) + 
  scale_x_continuous(breaks = seq(0, 20, 4), expand = expansion(add = c(0.2,0.4)), name = "event phase (h)") +
  scale_y_discrete(expand = expansion(add=c(0.4, 0.4)), name = "") + 
  scale_x_continuous(breaks = seq(0, 20, 4), expand = expansion(add = c(0.2,0.4)), name = "event phase (h)") +
  theme(plot.margin = margin(0, 0, 0, 0),legend.position="right",
                                               legend.justification="centre",
                                               legend.margin=margin(0,0,0,0), text=element_text(size=10),
                                               legend.box.margin=margin(0,5,0,5),
        legend.title = element_blank()) 

phase_events<- ggplot(df %>% group_by(phase_int) %>% summarise(counts= sum(counts)),aes(x=phase_int, y= counts)) + 
  geom_segment(aes(x=phase_int, xend=phase_int, y=0, yend=counts), color="black",size=5) +
  geom_line(size=1, color="black") +  
  theme_bw() +
  geom_text(aes(label=counts), position=position_dodge(width=0.9), vjust=-1.2, size=3) + 
  coord_cartesian(ylim=c(0,80)) +
  theme_void() + 
  theme(legend.position="none")

num_events<- ggplot(rhy_event_analysis %>% data.frame() %>%
                      filter(event<fdr_cutoff) %>% 
                      group_by(type) %>% mutate(counts = n()),aes(x=type, y= counts, fill=n))+geom_text(aes(label= counts), vjust=0) +geom_col(position="stack", width=0.3) + theme_void()+ scale_x_discrete(expand = expansion(add=c(0.4, 0.4)), name = "")+  theme(legend.position= "bottom")+ scale_fill_manual(values= c("black","#ffbf00"))+ 
  theme(plot.margin = margin(0, 0, 0, 0), legend.title = element_blank(), legend.text = element_text(size=8), legend.key.size = unit(0.1,"cm")) 

Fig3E <- plot_grid(phase_events,PSEA_Plot, nrow=2, rel_heights = c(1,2), align = "v", axis = "lr")

##----------Plotting Isoforms Structures-------------------##
anno_bed<-read.table(paste(flair_path,"/productivity.txt", sep=""), sep="\t") %>% 
  mutate(transcript_id= .$V4, productivity_flair= .$V13) %>% left_join(., annotation) %>% 
  mutate(productivity_flair = fct_relevel(productivity_flair, c("PRO","PTC","NGO","NST")),
         productivity_flair = fct_recode(productivity_flair, 
                                         !!!c("productive" = "PRO",
                                              "premature termination codon" = "PTC",
                                              "no start codon" = "NGO",
                                              "has start but no stop codon" = "NST"))) %>% filter(transcript_id %in% expressed_t$transcript_id) %>% dplyr::select(c(1:13,16))
colnames(anno_bed) <- c("chrom", "chromStart", "chromEnd", "transcript_id", "score", "strand", "thickStart", "thickEnd", "itemRgb", "exons_nb", "exons_size", "exons_start",  "productivity","gene_name")

######a function for making psigraph####
circadian_time<- c("CT0", "CT4","CT8","CT12","CT16","CT20")
circadian_time<- factor(circadian_time,levels=c("CT0", "CT4","CT8","CT12","CT16","CT20"))
psi_dataX <- map(files, function(f) read.delim(f, row.names = 1) %>% rownames_to_column(var = "eventID")) %>%
  bind_rows %>%
  filter(if_all(starts_with("CT"), ~ !is.na(.))) %>% dplyr::select(-c(14:17)) 
psi_dataX$Amplitude <- apply(psi_dataX[, 2:13], 1, function(x) diff(range(x)))


##function to get psi graph
psi_graph<- function(events){
  data <- psi_dataX %>% data.frame( )%>% filter(eventID %in% events) %>% dplyr::select(- Amplitude )%>%
    mutate(mCT_0=rowMeans(.[2:3]),mCT_4=rowMeans(.[10:11]),mCT_8=rowMeans(.[12:13]),mCT_12=rowMeans(.[4:5]),mCT_16=rowMeans(.[6:7]),mCT_20=rowMeans(.[8:9]) ) %>% 
    dplyr::select(eventID, mCT_0, mCT_4, mCT_8, mCT_12, mCT_16, mCT_20) %>%
    setNames(c("eventID", 0,4,8,12,16,20)) %>%
    pivot_longer(!eventID, names_to = "time", values_to = "psi")
  minpsi = min(data$psi)
  maxpsi = max(data$psi)
  ggplot(data,aes(x=factor(time,levels=c(0,4,8,12,16,20)), y= psi, group=eventID))+ geom_point(size=2) + 
    geom_smooth(method = "loess", col="black", se=F) + 
    theme(panel.grid.major = element_blank(), axis.title.y = element_blank(),
          panel.grid.minor = element_blank(),panel.background = element_blank(), 
          axis.line = element_line(colour = "black"), text = element_text(size = 8), 
          axis.ticks= element_blank(), axis.title.x = element_blank()) +
    scale_y_continuous(breaks = c(round(minpsi,2),round(maxpsi,2))) + 
    scale_x_discrete(breaks=c(0,12,20)) 
}
#Function for the isoform structure with box around it
###if you do not tell about the borders, it will give borders around all the exons
gtf <- rtracklayer::import(paste(flair_path,"/transcriptome_productivity.gtf", sep="")) %>% sortSeqlevels %>% sort %>% dplyr::as_tibble() %>% filter(transcript_id %in%  anno_bed$transcript_id) %>% left_join(.,(anno_bed %>% dplyr::select(c("transcript_id", "productivity"))), by="transcript_id")
#!!!!choosing exons can get tricky
gene_structure <- function(gene_name,eventInfo,boxes, exonnum1, exonnum2,exonnum3, exonnum4, ...){
  strand<- substr(eventInfo,nchar(eventInfo),nchar(eventInfo)) 
  strand<- substr(eventInfo,nchar(eventInfo),nchar(eventInfo))
  gene_annotation_from_gtf <- gtf %>% 
    filter(gene_id == !!gene_name) %>%  
    dplyr::select(seqnames, start, end, strand, type, gene_id, transcript_id, source,productivity, gene_name)
  name<- gene_annotation_from_gtf %>% pull(gene_name)
  t_list_alternate<- unlist(rhy_event_analysis %>% filter(eventID == eventInfo) %>% pull(alternative_transcripts)) #get one alternate
  t_list_primary<-unlist(rhy_event_analysis %>% filter(eventID == eventInfo) %>% pull(primary_transcripts)) #get  primary
  t_list<- list(t_list_alternate[t_list_alternate %in% expressed_t$transcript_id][1],t_list_primary[t_list_primary %in% expressed_t$transcript_id][1] ) #choosing the first one which ever it is for each category ?
  print(t_list)
  type<- rhy_event_analysis %>% filter(eventID == eventInfo) %>% pull(type)
  gene_rescaled <- shorten_gaps(
    exons = filter(gene_annotation_from_gtf, type == "exon"), 
    introns = to_intron(filter(gene_annotation_from_gtf, type == "exon"), "transcript_id"), 
    group_var = "transcript_id") %>% filter(transcript_id %in% t_list) 
  gene_rescaled_exons <- gene_rescaled %>% dplyr::filter(type == "exon") 
  gene_rescaled_introns <- gene_rescaled %>% dplyr::filter(type == "intron") 
  ymin<- 0.5
  ymax<- 2.5
  if(boxes==1){
    if(type %in% c("SE","RI")){
      xmin1<- (gene_rescaled_exons[exonnum1,2] %>% pull(start)) - 30 
      xmax1<- (gene_rescaled_exons[exonnum1,3] %>% pull(end)) + 30 
    }
    if(type %in% c("MX", "AL", "AF","A5", "A3") ){
      xmin1<- (gene_rescaled_exons[exonnum1,2] %>% pull(start)) - 30 
      xmax1<- (gene_rescaled_exons[exonnum2,3] %>% pull(end)) + 30 
    }
    str<- ggplot(gene_rescaled_exons, aes(xstart = start, xend = end, y = transcript_id)) +
      geom_rangeX(fill="grey90",color="black" ,height=0.6) +
      geom_intron(data = gene_rescaled_introns, aes(strand=strand), 
                  arrow=grid::arrow(ends = "last", length = grid::unit(0.005, "inches")),
                  color="black",
                  arrow.min.intron.length=2) +
      theme_void(base_size = base_size) +
      theme(legend.position = "none") +
      scale_fill_colorblind() + scale_color_brewer() +
      #ggtitle(ifelse(rhy_event_analysis %>% filter(eventID == eventInfo) %>% pull(novel)==TRUE ,paste(type, ":", name,"(", strand, ")", "*") ,paste(type, ":", name,"(", strand, ")"))) + #if it is novel, then add an asterik
      theme(plot.title = element_text(size=10))
    str<- str + geom_rect(aes(xmin = xmin1, xmax = xmax1, ymin = ymin, ymax = ymax), color = "black", fill= "transparent", size=0.5)
  }
  
  if(boxes==2){ # if you want 2 boxss
    xmin1<- (gene_rescaled_exons[exonnum1,2] %>% pull(start)) - 30  
    xmax1<- (gene_rescaled_exons[exonnum2,3] %>% pull(end)) + 30 
    xmin2<- (gene_rescaled_exons[exonnum3,2] %>% pull(start)) - 30 
    xmax2<- (gene_rescaled_exons[exonnum4,3] %>% pull(end)) + 30 
    
    str<- ggplot(gene_rescaled_exons, aes(xstart = start, xend = end, y = transcript_id)) +
      geom_rangeX(fill="black",color="black" ,height=0.6) +
      geom_intron(data = gene_rescaled_introns, aes(strand=strand), 
                  arrow=grid::arrow(ends = "last", length = grid::unit(0.005, "inches")),
                  color="black",
                  arrow.min.intron.length=2) +
      theme_void(base_size = base_size) +
      theme(legend.position = "none") +
      scale_fill_colorblind() + scale_color_brewer() +
      #ggtitle(ifelse(rhy_event_analysis %>% filter(eventID == eventInfo) %>% pull(novel)==TRUE ,paste(type, ":", name,"(", strand, ")", "*") ,paste(type, ":", name,"(", strand, ")"))) + #if it is novel, then add an asterik
      theme(plot.title = element_text(size=10))
    str<- str + geom_rect(aes(xmin = xmin1, xmax = xmax1, ymin = ymin, ymax = ymax), color = "black", fill= "transparent", size=0.2)
    str<- str + geom_rect(aes(xmin = xmin2, xmax = xmax2, ymin = ymin, ymax = ymax), color = "black", fill= "transparent", size=0.2)
  }
  return(str) 
}
###gene structures finalised (facets?)
##the number of exons is the count, but be sure which one you want to put as the list of transcript will be in alpabetical orer, better to have a look and put manually
theme_forbottom<-theme(axis.title.x = element_text(size=7),plot.title = element_text(size=7, hjust=-1, face="bold"), plot.subtitle = element_text(size=6, hjust=-1), axis.text.x= element_text(size=5), axis.text.y = element_text(size=8)) 
theme_fortops<-theme(axis.title.x = element_blank(), plot.title = element_text(size=6, hjust=-1,face="bold"), plot.subtitle = element_text(size=8, hjust=-1), axis.text.x= element_blank(),  axis.text.y = element_text(size=8)) 

###gathering of lower panel
 
SE_Arntl<- plot_grid(psi_graph("ENSMUSG00000055116;SE:7:112880009-112882464:112882502-112882728:+") +ggtitle("SE:Arntl(+)"), gene_structure("ENSMUSG00000055116","SE:7:112880009-112882464:112882502-112882728:+", 1, 6, 6), ncol=1, rel_heights = c(3,1), align = 'v', axis = 'l')
SE_Ythdc1<-plot_grid(psi_graph("ENSMUSG00000035851;SE:5:86968034-86968492:86968545-86969111:+") +ggtitle("SE:Ythdc1(+)"), gene_structure(	"ENSMUSG00000035851" ,"SE:5:86968034-86968492:86968545-86969111:+",1, 6, 6), ncol=1, rel_heights = c(3,1), align = 'v', axis = 'l')
MX_Lrrc28<- plot_grid(psi_graph("ENSMUSG00000030556;MX:7:67267878-67268822:67268859-67290852:67267878-67278043:67278083-67290852:-") +ggtitle("MX:Lrrc28(-)"), gene_structure("ENSMUSG00000030556","MX:7:67267878-67268822:67268859-67290852:67267878-67278043:67278083-67290852:-",1, 7, 16), ncol=1, rel_heights = c(3,1), align = 'v', axis = 'l')
MX_Nudt13<- plot_grid(psi_graph("ENSMUSG00000021809;MX:14:20344860-20345156:20345334-20354021:20344860-20350655:20350747-20354021:+" ) +ggtitle("MX:Nudt13(+)"), gene_structure("ENSMUSG00000021809","MX:14:20344860-20345156:20345334-20354021:20344860-20350655:20350747-20354021:+",1, 11, 2), ncol=1, rel_heights = c(3,1), align = 'v', axis = 'l')
RI_Cmtr1<- plot_grid(psi_graph("ENSMUSG00000024019;RI:17:29917612:29917703-29918194:29918262:+") +ggtitle("RI:Cmtr1(+)"), gene_structure("ENSMUSG00000024019","RI:17:29917612:29917703-29918194:29918262:+",1, 34, 34), ncol=1, rel_heights = c(3,1), align = 'v', axis = 'l')
RI_Ubxn1<- plot_grid(psi_graph("ENSMUSG00000071655;RI:19:8849372:8849478-8849582:8849651:+") +ggtitle("RI:Ubxn1(+)"), gene_structure("ENSMUSG00000071655","RI:19:8849372:8849478-8849582:8849651:+",1, 10, 10), ncol=1, rel_heights = c(3,1), align = 'v', axis = 'l')
AS5_Hmox2<-plot_grid(psi_graph("ENSMUSG00000004070;A5:16:4575034-4580423:4574991-4580423:+") +ggtitle("AS5:Hmox2(+)"), gene_structure("ENSMUSG00000004070","A5:16:4575034-4580423:4574991-4580423:+",1, 7, 1), ncol=1, rel_heights = c(3,1), align = 'v', axis = 'l')
AS5_Nmrk1<-plot_grid(psi_graph("ENSMUSG00000037847;A5:19:18609445-18613443:18609363-18613443:+") +ggtitle("AS5:Nmrk1(+)"), gene_structure("ENSMUSG00000037847","A5:19:18609445-18613443:18609363-18613443:+",1, 1, 10), ncol=1, rel_heights = c(3,1), align = 'v', axis = 'l')
AS3_Fmr1<-plot_grid(psi_graph("ENSMUSG00000000838;A3:X:67754385-67755837:67754385-67755912:+") +ggtitle("AS3Fmr1(+)"), gene_structure("ENSMUSG00000000838","A3:X:67754385-67755837:67754385-67755912:+",1, 15, 15), ncol=1, rel_heights = c(3,1), align = 'v', axis = 'l')
AS3_Wdr13<-plot_grid(psi_graph("ENSMUSG00000031166;A3:X:7997864-7998526:7997830-7998526:-") +ggtitle("AS3:Wdr13(-)"), gene_structure("ENSMUSG00000031166","A3:X:7997864-7998526:7997830-7998526:-",1, 8, 8), ncol=1, rel_heights = c(3,1), align = 'v', axis = 'l')
AF_Atp5g1<-plot_grid(psi_graph("ENSMUSG00000006057;AF:11:95965859-95966383:95966445:95965859-95966453:95966496:-") +ggtitle("AF:Atp5g1(-)"), gene_structure("ENSMUSG00000006057","AF:11:95965859-95966383:95966445:95965859-95966453:95966496:-",1, 10, 5), ncol=1, rel_heights = c(3,1), align = 'v', axis = 'l')
AF_Tef<-plot_grid(psi_graph("ENSMUSG00000022389;AF:15:81686622:81687087-81699145:81695615:81695921-81699145:+") +ggtitle("AF:Tef(+)"), gene_structure("ENSMUSG00000022389","AF:15:81686622:81687087-81699145:81695615:81695921-81699145:+",1, 5, 1), ncol=1, rel_heights = c(3,1), align = 'v', axis = 'l')
AL_Paxx<-plot_grid(psi_graph("ENSMUSG00000047617;AL:2:25345153:25345660-25349787:25349494:25349713-25349787:-") +ggtitle("AL:Paxx(-)"), gene_structure("ENSMUSG00000047617","AL:2:25345153:25345660-25349787:25349494:25349713-25349787:-",1, 8, 1), ncol=1, rel_heights = c(3,1), align = 'v', axis = 'l')
AL_Ddo<-plot_grid(psi_graph("ENSMUSG00000063428;AL:10:40513517-40515664:40517564:40513517-40550582:40550911:+") +ggtitle("AL:Ddo(+)"), gene_structure("ENSMUSG00000063428","AL:10:40513517-40515664:40517564:40513517-40550582:40550911:+",1, 5, 10), ncol=1, rel_heights = c(3,1), align = 'v', axis = 'l')





ggsave(paste(fig_dir,"/cyclingpsi_SE.svg",sep=""), plot=plot_grid(SE_Arntl,SE_Ythdc1, nrow=2), width=1.5, height = 3)
ggsave(paste(fig_dir,"/cyclingpsi_RI.svg",sep=""), plot=plot_grid(RI_Cmtr1,RI_Ubxn1, nrow=2), width=1.5, height = 3)
ggsave(paste(fig_dir,"/cyclingpsi_MX.svg",sep=""), plot=plot_grid(MX_Lrrc28,MX_Nudt13, nrow=2), width=1.5, height = 3)
ggsave(paste(fig_dir,"/cyclingpsi_A5.svg",sep=""), plot=plot_grid(AS5_Hmox2,AS5_Nmrk1, nrow=2), width=1.5, height = 3)
ggsave(paste(fig_dir,"/cyclingpsi_A3.svg",sep=""), plot=plot_grid(AS3_Fmr1,AS3_Wdr13, nrow=2), width=1.5, height = 3)
ggsave(paste(fig_dir,"/cyclingpsi_AL.svg",sep=""), plot=plot_grid(AL_Paxx,AL_Ddo, nrow=2), width=1.5, height = 3)
ggsave(paste(fig_dir,"/cyclingpsi_AF.svg",sep=""), plot=plot_grid(AF_Atp5g1,AF_Tef, nrow=2), width=1.5, height = 3)

rhythmic_tx <- read.csv(paste(list_dir,"/circadian_transcripts.csv", sep=""), header=TRUE) %>% dplyr::select("gene_id")
rhythmic_tx<-unique(rhythmic_tx$gene_id)
rhythmic_AS <- unique(cycling_AS$geneID)

library("ggvenn")
venn_rR <- ggvenn(list(rhythmic_tx=rhythmic_tx, rhythmic_AS=rhythmic_AS), columns = c('rhythmic_tx', 'rhythmic_AS'),
                  show_percentage = FALSE, fill_color = c("#434299", "#822d30"), fill_alpha = 0.3, auto_scale = FALSE, text_size = 4, set_name_size=4) 

#############
### Lies ####
#############

##############################
# barplot nb of cycling events


Fig3D <- cycling_AS %>%
  mutate(eventID=paste(geneID, eventID, sep=";")) %>%
  left_join(., psi_dataX[,c(1,14)]) %>%
  mutate(amplitude = ifelse(Amplitude > 0.1, "strong", "weak")) %>%
  filter(!is.na(Amplitude)) %>%
  ggplot(aes(y=factor(type, levels = event_order), group=factor(amplitude, levels=c("weak", "strong")), fill=factor(amplitude, levels=c("weak", "strong")))) +
  geom_bar(stat='count', width = 0.5, col="black") + 
  theme_classic() +
  scale_x_continuous(position = "top", expand = c(0,0), limits = c(0,150)) +
  theme(axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), text = element_text(size=base_size)) +
  xlab("number of rhythmic AS events") +
  scale_fill_manual(values=c("grey90", 'black')) + #scale_fill_manual(values=c("grey90", 'black'))
  labs(fill = "Amplitude") #+ theme(axis.text.y = element_blank())

cycling_AS %>%
  mutate(eventID=paste(geneID, eventID, sep=";")) %>%
  left_join(., psi_dataX[,c(1,14)]) %>%
  mutate(amplitude = ifelse(Amplitude > 0.1, "strong", "weak")) %>%
  filter(!is.na(Amplitude)) %>% 
  as.data.frame() %>% 
  dplyr::select(-c(alternative_transcripts, total_transcripts, primary_transcripts)) %>% 
  write.table(., "../lists/rhythmic_AS_events.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)



###########################
# calculate splicing entropy 
#for each gene, calculate the entropy defined as : Entropy = -Σi Ψi * log2(Ψi) where i is every isoform of this gene (or every event?)

entropy_df <- psi_data[,-5:-6] %>%
  group_by(geneID, eventID, time) %>%
  mutate(psi_av = mean(psi)) %>%
  ungroup() %>%
  dplyr::select(-psi) %>% 
  filter(psi_av > 0) %>%
  distinct() %>%
  mutate(psi_x_log2psi = psi_av*log2(psi_av)) %>%
  group_by(geneID, time) %>%
  mutate(entropy=-sum(psi_x_log2psi)) %>%
  ungroup()

################################
# get splicng entropy per GO term

# first get the go terms for genes in our list
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
GO_df <- getBM(attributes = c('ensembl_gene_id', 'go_id', 'name_1006', 'definition_1006', 'go_linkage_type', 'namespace_1003'), filters = 'ensembl_gene_id', values = unique(psi_data$geneID), mart = mart) %>%
  setNames(c('ensembl_gene_id', 'GO_term_accession', 'GO_term_name', 'GO_term_definition', 'GO_term_evidence_code', 'GO_domain'))

# keep all cellular component go terms that gather at list 50 genes 
nb_go_df <- data.frame()
for (i in unique(GO_df$GO_term_name)){
  nb_go_df <- rbind(nb_go_df, data.frame(GO = i, nb_genes = length(unique(filter(GO_df, GO_term_name == i)$ensembl_gene_id)))) #get the nb of genes involved in the given GO term
}

# create the df
list_GO_min <- filter(nb_go_df, nb_genes > 50)$GO # list of GO that contains min 50 genes
GO_df_filtered <- data.frame()
for (i in list_GO_min){
  df <- filter(psi_data[,-5:-6], geneID %in% filter(GO_df, GO_term_name == i)$ensembl_gene_id) %>% mutate(GO = i, GO_domain = unique(filter(GO_df, GO_term_name == i)$GO_domain))
  GO_df_filtered <- rbind(GO_df_filtered ,df)
}
GO_df_filtered <- rbind(GO_df_filtered %>% left_join(., nb_go_df, by = 'GO'), mutate(psi_data[,-5:-6], GO = "All", nb_genes = length(unique(psi_data$geneID)), GO_domain = "All")) 

G0_df_entropy <- GO_df_filtered %>%
  group_by(geneID, eventID, time) %>%
  mutate(psi_av = mean(psi)) %>%
  ungroup() %>%
  dplyr::select(-psi) %>% 
  filter(psi_av > 0) %>%
  distinct() %>%
  mutate(psi_x_log2psi = psi_av*log2(psi_av)) %>%
  group_by(GO_domain, GO, geneID,  time) %>%
  mutate(entropy=-sum(psi_x_log2psi)) %>%
  ungroup() %>%
  group_by(GO_domain, GO) %>%
  mutate(median_GO = median(entropy)) %>%
  ungroup()

G0_df_entropy_table <- dplyr::select(G0_df_entropy, c(GO_domain, GO, median_GO)) %>% distinct() %>% setNames(c("GO_domain", "GO", "median_entropy"))

supp_entropy_plot <- ggplot(filter(G0_df_entropy, GO %in% c("All", "circadian rhythm", "extracellular exosome")), aes(x=GO, y= entropy)) + 
  geom_boxplot(outlier.shape=NA) + 
  ylim(0,250) + 
  theme_classic()

GO_entropy_distribution <- ggplot(G0_df_entropy_table, aes(x=median_entropy)) +
  geom_density() + 
  theme_classic() +
  annotate(geom='line', x= filter(G0_df_entropy_table, GO == "All")$median_entropy, y=c(-Inf, Inf), linetype="dashed") +
  xlab("entropy") +
  ylab("fraction") +
  ggtitle(label = "Distribution of entropy per GO family")

# calculate rhythmicity of median entropy of each GO

G0_entropy_circadian_df <- G0_df_entropy %>%
  group_by(time, GO) %>%
  mutate(median_GO_time = median(entropy)) %>%
  ungroup() %>%
  mutate(inphase = cos(2*pi*time/24),
         outphase = sin(2*pi*time/24)) %>% 
  dplyr::select(c(GO_domain, GO, time, median_GO_time, inphase, outphase)) %>%
  nest_by(GO) %>%
  mutate(glm_pval = anova(glm(median_GO_time~inphase + outphase, data = data), glm(median_GO_time~1, data = data) , test = "LRT")[2,5]) %>%
  collect() %>%
  unnest(c(data)) %>%
  ungroup %>%
  distinct() %>%
  group_by(GO) %>%
  mutate(padj = p.adjust(glm_pval),
         amplitude=max(median_GO_time)-min(median_GO_time),
         r_amp = amplitude/mean(median_GO_time)) %>%
  ungroup()

###

interesting <- c("All", "regulation of RNA splicing","ribosome binding ","pheromone binding","odorant binding","mitochondrion organization","insulin receptor signaling pathway", "cellular response to insulin stimulus", "cell cycle","cholesterol metabolic process", "glucose homeostasis")

Fig3Ba <- filter(G0_df_entropy, GO %in% interesting) %>%
  mutate(GO = str_replace(GO, "cellular response to insulin stimulus", "cellular response\nto insulin stimulus"),
         GO = str_replace(GO, "insulin receptor signaling pathway",  "insulin receptor\nsignaling pathway")) %>%
  ggplot() + 
  geom_boxplot(aes(x=as.factor(time), y= entropy), outlier.shape=NA, fill = "grey") + 
  facet_wrap('GO', nrow=1, scales = "free_y", strip.position="bottom") +
  theme_classic() +
  theme(axis.title.x = element_blank(), strip.background = element_rect(fill = NA, colour = NA), strip.placement = "outside", strip.text = element_text(size=7)) 

#####
entropy_delta_df <- filter(G0_df_entropy, GO %in% interesting) %>%
  mutate(GO = str_replace(GO, "cellular response to insulin stimulus", "cellular response\nto insulin stimulus"),
         GO = str_replace(GO, "insulin receptor signaling pathway",  "insulin receptor\nsignaling pathway"),
         GO = str_replace(GO, "glucose homeostasis",  "glucose homeostasis"),
         GO = str_replace(GO, "cholesterol metabolic process",  "cholesterol\nmetabolic process")) %>%
  dplyr::select(GO, time, entropy, median_GO) %>%
  group_by(GO, time) %>%
  mutate(median_GO = median(entropy)) %>%
  ungroup() %>%
  dplyr::select(-entropy) %>%
  group_by( GO) %>%
  mutate(min_entropy = min(median_GO)) %>%
  ungroup() %>%
  mutate(relative_entropy = median_GO/max(median_GO)) %>%
  dplyr::select(GO, time, relative_entropy) %>%
  distinct() %>%
  #mutate(name2plot = ifelse(GO == "insulin receptor\nsignaling pathway" & time == 20, GO, 
  #                          ifelse(GO == "cellular response\nto insulin stimulus" & time == 20, GO,
  #                                ifelse(GO == "All" & time == 20, GO,
  #                                       ifelse(GO == "cell cycle" & time == 20, GO, NA)))))
  mutate(name2plot=ifelse(time == 20, GO,NA) )

write.csv(entropy_delta_df, paste(list_dir, "/entropy_delta_df.csv", sep=""), row.names = FALSE)
  
Fig3Bb <- ggplot(entropy_delta_df, aes(x=time, y=relative_entropy, group=GO, col=GO, label= name2plot)) +
  geom_point()  +
  geom_smooth(method = "loess",se=F) +
  ylab("Relative entropy") +
  theme_classic() +
  scale_color_manual(values=c( "grey70", "black", "#434299", "#822d30", "#730864", "#014F05", "purple")) +
  geom_text(nudge_x = 1.5, nudge_y = ifelse(entropy_delta_df$GO == "All", -0.03, ifelse(entropy_delta_df$GO == "cholesterol metabolic process", -1.2, 0.03)), size = 2.5) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0,24), breaks = c(0, 4, 8, 12, 16, 20)) +
  xlab("Time (CT)")

#Fig3B <- plot_grid(Fig3Ba, Fig3Bb, nrow=1, rel_widths = c(2,1))
#Fig3B
#############
####test#####
#############
entropy_delta_dfx <- G0_df_entropy %>%
  mutate(GO = str_replace(GO, "cellular response to insulin stimulus", "cellular response\nto insulin stimulus"),
         GO = str_replace(GO, "insulin receptor signaling pathway",  "insulin receptor\nsignaling pathway")) %>%
  dplyr::select(GO, time, entropy, median_GO) %>%
  group_by(GO, time) %>%
  mutate(median_GO = median(entropy)) %>%
  ungroup() %>%
  dplyr::select(-entropy) %>%
  group_by(GO) %>%
  mutate(min_entropy = min(median_GO)) %>%
  ungroup() %>%
  mutate(relative_entropy = median_GO/max(median_GO)) %>%
  dplyr::select(GO, time, relative_entropy) %>%
  distinct() 
 
Fig3Bbx <- ggplot(entropy_delta_dfx, aes(x=time, y=relative_entropy, group=GO, col=GO, label=GO)) +
  geom_point()  +
  geom_smooth(method = "loess",se=F) +
  ylab("Relative entropy") +
  theme_classic() +
  #scale_color_manual(values=c( "grey70", "black", "#434299", "#822d30")) +
  geom_text(nudge_x = .8, nudge_y = ifelse(entropy_delta_df$GO == "All", 0, 0), size = 2.5) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0,24), breaks = c(0, 4, 8, 12, 16, 20)) +
  xlab("Time (CT)") + facet_wrap('GO')

ggsave(paste(fig_dir,"/test_all_entropy_fixedY.svg",sep=""), plot=Fig3Bbx, width=100, height = 50, limitsize=F)

###
label_data <- dplyr::select(GO_df_filtered, c("GO", "nb_genes")) %>%
  distinct() %>%
  left_join(G0_df_entropy_table, ., by='GO', keep=FALSE) %>%
  filter(median_entropy > 90 | median_entropy < 10) %>%
  filter(!GO %in% c("molecular_function", "cellular_component", "biological_process")) %>%
  mutate(GO = str_replace(GO, "regulation", "regulation\n"),
         GO = str_replace(GO, "response", "response\n"),
         GO = str_replace(GO, "receptor", "receptor\n"),
         GO = str_replace(GO, "phosphorylation", "\nphosphorylation"),
         GO = str_replace(GO, "assembly", "\nassembly"),
         GO = str_replace(GO, "recombination", "\nrecombination"),
         GO = str_replace(GO, "exosome", "\nexosome"),
         GO = str_replace(GO, "processing", "\nprocessing"),
         GO = str_replace(GO, "kinase B", "kinase B\n"),
         GO = str_replace(GO, ",", ",\n"),
         id=as.numeric(row.names(.))) %>%
  mutate(nb_genes = ifelse(GO %in% c("hydrolase activity,  acting on glycosyl bonds", "positive regulation  of peptidyl-tyrosine  phosphorylation"), nb_genes+30, nb_genes+10))


number_of_bar <- nrow(label_data)
angle <-  90 - 360 * (label_data$id-0.5) /number_of_bar
label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$angle<-ifelse(angle < -90, angle+180, angle)
label_data$ypos <- log2(label_data$median_entropy) +1
label_data$vjust <- (rep(0.5, nrow(label_data)))
label_data$vjust[5] <- label_data$vjust[5] - 2 
label_data$ypos[5] <- label_data$ypos[5] + 10
label_data$ypos[6] <- label_data$ypos[6] + 8
label_data$vjust[6] <- label_data$vjust[6] + 0.5
label_data$ypos[10] <- label_data$ypos[10] + 1
label_data$vjust[10] <- label_data$vjust[10] - 0.5
label_data$ypos[4] <- label_data$ypos[4] - 1
label_data$vjust[4] <- label_data$vjust[4] - 1
label_data$ypos[12] <- label_data$ypos[12] + 5
label_data$vjust[12] <- label_data$vjust[12] - 1
label_data$ypos[8] <- label_data$ypos[8] + 3
label_data$vjust[8] <- label_data$vjust[8] + 3


Fig3A <- dplyr::select(GO_df_filtered, c("GO", "nb_genes")) %>%
  distinct() %>%
  left_join(G0_df_entropy_table, ., by='GO', keep=FALSE) %>%
  filter(median_entropy > 90 | median_entropy < 10) %>%
  filter(!GO %in% c("molecular_function", "cellular_component", "biological_process")) %>%
  mutate(GO = str_replace(GO, "regulation", "regulation\n"),
         GO = str_replace(GO, "response", "response\n"),
         GO = str_replace(GO, "receptor", "receptor\n"),
         GO = str_replace(GO, "phosphorylation", "\nphosphorylation"),
         GO = str_replace(GO, "assembly", "\nassembly"),
         GO = str_replace(GO, "recombination", "\nrecombination"),
         GO = str_replace(GO, "exosome", "\nexosome"),
         GO = str_replace(GO, "processing", "\nprocessing"),
         GO = str_replace(GO, "kinase B", "kinase B\n"),
         GO = str_replace(GO, ",", ",\n")) %>%
  ggplot(aes(x=GO, fill= log2(median_entropy), y=log2(median_entropy))) +
  geom_bar(stat='identity', width=0.5, color="black") +
  coord_polar() +
  ylim(-5,20) +
  theme_classic() +
  theme(axis.line = element_blank(), axis.title=element_blank(), axis.text.y=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), plot.margin = unit(rep(-1,4), "cm"), legend.position="none") +
  geom_text(data=label_data, aes(x=GO, y=ypos, label=GO, hjust=hjust, vjust=vjust), color="black", size=2.5, inherit.aes = FALSE ) +
  scale_fill_gradient2(low="#822d30", high="#434299", midpoint = log2(17))

#####unrolled
interesting <- c("All","cell morphogenesis","RNA polymerase II transcription regulator complex","RNA splicing","cholesterol binding","peptidyl-serine phosphorylation","positive regulation of cell growth","cell cycle",
                  "pheromone binding", "odorant binding","mitochondrion organization","insulin receptor signaling pathway", "cellular response to insulin stimulus","cholesterol metabolic process", "glucose homeostasis")

Fig3A <- dplyr::select(GO_df_filtered, c("GO", "nb_genes")) %>%
  distinct() %>%
  left_join(G0_df_entropy_table, ., by='GO', keep=FALSE) %>%
  #filter(median_entropy > 90 | median_entropy < 10) %>%
  filter(GO  %in% interesting) %>%
  filter(!GO %in% c("molecular_function", "cellular_component", "biological_process")) %>%
  mutate(
         GO = str_replace(GO, "tumor", "tumor\n"),
         GO = str_replace(GO, "peptidyl", "peptidyl\n"),
         #GO = str_replace(GO, "ion", "ion\n"),
         #GO = str_replace(GO, "receptor", "receptor\n"),
         #GO = str_replace(GO, "phosphorylation", "\nphosphorylation"),
         #GO = str_replace(GO, "assembly", "\nassembly"),
         #GO = str_replace(GO, "recombination", "\nrecombination"),
         #GO = str_replace(GO, "exosome", "\nexosome"),
         #GO = str_replace(GO, "processing", "\nprocessing"),
         #GO = str_replace(GO, "kinase B", "kinase B\n"),
         #GO = str_replace(GO, "ion", "ion\n"),
         GO = str_replace(GO, ",", ",\n")) 
         
Fig3A_plot<- Fig3A %>% mutate(entropy_by_gene= median_entropy/nb_genes) %>%
  ggplot(aes(x=entropy_by_gene,  y=GO)) +
  geom_bar(stat='identity', color="black", fill="#165488", width=0.6) +
  #coord_polar() +
  #ylim(-5,20) +
  theme_classic() +
  xlab("Relative entropy(GO median entopy/nb_genes)")
  #theme(axis.line = element_blank(), axis.title=element_blank(), axis.text.y=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), plot.margin = unit(rep(-1,4), "cm"), legend.position="none") +
  #geom_text(data=label_data, aes(x=GO, y=ypos, label=GO, hjust=hjust, vjust=vjust), color="black", size=2.5, inherit.aes = FALSE ) +
  #scale_fill_gradient2(low="#822d30", high="#434299", midpoint = 0.5 ) #+ 
  #theme(axis.text. =element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste(fig_dir,"/entropy_unrolled.svg",sep=""), plot=Fig3A_plot, width=7, height = 5, limitsize= FALSE) 
dplyr::select(GO_df_filtered, c("GO", "nb_genes")) %>%
  distinct() %>%
  left_join(G0_df_entropy_table, ., by='GO', keep=FALSE) %>%
  write.table(., "../lists/GO_entropy_barplots.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

#### entropy distribution rhythmic expr
#######################################

circadian_expression <- read.csv(paste(list_dir,"/circadian_transcripts.csv", sep=""), header = TRUE)

Fig3C <- dplyr::select(entropy_df, -c(psi_av, psi_x_log2psi, eventID)) %>%
  distinct() %>%
  group_by(geneID) %>%
  mutate(entropy_per_gene = mean(entropy)) %>%
  dplyr::select(-c(time, entropy)) %>%
  distinct() %>%
  mutate(circadian = ifelse(geneID %in% circadian_expression$gene_id, 'rhythmic', 'NR')) %>%
  ggplot(aes(y=entropy_per_gene, x=circadian)) +
  geom_violin(alpha=0.5, aes(fill=circadian)) + 
  geom_boxplot(width=0.1, outlier.shape=NA) + 
  theme_classic() +
  xlab("entropy per gene") +
  ylab("fraction") +
  scale_fill_manual(values=c("grey", "#434299")) +
  ylim(0,8) +
  theme(legend.position="bottom", legend.title = element_blank()) + 
  stat_compare_means(method = "t.test")

dplyr::select(entropy_df, -c(psi_av, psi_x_log2psi, eventID)) %>%
  distinct() %>%
  group_by(geneID) %>%
  mutate(entropy_per_gene = mean(entropy)) %>%
  dplyr::select(-c(time, entropy)) %>%
  distinct() %>%
  setNames(c("ensembl_gene_id", "entropy")) %>%
  mutate(circadian = ifelse(ensembl_gene_id %in% circadian_expression$gene_id, 'rhythmic', 'NR')) %>%
  left_join(., getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "mgi_symbol"), values = .$ensembl_gene_id, mart = mart)) %>%
  write.table(., "../lists/3_entropy_R_NR_Fig3C.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
  
  
####

fdr_cutoff <- 0.05
count_threshold <- 5
base_size<-9
retrieve_from_biomart <- FALSE
tx_names_with_version <- TRUE

files <- list.files(path = salmon_path, pattern = ".sf", 
                    full.names = TRUE, include.dirs = FALSE, no.. = TRUE, recursive = TRUE)
names(files) <- str_extract(files, "([CK][TO]\\d+_\\d)(\\.sf)$", group = 1) 
txi.salmon <- tximport(files, type = "salmon", txOut = TRUE, countsFromAbundance = "no", ignoreTxVersion = TRUE, 
                       importer = function(x) read.delim(x))
exp_design <- data.frame(batch = factor(str_extract(colnames(txi.salmon$abundance), "(?<=_)\\d")),
                         time = as.numeric(str_extract(colnames(txi.salmon$abundance), "\\d+")),
                         cond = factor(str_extract(colnames(txi.salmon$abundance), "^[A-Z]+")),
                         sample_id= colnames(txi.salmon$abundance))
exp_design[exp_design$cond == "KO",]$batch = 2

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

keep <- rowSums(counts > count_threshold) >= 2   ##this is a universal threshold cosidering 16 samples

dds <- dds[keep, ]

abundance_no_batch_effect <- txi.salmon$abundance %>% data.frame %>%
  magrittr::set_rownames(str_extract(rownames(.), "([A-Za-z0-9-]+)")) %>% 
  {log2(. + 0.001)} %>%
  {limma::removeBatchEffect(.[rownames(dds), ], batch = colData(dds)$batch)}
  

Wang <- readxl::read_excel("../data/Wang_2017.xlsx", sheet="Complete dataset") %>% 
  dplyr::select(c(Gene.names, phase.WT, amp.WT, relamp.WT, qv.WT)) %>% 
  setNames(c("mgi_symbol", "phase", "amp", "relamp", "q_value")) %>% 
  filter(q_value < 0.1) %>%
  mutate(dataset = "proteome") 


Robles_phospho <- readxl::read_excel("../data/circadian_phosphoproteome_Robles2017_supp.xlsx", sheet="B - Cycling q<0.1", na = "NaN") %>%
  mutate(CT0 = rowMeans(dplyr::select(., contains("Intensity_CT0")), na.rm = TRUE),
         CT3 = rowMeans(dplyr::select(., contains("Intensity_CT3")), na.rm = TRUE),
         CT6 = rowMeans(dplyr::select(., contains("Intensity_CT6")), na.rm = TRUE),
         CT9 = rowMeans(dplyr::select(., contains("Intensity_CT9")), na.rm = TRUE),
         CT12 = rowMeans(dplyr::select(., contains("Intensity_CT12")), na.rm = TRUE),
         CT15 = rowMeans(dplyr::select(., contains("Intensity_CT15")), na.rm = TRUE),
         CT18 = rowMeans(dplyr::select(., contains("Intensity_CT18")), na.rm = TRUE),
         CT21 = rowMeans(dplyr::select(., contains("Intensity_CT21")), na.rm = TRUE),
         CT24 = rowMeans(dplyr::select(., contains("Intensity_CT24")), na.rm = TRUE),
         CT27 = rowMeans(dplyr::select(., contains("Intensity_CT27")), na.rm = TRUE),
         CT30 = rowMeans(dplyr::select(., contains("Intensity_CT30")), na.rm = TRUE),
         CT33 = rowMeans(dplyr::select(., contains("Intensity_CT33")), na.rm = TRUE),
         CT36 = rowMeans(dplyr::select(., contains("Intensity_CT36")), na.rm = TRUE),
         CT39 = rowMeans(dplyr::select(., contains("Intensity_CT39")), na.rm = TRUE),
         CT42 = rowMeans(dplyr::select(., contains("Intensity_CT42")), na.rm = TRUE),
         CT45 = rowMeans(dplyr::select(., contains("Intensity_CT45")), na.rm = TRUE)) %>%
  dplyr::select(!contains("Intensity")) %>%
  dplyr::select(c(Gene_names, Phase, q_value, contains("CT"))) %>%
  mutate(mean_D1 = rowMeans(.[, 4:11], na.rm = TRUE),
         mean_D2 = rowMeans(.[, 12:19], na.rm = TRUE)) %>%
  rowwise() %>%
  mutate(amp = mean(max(c_across(CT0:CT21), na.rm = TRUE)-min(c_across(CT0:CT21), na.rm = TRUE), max(c_across(CT24:CT45), na.rm = TRUE)-min(c_across(CT24:CT45), na.rm = TRUE)),
         ramp = mean(c(abs((max(c_across(CT0:CT21), na.rm = TRUE)-mean_D1)/mean_D1), abs((max(c_across(CT24:CT45), na.rm = TRUE)-mean_D2)/mean_D2)))) %>% # ramp seems not to be the best metric (since there are negative means)
  dplyr::select(c(Gene_names, Phase, amp, ramp, q_value)) %>%
  setNames(c("mgi_symbol", "phase", "amp", "relamp", "q_value")) %>%
  mutate(dataset = "phospho_proteome") 


list_SF_expr <- getBM(filters = "ensembl_gene_id", attributes = "mgi_symbol", values = unique(filter(GO_df, GO_term_name == "RNA splicing")$ensembl_gene_id), mart = mart)$mgi_symbol

SF_plot <- rbind(Wang, Robles_phospho) %>%
  mutate(mgi_symbol = str_replace(mgi_symbol, ".*;", "")) %>%
  rowwise() %>%
  mutate(relamp_dataset = ifelse(dataset == "proteome", amp/median(Wang$amp), amp/median(Robles_phospho$amp))) %>% 
  filter(mgi_symbol %in% list_SF_expr) %>% 
  group_by(phase, amp, mgi_symbol) %>%
  filter(row_number()==1) %>%
  ungroup() %>%
  ggplot(aes(x=phase, y=mgi_symbol, size=relamp_dataset, fill=dataset)) +
  geom_point(alpha=0.4, shape=21, color="black") + 
  theme_classic() + 
  scale_size(range = c(0.1, 15), breaks = c(1,3), name = "Amplitude") +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(), legend.position = "right", 
        panel.grid.major.y = element_line(color = "grey90", size = 0.25), 
        legend.box="vertical", legend.key.size = unit(0.1, 'cm'), 
        legend.key.height = unit(10, 'points')) +
  scale_fill_manual(values = c("#822d30", "#434299"), name = "") 

list_SF_plotted <- rbind(Wang, Robles_phospho) %>%
  mutate(mgi_symbol = str_replace(mgi_symbol, ".*;", "")) %>%
  rowwise() %>%
  mutate(relamp_dataset = ifelse(dataset == "proteome", amp/max(Wang$amp), amp/max(Robles_phospho$amp))) %>% 
  filter(mgi_symbol %in% list_SF_expr) %>% 
  group_by(phase, amp, mgi_symbol) %>%
  filter(row_number()==1) %>%
  ungroup() %>%
  .$mgi_symbol

Nanopore_boxplot <- abundance_no_batch_effect %>%
  as.data.frame() %>%
  mutate(ensembl_transcript_id = rownames(.)) %>%
  left_join(., getBM(filters = "ensembl_transcript_id", attributes = c("mgi_symbol", "ensembl_transcript_id"), values = .$ensembl_transcript_id, mart = mart)) %>%
  dplyr::select(!contains("KO")) %>%
  group_by(mgi_symbol) %>%
  mutate(across(CT0_1:CT8_2, ~ sum(.))) %>%
  ungroup() %>%
  dplyr::select(-ensembl_transcript_id) %>%
  distinct() %>%
  pivot_longer(!mgi_symbol, values_to = "expression", names_to = "timepoint") %>%
  filter(mgi_symbol %in% list_SF_plotted) %>% 
  ggplot(aes(y=mgi_symbol, x=expression)) +
  geom_boxplot(outlier.shape = NA) +
  #scale_y_discrete(position = "right") +
  theme_classic() +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size=6), panel.grid.major.y = element_line(color = "grey90", size = 0.25)) + 
  xlim(0, 25)

Fig3G <- plot_grid(Nanopore_boxplot, SF_plot, align = 'h', axis = "b", rel_widths = c(1,2))

abundance_no_batch_effect %>%
  as.data.frame() %>%
  mutate(ensembl_transcript_id = rownames(.)) %>%
  left_join(., getBM(filters = "ensembl_transcript_id", attributes = c("mgi_symbol", "ensembl_transcript_id"), values = .$ensembl_transcript_id, mart = mart)) %>%
  dplyr::select(!contains("KO")) %>%
  group_by(mgi_symbol) %>%
  mutate(across(CT0_1:CT8_2, ~ sum(.))) %>%
  ungroup() %>%
  dplyr::select(-ensembl_transcript_id) %>%
  distinct() %>%
  pivot_longer(!mgi_symbol, values_to = "expression", names_to = "timepoint") %>%
  filter(mgi_symbol %in% list_SF_plotted) %>%
  pivot_wider(names_from = timepoint, values_from = expression) %>% 
  write.table(., "../lists/phospho_proteome.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

ggsave(paste(fig_dir,"/3_rythmic_cycling_events.svg",sep=""), plot=Fig3D, width=3, height = 2)
ggsave(paste(fig_dir,"/3_cycling_AS.svg",sep=""), plot=Fig3E, width=4, height = 2.5)
ggsave(paste(fig_dir,"/3_GO_entropy.svg",sep=""), plot=Fig3A_plot, width=7, height = 3)
ggsave(paste(fig_dir,"/3_entropy_R_NR.svg",sep=""), plot=Fig3C, width=2.5, height = 4)
ggsave(paste(fig_dir,"/3_Fig3Bb.svg",sep=""), plot=Fig3Bb, width=2.5, height = 3.5)
ggsave(paste(fig_dir,"/3_venn2.svg",sep=""), plot=venn_rR, width=2.5, height = 3.5)

####################
### Supplemental ### 
####################

# overall entropy per gene accross time 
Overall_entropy <- ggplot(entropy_df, aes(x=time, y=entropy, group=time)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  ylim(0,3)
### Koike ChIP enrichment ###
#############################
library(ggpubr)

binding_factors <- c('Per1', 'Per2', 'Cry1', 'Cry2', 'Bmal1', 'Clock', 'Npas2')
Koike_df <- data.frame()
for (factor in binding_factors){
  table <- read.table(paste("../data/Koike_", factor, "_mm39.bed", sep="")) %>% 
    setNames(c('chr', 'start', 'end')) %>% 
    mutate(mid_position=(as.numeric(start)+as.numeric(end))/2, 
           binding_protein=factor)
  Koike_df <- rbind(Koike_df, table)
}

## get the "splicing choice site" (first bondary where machinery needs to make the decision)
Koike_psi_df <- cycling_AS[,c(1,2,5,7)] %>%
  mutate(eventID = paste(geneID, eventID, sep=";")) %>%
  left_join(psi_dataX[,c(1,14)], ., by='eventID') %>%
  separate(eventID, into = c("gene", "event"), sep=";") %>%
  mutate(strand = substr(event, nchar(event), nchar(event)),
         type = sub(":.*", "", event),
         circadian = ifelse(!is.na(mgi_symbol), "rhythmic", "NR"),
         event = substr(event, 1, nchar(event) - 2),
         chr = paste("chr", sub(".*:(.*?):.*", "\\1", event), sep="")) %>%
  filter(!type == "AF") %>% # I remove AF because it's not splicing
  mutate(position = case_when(
    strand == "+" & type %in% c("SE", "MX", "A3", "AL") ~ sub(".*:.*:(.*?)-(.*?):.*", "\\1", event), 
    strand == "+" & type %in% c("RI", "A5") ~ sub(".*?:.*?:.*?:(\\d+)-.*", "\\1", event),  
    strand == "-" & type %in% c("AL", "A3", "A5", "MX", "SE") ~ sub(".*-(.*)", "\\1", event), 
    strand == "-" & type == "RI" ~sub(".*-(\\d+):.*$", "\\1", event)),
    position=as.numeric(position)) %>%
  dplyr::select(-geneID)


BS_df <- data.frame()

for (i in 1:nrow(Koike_psi_df)){
  list_same_chr <- filter(Koike_df, chr == Koike_psi_df[i,8]) # filter positions on the same chr
  df <- data.frame(data.frame(matrix(NA, nrow = 1)))
  for (factor in binding_factors){
    x <- filter(list_same_chr, binding_protein == factor)
    closest_BS <- x$mid_position[which.min(abs(x$mid_position - Koike_psi_df[i,9]))] # get the closest binding site to the given TSS
    closest_BS <- ifelse(nrow(list_same_chr) == 0, "NA", closest_BS)  
    df <- cbind(df, data.frame(factor=closest_BS) %>% setNames(factor))
  }
  BS_df <- rbind(BS_df, df)
}

Koike_psi_df <- cbind(Koike_psi_df, BS_df[,-1]) %>%
  mutate(across(Per1:Npas2, ~ as.numeric(.)),
         across(Per1:Npas2, ~ abs(. - position)))

distance_Koike_event_plot <- Koike_psi_df %>%
  pivot_longer(Per1:Npas2, values_to = 'distance', names_to = 'binding_factor') %>%
  ggplot(aes(x=circadian, y=distance)) +
  geom_boxplot(outlier.shape = NA, aes(fill = circadian), alpha=0.5) +
  facet_wrap('binding_factor', nrow=1) +
  theme_classic() + 
  stat_compare_means(method = "t.test", size=2.5) +
  theme(axis.title.x = element_blank()) +
  scale_fill_manual(values=c("grey", "#434299"))


#### Overlap splicing/expression
################################

circadian_expression <- read.csv(paste(list_dir,"/circadian_transcripts.csv", sep=""), header = TRUE)

library("ggvenn")
venn_expr_splicing <- ggvenn(list(rhythmic_transcript=circadian_expression$gene_id, rhythmic_splicing=cycling_AS$geneID), 
                             columns = c('rhythmic_transcript', 'rhythmic_splicing'), show_percentage = FALSE, 
                             fill_color = c("#434299", "#822d30"), 
                             fill_alpha = 0.3, 
                             auto_scale = FALSE, 
                             text_size = 4,
                             set_name_size=4) 

###########################################
###### Correlation expression/entropy #####
###########################################

# replace negative values by 0 

abundance_no_batch_effect2 <- abundance_no_batch_effect
abundance_no_batch_effect2[abundance_no_batch_effect2 < 0] <- 0

# correlation entropy/expression per GO
entropy_expr_GO_df <- as.data.frame(abundance_no_batch_effect2) %>% 
  mutate(ensembl_transcript_id = row.names(.)) %>%
  rowwise() %>%
  mutate(mean_expr_WT = mean(CT0_1:CT8_2)) %>%
  dplyr::select(ensembl_transcript_id, mean_expr_WT) %>%
  left_join(., getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id'), filters = 'ensembl_transcript_id', values = .$ensembl_transcript_id, mart = mart)) %>%
  group_by(ensembl_gene_id) %>%
  mutate(mean_expr_WT = base::mean(mean_expr_WT, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(ensembl_gene_id, mean_expr_WT) %>%
  distinct() %>%
  left_join(dplyr::rename(GO_df_filtered, 'ensembl_gene_id'='geneID'), .) %>%
  group_by(GO) %>%
  mutate(median_expr_WT = median(mean_expr_WT, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(GO, median_expr_WT) %>%
  distinct() %>%
  left_join(., G0_df_entropy_table, by='GO')

correlation_entropy_expr_GO <- ggplot(entropy_expr_GO_df, aes(x=median_expr_WT, y=median_entropy)) +
  geom_point(size=0.7) +
  theme_classic() + 
  stat_cor(method = "pearson", label.y = 149) +
  ylab("median entropy per GO in WT") +
  xlab("median expression per GO in WT")

# correlation entropy/expression per gene
entropy_expr_gene_df <- as.data.frame(abundance_no_batch_effect2) %>% 
  mutate(ensembl_transcript_id = row.names(.)) %>%
  rowwise() %>%
  mutate(mean_expr_WT = mean(CT0_1:CT8_2)) %>%
  dplyr::select(ensembl_transcript_id, mean_expr_WT) %>%
  left_join(., getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id'), filters = 'ensembl_transcript_id', values = .$ensembl_transcript_id, mart = mart)) %>%
  group_by(ensembl_gene_id) %>%
  mutate(mean_expr_WT = base::mean(mean_expr_WT, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(ensembl_gene_id, mean_expr_WT) %>%
  distinct() %>%
  left_join(dplyr::rename(GO_df_filtered, 'ensembl_gene_id'='geneID'), .) %>%
  dplyr::select(ensembl_gene_id, mean_expr_WT) %>%
  distinct() %>%
  left_join(., dplyr::select(entropy_df, c(geneID, time, entropy)) %>% distinct() %>% group_by(geneID) %>% mutate(mean_entropy = mean(entropy, na.rm = TRUE)) %>% ungroup() %>% dplyr::select(geneID, mean_entropy) %>% distinct() %>% setNames(c('ensembl_gene_id', 'mean_entropy')))
  
correlation_entropy_expr_gene <- ggplot(entropy_expr_gene_df, aes(x=mean_expr_WT, y=mean_entropy)) +
  geom_point(size=0.7) +
  theme_classic() + 
  stat_cor(method = "pearson", label.y = 35) +
  ylab("mean entropy per gene in WT") +
  xlab("mean expression per gene in WT")

# also test with nb of isoform expressed

nb_isoforms_df <- read.csv("../results/expr_tx_list.csv") %>%
  group_by(gene_id) %>%
  mutate(nb_isoforms = n()) %>%
  ungroup() %>%
  dplyr::select(gene_id, nb_isoforms) %>%
  setNames(c('ensembl_gene_id', 'nb_isoforms')) %>%
  distinct()

correlation_entropy_nbiso_gene <- left_join(entropy_expr_gene_df, nb_isoforms_df) %>%
  ggplot(aes(x=nb_isoforms, y=mean_entropy)) +
  geom_point(size=0.7) +
  theme_classic() + 
  stat_cor(method = "pearson", label.y = 35) +
  ylab("mean entropy per gene in WT") +
  xlab("number of isoform per gene")

# per gene entropy rhythmicity analysis 

#gene_entropy_circadian_df <- dplyr::select(entropy_df, c(geneID, time, entropy)) %>% 
#  distinct() %>%
#  mutate(inphase = cos(2*pi*time/24),
#         outphase = sin(2*pi*time/24)) %>%
#  nest_by(geneID) %>%
#  mutate(glm_pval = anova(glm(entropy~inphase + outphase, data = data), glm(entropy~1, data = data) , test = "LRT")[2,5]) %>%
#  collect() %>%
#  unnest(c(data)) %>%
#  ungroup() %>%
#  distinct() %>%
#  group_by(geneID) %>%
#  mutate(padj = p.adjust(glm_pval),
#         amplitude=max(entropy)-min(entropy),
#         r_amp = amplitude/mean(entropy)) %>%
#  ungroup() %>%
#  dplyr::select(geneID, padj, amplitude, r_amp) %>%
#  distinct() %>%
#  rename('geneID' = 'ensembl_gene_id') %>%
#  left_join(., getBM(attributes = c('ensembl_gene_id', 'mgi_symbol'), filters = 'ensembl_gene_id', values = .$ensembl_gene_id, mart = mart))

ggsave(paste(fig_dir,"/3_correlation_entropy_expr_gene.svg",sep=""), plot=correlation_entropy_expr_gene, width=3, height = 3)
ggsave(paste(fig_dir,"/3_correlation_entropy_expr_GO.svg",sep=""), plot=correlation_entropy_expr_GO, width=3, height = 3)
ggsave(paste(fig_dir,"3_venn_exp_splicing.svg",sep=""), plot=venn_expr_splicing, width=4, height = 3)
ggsave(paste(fig_dir,"/3_overall_entropy.svg",sep=""), plot=Overall_entropy, width=4, height = 3)
ggsave(paste(fig_dir,"/3_phosphoproteome.svg",sep=""), plot=Fig3G, width=7, height = 7)
ggsave(paste(fig_dir,"/3_entropy_cycling_normlaised.svg",sep=""), plot=Fig3Bb, width=3.5, height = 3)

###gathering and saving

#F3<- plot_grid(
#  plot_grid(Fig3A, Fig3B, Fig3C, rel_widths = c(1,2,1), nrow=1),
#  NULL,
#  plot_grid(NULL, Fig3D, Fig3E, rel_widths = c(1,2,2), nrow=1),
#  NULL,
#  plot_grid(Fig3F, Fig3G, rel_widths = c(2,1.5), nrow=1),
#  rel_heights = c(1,0.1,1,0.1,2.5), ncol=1)

#ggsave("../results/figure3.png", F3, width=10, height=14)
#ggsave("../results/figure3.pdf", F3, width=11, height=14)
