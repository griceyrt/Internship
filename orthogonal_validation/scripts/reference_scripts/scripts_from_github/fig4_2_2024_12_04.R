
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
library(UpSetR)
library(grid)
library(gridExtra)
library(org.Mm.eg.db)
library(GenomicFeatures)
library(clusterProfiler)
library(AnnotationHub)
library(ggrepel)
library(ggh4x)
library(ReactomePA)
library(GenomicFeatures)
library(ggtranscript)

setwd("C:/Users/kmamgain/prettycoolPores_new/scripts")   #change here

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
mycolors_class <- c("Novel-novel"= "#fc8d62","Novel-Annotated"= "#8da0cb","Annotated"= "#66c2a5")
scale_fill <- scale_fill_manual(values = c("#BDD7E7", "#6BAED6", "#3182BD", "#08519C", "#EFF3FF"))
scale_color <- scale_color_manual(values = c("#BDD7E7", "#6BAED6", "#3182BD", "#08519C", "#EFF3FF"))
#scale_fill <- scale_fill_grey()
#scale_color <- scale_color_grey()
x_theme<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
event_order <- c("AF", "AL", "A3", "A5", "RI", "MX", "SE")

#####################
##### Functions ##### 
#####################

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

#GeomRange <- ggplot2::ggproto("GeomRange", ggplot2::GeomTile,
#                              required_aes = c("xstart", "xend", "y"),
 #                             default_aes = aes(
 #                               fill = "grey",
 ##                               colour = "black",
 #                               size = 0.25,
 #                               linetype = 1,
 #                               alpha = NA,
 #  ##                             height = NA
 #                             ),
  #                            setup_data = function(data, params) {
  #                              # modified from ggplot2::GeomTile
  #                              data$height <- data$height %||% params$height %||% 0.5
  #                              
  #                              transform(
  #                                data,
  #                                xmin = xstart,
  #                                xmax = xend,
  #                                ymin = y - height / 2,
   #                               ymax = y + height / 2,
   #                               height = NULL
   ##                             )
   #                           },
    #                          draw_panel = function(self,
    ##                                                data,
    #                                                panel_params,
     #                                               coord,   radius = grid::unit(6, "pt"),
     #                                               vjust = NULL,
     ##                                               lineend = "butt",
     #                                               linejoin = "mitre") {
     #                           if (!coord$is_linear()) {
      #                            # prefer to match geom_curve and warn
      ##                            # rather than copy the implementation from GeomRect for simplicity
      #                            # also don'think geom_range would be used for non-linear coords
      #                            warn("geom_ is not implemented for non-linear coordinates")
      #                          }
                                
                                
                                
   #                             coords <- coord$transform(data, panel_params)
                                # K added part here inspired from https://stackoverflow.com/questions/64355877/round-corners-in-ggplots-geom-tile-possible
   #                             lapply(1:length(coords$xmin), funion(i)) {      
   ##                               
    #                              grid::roundrectGrob(
    #                                coords$xmin[i], coords$ymax[i],
    ##                                width = (coords$xmax[i] - coords$xmin[i]),
    #                                height = (coords$ymax[i] - coords$ymin)[i],
     #                               r = radius,
     ##                               default.units = "native",
     #                               just = c("left", "top"),
      ##                              gp = grid::gpar(
      #                                col = coords$colour[i],
       ##                               fill = alpha(coords$fill[i], coords$alpha[i]),
       #                               lwd = coords$size[i] * .pt,
       #                               lty = coords$linetype[i],
        #                              lineend = "butt"
        ##                            )
        #                          )
         #                         
         ##                       }) -> gl
         #                       
         #                       grobs <- do.call(grid::gList, gl)
          #                      
          #                      ggname("geom_rrect", grid::grobTree(children = grobs))
          ##                      
          #                    }
                              
#)

ggname <- function(prefix, grob) {
  grob$name <- grid::grobName(grob, prefix)
  grob
}

library(ggplot2)

#####################
#####################


expressed_list<- read.csv("../results/expr_tx_list.csv")

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

results<- read.csv("../results/DE_all.csv")
dpsi_suppa <- 0.1
pval_suppa <- 0.05

files_suppa <- list.files(path = diff_suppa_path, pattern = ".dpsi", full.names = TRUE, include.dirs = TRUE, no.. = TRUE)
files_suppa <- files_suppa[-8]
names(files_suppa) <- str_extract(files_suppa, "...dpsi") %>% gsub(".dpsi", "", .)

df <- data.frame()

for (i in files_suppa) {
  event <- str_extract(i, "...dpsi") %>% gsub(".dpsi", "", .)
  data <- read.table(i, header = TRUE, sep="\t", strip.white=TRUE, na.strings="NA")
  colnames(data) <- c("X","dPSI", "pval")
  sig_up_in_KO <- dplyr::filter(data, dPSI < -dpsi_suppa) %>% dplyr::filter(pval < pval_suppa)
  sig_down_in_KO <- dplyr::filter(data, dPSI > dpsi_suppa) %>% dplyr::filter(pval < pval_suppa)
  df <- rbind(df, c(event, nrow(sig_up_in_KO), -nrow(sig_down_in_KO)))
}

colnames(df) <- c("event","up", "down")
df$event <- factor(df$event, levels=c("SE", "MX", "A5", "A3", "RI", "AF", "AL"))

suppa_up <- ggplot(df[,1:2], aes(x=event, y=as.numeric(up))) + geom_bar(stat="identity", width=0.5, color="black", fill="#d5e6ff") + theme_minimal() + theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=8), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank()) + coord_flip()
suppa_down <- ggplot(df[,1:3], aes(x=event, y=-as.numeric(down))) + geom_bar(stat="identity", width=0.5, color="black", fill="#d5e6ff") + theme_minimal() + theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=8), axis.title.y = element_blank(), axis.title.x = element_blank()) + coord_flip() + scale_x_discrete(position = "top") + scale_y_reverse()

dfx <- data.frame()
for (i in files_suppa) {
  event <- str_extract(i, "...dpsi") %>% gsub(".dpsi", "", .)
  data <- read.table(i, header = TRUE, sep="\t", strip.white=TRUE, na.strings="NA") %>% mutate(event = event)
  colnames(data) <- c("X","dPSI", "pval", "event")
  data<- data %>% filter(pval< 0.05) %>%mutate(diff= ifelse(dPSI < - dpsi_suppa , "up", ifelse(dPSI > dpsi_suppa,"down","no"))) %>% 
    filter(diff != "no") %>% group_by(event, diff) %>% summarise(n= n()) %>% mutate(n = if_else(diff == "down",-n, n))
  dfx <- rbind(dfx, data)
}
misreg <-ggplot(dfx, aes(x=factor(event, levels = rev(event_order)), y=n, group=diff, fill= diff))+ geom_bar(stat= "identity", color="black", width=0.6)+ ylab("Number of events") +
   scale_y_continuous(breaks = seq(-150, 150, 50),labels = abs(seq(-150, 150, 50))) +
    scale_fill_manual(values= c(rgb(0.7,0.2,0.1,0.5),rgb(0.2,0.7,0.1,0.5))) + 
    theme_classic() +
    theme( panel.background = element_rect("white"), panel.grid = element_blank(),legend.title = element_blank(),
           axis.ticks= element_blank(),axis.text.x=element_text(vjust= 0, colour="black"), legend.position = "top",
           legend.key.size = unit(0.3, "cm"),
           axis.title.x = element_blank(),
           axis.line.x = element_blank())


# Start on the same basis as circadian analysis to be consistent #
   files <- list.files(path = salmon_path, pattern = ".sf", 
                       full.names = TRUE, include.dirs = FALSE, no.. = TRUE, recursive = TRUE)
   names(files) <- str_extract(files, "([CT|KO]\\d+_\\d)(\\.sf)$", group = 1) %>% gsub("T","CT",.) %>% gsub( "O", "KO", .)
   
   #tx2gene <- read.table("../biomart_annotation_98bis2.txt", sep = "\t", header=TRUE) %>%
   # rename_with(~ tolower(gsub(".", "_", .x, fixed = TRUE))) %>% relocate(transcript_id, gene_id) %>%
   # distinct(transcript_id, .keep_all = TRUE)
   
   
   annotation$gene_id<-  sub(':','loc',annotation$gene_id) #why not above?
   tx2gene<- annotation[c("transcript_id", "gene_id")]
   
   BPPARAM <-  BiocParallel::MulticoreParam(2)
   
   txi.salmon <- tximport(files, type = "salmon", txOut = TRUE,
                          countsFromAbundance = "no",
                          ignoreTxVersion = TRUE,
                          importer = function(x) read.delim(x))
   
   
   exp_design <- data.frame(batch = factor(str_extract(colnames(txi.salmon$abundance), "(?<=_)\\d")),
                            time = as.numeric(str_extract(colnames(txi.salmon$abundance), "\\d+")),
                            cond = factor(str_extract(colnames(txi.salmon$abundance), "^[A-Z]+")),
                            sample_id= colnames(txi.salmon$abundance))
   exp_design[exp_design$cond == "KO",]$batch = 2
   #rownames(exp_design) <- colnames(counts)
   exp_design$cond <- fct_recode(exp_design$cond, WT = "CT", "Per_KO" = "KO")
   
   exp_design$inphase = cos(2 * pi * exp_design$time / 24)
   exp_design$outphase = sin(2 * pi * exp_design$time / 24)
   
   counts <- magrittr::set_colnames(txi.salmon$counts, exp_design$sample_id)
   mode(counts) <- "integer"
   
   dds <- DESeqDataSetFromMatrix(countData = counts,
                                 colData = exp_design,
                                 design= ~ batch + cond + inphase + outphase)
   
   row.names(tx2gene) <- tx2gene$transcript_id
   
   rowData(dds) <- tx2gene[rownames(dds),]
   
   keep <- rownames(dds) %in% expressed_list$transcript_id
   
   dds <- dds[keep, ]
   
   dds <- DESeq(dds, reduced = ~batch + cond, test = "LRT", fitType = "local")
   
   abundance = varianceStabilizingTransformation(dds, fitType = "local")
   abundance_no_batch_effect <- limma::removeBatchEffect(assay(abundance),
                                                         batch = colData(dds)$batch,
                                                         design = model.matrix(~cond + inphase + outphase, data = colData(dds)))
   
   cts2 <- data.frame(gene_id = rowData(dds)$gene_id, feature_id = rowData(dds)$transcript_id, counts(dds)[, exp_design$time==20], check.names = FALSE)
   
   d <- DRIMSeq::dmDSdata(counts = cts2, samples = exp_design[exp_design$time==20, ])
   
   n <- 4
   n.small <- 2
   
   d <- DRIMSeq::dmFilter(d,
                          min_samps_feature_expr=n.small, min_feature_expr=10,
                          min_samps_feature_prop=n.small, min_feature_prop=0.1,
                          min_samps_gene_expr=n-n.small, min_gene_expr=10)
   
   sampleData = samples(d)
   
   sampleData$time <- factor(sampleData$time)
   dxko <- DEXSeqDataSet(countData = as.matrix(counts(d)[,-c(1,2)]),
                         sampleData = sampleData,
                         design =~sample_id + exon + batch:exon + cond:exon,
                         featureID = counts(d)$feature_id,
                         groupID = counts(d)$gene_id %>% sub(':','loc',.))
   
   dxko <- estimateSizeFactors(dxko)
   dxko <- estimateDispersions(dxko, quiet=TRUE, fitType = 'local', BPPARAM=BPPARAM)
   dxko <- testForDEU(dxko, reducedModel=~sample_id + exon + batch:exon, BPPARAM=BPPARAM)
   dxkor <- DEXSeqResults(dxko, independentFiltering=TRUE)
   
   qval <- DEXSeq::perGeneQValue(dxkor)
   
   dxko.g <- data.frame(gene=names(qval),qval)
   
   pScreen <- qval
   
   pConfirmation <- matrix(dxkor$pvalue, ncol = 1)
   dimnames(pConfirmation) <- list(dxkor$featureID, "transcript")
   
   stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                         pScreenAdjusted=TRUE, tx2gene=tx2gene)
   stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05, allowNA=TRUE)
   suppressWarnings({
     dexko.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                      onlySignificantGenes=TRUE)
   })
   
   ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
   
   gene_symbols <- getBM(attributes=c('mgi_symbol', 'ensembl_gene_id'),
                         filters = 'ensembl_gene_id',
                         values = dexko.padj$geneID,
                         mart = ensembl)
   
   df_ko <- dexko.padj %>%
     merge(gene_symbols, by.x = "geneID", by.y = "ensembl_gene_id", all.x=TRUE) %>%
     cbind(abundance_no_batch_effect[dexko.padj$txID,]) %>%
     group_by(geneID) %>%
     mutate(dummy_var = as.character(row_number())) %>%
     pivot_longer(starts_with(c("CT", "KO")), names_to = c("sample", "rep"), names_sep = "_") %>%
     mutate(time = as.numeric(sub("CT|KO", "", sample)))
   
   abundance_no_batch_effect_no_log <- 2^abundance_no_batch_effect
   
   df_ko_new <- 2^abundance_no_batch_effect[,c("CT20_1","CT20_2", "KO20_1", "KO20_2")] |>
     data.frame() |>
     rownames_to_column("txID") |>
     right_join(dexko.padj) |>
     merge(gene_symbols, by.x = "geneID", by.y = "ensembl_gene_id", all.x=TRUE)
   cts <- df_ko_new[,c("CT20_1","CT20_2", "KO20_1", "KO20_2")]
   gene.cts <- rowsum(cts, df_ko_new$geneID)
   total.cts <- gene.cts[match(df_ko_new$geneID, rownames(gene.cts)),]
   prop <- cts/ total.cts
   colnames(prop) <- paste0("prop_", colnames(prop))
   df_ko_new <- cbind(df_ko_new, prop, propSD = sqrt(rowVars(data.matrix(prop)))) |>
     arrange(desc(propSD)) 
   
   df_DTU_plot<- df_ko_new %>%
     rowwise() %>%
     mutate(CT_mean = round(mean(c_across(c('prop_CT20_1', 'prop_CT20_2')), na.rm=TRUE),2),
            KO_mean=round(mean(c_across(c('prop_KO20_1', 'prop_KO20_2')), na.rm=TRUE),2)) %>% 
     dplyr::select(c(geneID, mgi_symbol ,txID,CT_mean,KO_mean )) %>%
     group_by(geneID) %>% 
     mutate(number= as.character(row_number())) %>% pivot_longer(cols=c('CT_mean', 'KO_mean'),
                       names_to='phenotype',
                       values_to='prop')
   
   write.csv(df_ko_new, paste(list_dir, "/df_ko_DTU_wtko.csv", sep=""))
   
   interesting<- c("Abhd14a", "Acox2", "Ankrd33b","Comt","Cox16","Ligp1", "Eef1g", "Hyi", "Ddx5", "Egfr")
   
   dtu_plot_new<-ggplot(data = df_DTU_plot %>% filter(mgi_symbol %in% interesting) , 
                        aes(x = phenotype, y = prop, fill = number)) + 
     facet_wrap(~mgi_symbol) +
     geom_col() + 
     scale_fill +
     theme_minimal() +
     scale_x_discrete(labels=c("CT_mean" = "WT", "KO_mean" = "KO")) +
     theme(legend.position = "none") +
     theme(strip.text = element_text(size=7), panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(), legend.position = "none", 
           axis.line=element_line(size = 0.5, colour = "black"), axis.title.x=element_blank()) 

      
   
   df_ko_no_log <- dexko.padj %>%
     merge(gene_symbols, by.x = "geneID", by.y = "ensembl_gene_id", all.x=TRUE) %>%
     cbind(abundance_no_batch_effect_no_log[dexko.padj$txID,]) %>%
     group_by(geneID) %>%
     mutate(dummy_var = as.character(row_number())) %>%
     pivot_longer(starts_with(c("CT", "KO")), names_to = c("sample", "rep"), names_sep = "_") %>%
     mutate(time = as.numeric(sub("CT|KO", "", sample)))
   
   #selection <- filter(df_ko_no_log, mgi_symbol %in% c("Bdh1", "Ankrd33b" ,"Commd6","H2-K1","Gib1", "Mup12","Hyi", "Npm1", "Calr", "Arhgap26", "Egfr", "Herpud1", "Etfrf1",  "Gnai2", "Hsp90ab1", "Ivd", "Kng1", "Mcm10", "Nrp1", "Pigr", "Steap3", "Slco1b2", "Kng1", "Etfrf1", "Rrbp1", "Hnrnpdl", "Lifr", "Rnf128", "Tkfc", "Slc25a42", "Canx", "Glo1", "Mup15"))
   #selection_2 <-selection %>% filter(!(mgi_symbol %in% c("Hyi", "Npm1", "Canx", "Calr", "Mup15", "Glo1","Egfr", "Kng1", "Ivd", "Rrbp1") ))
   #DTU_WT20vsKO20 <- ggplot(selection_2[selection_2$time==20,], aes(fill=dummy_var, y=value, x=sample, col= dummy_var), alpha=0.7) + 
   #  geom_bar(position="fill", stat="identity", width=0.5) + 
  #   facet_wrap(~mgi_symbol, scales = "free_y", nrow=3) + 
  #   theme_bw(base_size = 10) + 
  #   theme_minimal() + 
  #   theme(strip.text = element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.line=element_line(size = 0.5, colour = "black"), axis.title.x=element_blank()) + 
  #   scale_fill + 
  #   annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) + 
  #   scale_y_continuous(breaks = scales::breaks_extended(n=3)) + 
  #   ylab("Fraction") + 
  #   scale_x_discrete(labels=c("CT20" = "WT", "KO20" = "KO")) + 
  #   scale_color

#separated caus emanual selection 
#DTU_WT20vsKO20_a <- ggplot(df_ko_no_log[df_ko_no_log$time==20,]%>% subset(., mgi_symbol %in%  c("Hyi")), aes(fill=dummy_var, y=value, x=sample, col= dummy_var), alpha=0.7) + geom_bar(position="fill", stat="identity", width=0.8) +facet_wrap(~mgi_symbol, scales = "free_y", nrow=3) + theme_bw(base_size = 10) + theme_minimal() + theme(strip.text = element_text(size=7), panel.grid.major = element_blank(),axis.title.y=element_text(size=6), axis.text= element_text(size=6), panel.grid.minor = element_blank(), legend.position = "none", axis.line=element_line(size = 0.5, colour = "black"), axis.title.x=element_blank()) + scale_fill + annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) + scale_y_continuous(breaks = scales::breaks_extended(n=3)) + ylab("Fraction") + scale_x_discrete(labels=c("CT20" = "WT", "KO20" = "KO")) + scale_color
#DTU_WT20vsKO20_b <- ggplot(df_ko_no_log[df_ko_no_log$time==20,]%>% subset(., mgi_symbol %in%  c("Npm1")), aes(fill=dummy_var, y=value, x=sample, col= dummy_var), alpha=0.7) + geom_bar(position="fill", stat="identity", width=0.8) +facet_wrap(~mgi_symbol, scales = "free_y", nrow=3) + theme_bw(base_size = 10) + theme_minimal() + theme(strip.text = element_text(size=7), panel.grid.major = element_blank(),axis.title.y=element_text(size=6), axis.text= element_text(size=6), panel.grid.minor = element_blank(), legend.position = "none", axis.line=element_line(size = 0.5, colour = "black"), axis.title.x=element_blank()) + scale_fill + annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) + scale_y_continuous(breaks = scales::breaks_extended(n=3)) + ylab("Fraction") + scale_x_discrete(labels=c("CT20" = "WT", "KO20" = "KO")) + scale_color
#DTU_WT20vsKO20_c <- ggplot(df_ko_no_log[df_ko_no_log$time==20,]%>% subset(., mgi_symbol %in%  c("Canx")), aes(fill=dummy_var, y=value, x=sample, col= dummy_var), alpha=0.7) + geom_bar(position="fill", stat="identity", width=0.8) +facet_wrap(~mgi_symbol, scales = "free_y", nrow=3) + theme_bw(base_size = 10) + theme_minimal() + theme(strip.text = element_text(size=7), panel.grid.major = element_blank(),axis.title.y=element_text(size=6), axis.text= element_text(size=6), panel.grid.minor = element_blank(), legend.position = "none", axis.line=element_line(size = 0.5, colour = "black"), axis.title.x=element_blank()) + scale_fill + annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) + scale_y_continuous(breaks = scales::breaks_extended(n=3)) + ylab("Fraction") + scale_x_discrete(labels=c("CT20" = "WT", "KO20" = "KO")) + scale_color
#DTU_WT20vsKO20_d <- ggplot(df_ko_no_log[df_ko_no_log$time==20,]%>% subset(., mgi_symbol %in%  c("Calr")), aes(fill=dummy_var, y=value, x=sample, col= dummy_var), alpha=0.7) + geom_bar(position="fill", stat="identity", width=0.8) +facet_wrap(~mgi_symbol, scales = "free_y", nrow=3) + theme_bw(base_size = 10) + theme_minimal() + theme(strip.text = element_text(size=7), panel.grid.major = element_blank(),axis.title.y=element_text(size=6), axis.text= element_text(size=6), panel.grid.minor = element_blank(), legend.position = "none", axis.line=element_line(size = 0.5, colour = "black"), axis.title.x=element_blank()) + scale_fill + annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) + scale_y_continuous(breaks = scales::breaks_extended(n=3)) + ylab("Fraction") + scale_x_discrete(labels=c("CT20" = "WT", "KO20" = "KO")) + scale_color


# Correlation splicing gene expression #
 ab <- as.data.frame(abundance_no_batch_effect_no_log) %>% dplyr::select(CT20_1, CT20_2, KO20_1, KO20_2) %>% mutate(transcript_id = row.names(.)) %>% merge(., tx2gene, by="transcript_id") %>% group_by(gene_id) %>% mutate(log2FC_gene = log2(sum(CT20_1+CT20_2)/sum(KO20_1 + KO20_2))) %>% ungroup()
 ab$transcript_id <- gsub("\\..", "", as.character(ab$transcript_id))
   
# create whole misreg event list:
   
   misreg_event_list <- data.frame()
   
   for (i in files_suppa) {
     data <- read.table(i, header = TRUE, sep="\t", strip.white=TRUE, na.strings="NA")
     colnames(data) <- c("X", "dPSI", "pval")
     
     data <- data %>% mutate(Event = .$X) %>% dplyr::filter(pval < pval_suppa) %>% dplyr::filter(dPSI > dpsi_suppa | dPSI < -dpsi_suppa) %>% dplyr::select(dPSI, Event)
     misreg_event_list <- rbind(misreg_event_list, data)
   }
   
   # get their associated transcripts:
   
   files_suppa_event <- list.files(path =suppa_path, pattern = ".ioe", full.names = TRUE, include.dirs = TRUE, no.. = TRUE)
   names(files_suppa_event) <- str_extract(files_suppa_event, ".._strict.ioe") %>% gsub("_strict.ioe", "", .)
   
   all_event_list <- data.frame()
   
   for (i in files_suppa_event) {
     data <- read.table(i, header = TRUE, sep="\t", strip.white=TRUE, na.strings="NA") %>% dplyr::filter(event_id %in% misreg_event_list$Event) %>% dplyr::select(event_id, alternative_transcripts, gene_id)
     all_event_list <- rbind(all_event_list, data)
   }
   
   colnames(all_event_list) <- c("Event", "alternative_transcripts", "gene_id")
   
   all_event_list <- merge(all_event_list, misreg_event_list, by.x = "Event")
   
   all_event_list2 <- all_event_list %>% separate(alternative_transcripts, as.character(seq(1,24)), ",") %>% pivot_longer(!c(Event, dPSI), names_to = "tx_nb", values_to = "transcript_id") %>% na.omit(.) %>% merge(., ab, by="transcript_id")
   
   correlation_exp_splicing <- ggplot(all_event_list2, aes(y= abs(dPSI), x= abs(log2FC_gene))) +
     geom_point( color="black", size=0.3) +
     geom_smooth(method=lm , color="#fd9997", fill="#d5e6ff", se=TRUE) + 
     theme_ipsum() + 
     stat_cor() + 
     theme_minimal() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom", axis.line=element_line(size = 0.5, colour = "black"), legend.title = element_blank())
   
   ## Upset misreg only ##
   
# modify upset script to add black lines around the plots:
#  trace(UpSetR:::Make_main_bar,edit=TRUE) #add in line 63, 100, 105 & 110: ', col = "black"' and change in line 74 'colour = Main_bar_data$color' into 'colour = "black"'
#   trace(UpSetR:::Make_size_plot,edit=TRUE) #change in line 30 : 'colour = sbar_color' into 'colour = "black"'
   
   library("ggvenn")
   venn_diffexp_splice <- ggvenn(list(diff_AS=misreg_T$gene_id, diffExp=results %>% filter(diffexp !="No") %>% pull(gene_id)), columns = c('diff_AS', 'diffExp'),
                     show_percentage = FALSE, fill_color = c("#434299", "#822d30"), 
                     fill_alpha = 0.3, auto_scale = FALSE, 
                     text_size = 5, set_name_size=3) 
   

   suppa_files <- c(AS3 = paste(suppa_path, "/transcriptome_ext_A3_strict.ioe", sep=""),
                    AS5 = paste(suppa_path, "/transcriptome_ext_A5_strict.ioe", sep=""),
                    AF = paste(suppa_path, "/transcriptome_ext_AF_strict.ioe", sep=""),
                    AL =paste(suppa_path, "/transcriptome_ext_AL_strict.ioe", sep=""),
                    MX = paste(suppa_path, "/transcriptome_ext_MX_strict.ioe", sep=""),
                    RI = paste(suppa_path, "/transcriptome_ext_RI_strict.ioe", sep=""),
                    ES = paste(suppa_path, "/transcriptome_ext_SE_strict.ioe", sep=""))
   
   expr_tx_list <- rownames(abundance_no_batch_effect)
   
   suppa_data <- imap(suppa_files, function(f, x) {
     read.table(f, header = TRUE) %>% rowwise() %>%
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
   
   ##input for upset
   upset_input<- data.frame(transcript_id= unique(suppa_data$tx_id)) %>%  mutate(class=case_when(grepl("^ENS", transcript_id) ~ "Annotated", 
                                                                                                 .default = "Novel"))
   
   for (AS_type in unique(suppa_data$event_type)) {
     # Create a new column with x value
     upset_input[[AS_type]] <- ifelse(upset_input$transcript_id %in% suppa_data$tx_id[suppa_data$event_type == AS_type], 1, 0)
   }
   
   #Upset_plot_splicing <- upset(filter(upset_input, transcript_id %in% all_event_list2$transcript_id), main.bar.color = "#d0ecb6", sets.x.label = "Total Number of Event", mainbar.y.label = "Number of Isoforms" , text.scale = c(1.3, 1.3, 1, 1, 2, 1), query.legend = "top", nsets = 7, sets.bar.color = "black", empty.intersections = "off", queries = list(list(query = elements, params = list("class", "Annotated"), color = "#d5e6ff", active = T, query.name = "Annotated")), show.numbers = "no") 
   #svg(paste(fig_dir,"/Upset_plot_splicing.svg", sep=""), height=5, width=10.5)
   ##Upset_plot_splicing
   #dev.off()
   
   ## dpsi scatter plots
   dpsi_files <- list.files(path = diff_suppa_path, pattern = ".dpsi", full.names = TRUE, include.dirs = TRUE, no.. = TRUE) %>% .[-8]
   psivec_files <- list.files(path = diff_suppa_path, pattern = ".psivec", full.names = TRUE, include.dirs = TRUE, no.. = TRUE) %>% .[-8]
   names(dpsi_files) <- str_extract(dpsi_files, ".{2}.dpsi") %>% gsub(".dpsi", "", .)
   names(psivec_files) <- str_extract(psivec_files, ".{2}.psivec") %>% gsub(".psivec", "", .)
   n=1
   
   ###K change> dpsi and psi column value###
   for (i in dpsi_files) {
     event <- str_extract(i, ".{2}.dpsi") %>% gsub(".dpsi", "", .)
     dpsi <- read.table(i, header = TRUE)
     psivec <-read.table(psivec_files[n], header = TRUE) %>% mutate(avgWT = rowMeans(.[,1:2]), avgKO = rowMeans(.[,3:4]), Dpsi = dpsi[,2], pval = dpsi[,3], color = ifelse(Dpsi < -0.1 & pval <0.05, rgb(0.7,0.2,0.1,0.5), ifelse(Dpsi > 0.1 & pval <0.05, rgb(0.2,0.7,0.1,0.5), "grey90")))
     scatter <- ggplot(psivec, aes(x=avgWT, y=avgKO)) + 
                geom_point(colour = psivec$color, size=0.3) +
                theme_classic() + 
                ggtitle(event)+ylab("PSI in KO")+
                xlab("PSI in WT") + 
                theme(aspect.ratio=1,axis.text = element_text(size = 7), axis.title.y=element_text(size=7), axis.ticks.y = element_blank(), axis.title.x=element_text(size=7), plot.title = element_text(hjust = 0.5, size=12), panel.border = element_rect(colour = "black", fill=NA, linewidth =1.5))
     assign(paste("scatter", event, sep = "_"), scatter)
     n=n+1
   }
   scatter<-plot_grid(scatter_SE+ scale_x_continuous(breaks=c(0,  0.5,  1), labels = c(0,  0.5,  1))+ scale_y_continuous(breaks=c(0,  0.5,  1), labels = c(0,  0.5,  1)), 
                      scatter_MX+theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + scale_x_continuous(breaks=c(0,  0.5,  1), labels = c(0, 0.5,  1)), 
                      scatter_RI+theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + scale_x_continuous(breaks=c(0, 0.5,1), labels = c(0, 0.5, 1)), 
                      scatter_A5+theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + scale_x_continuous(breaks=c(0, 0.5,  1), labels = c(0, 0.5, 1)), 
                      scatter_A3+theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + scale_x_continuous(breaks=c(0, 0.5,1), labels = c(0,  0.5,  1)),
                      scatter_AL+theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + scale_x_continuous(breaks=c(0, 0.5,  1), labels = c(0,  0.5, 1)),
                      scatter_AF+theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + scale_x_continuous(breaks=c(0, 0.5, 1), labels = c(0, 0.5, 1)), 
                      nrow = 1, rel_widths = c(1.3, rep(1,6)))
   
   write.table(psivec, paste(list_dir, "/scatter_Fig5A.txt", sep=""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
   
   ### correlation transcript expression and splicing ###
   ab <- as.data.frame(abundance_no_batch_effect_no_log) %>% dplyr::select(CT20_1, CT20_2, KO20_1, KO20_2) %>% mutate(log2FC = log2(((CT20_1+CT20_2)/2)/((KO20_1+KO20_2)/2))) %>% mutate(txID = row.names(.))
   ab$txID <- gsub("\\..", "", as.character(ab$txID))
   ## get transcripts associated with a given event :
   misreg_event_list <- data.frame()
   for (i in files_suppa) {
     data <- read.table(i, header = TRUE, sep="\t", strip.white=TRUE, na.strings="NA")
     colnames(data) <- c("X","dPSI", "pval")
     data <- data %>% mutate(Event = .$X) %>% dplyr::filter(pval < pval_suppa) %>% dplyr::filter(dPSI > dpsi_suppa | dPSI < -dpsi_suppa) %>% dplyr::select(dPSI, Event)
     misreg_event_list <- rbind(misreg_event_list, data)
   }
   # get their associated transcripts:
   files_suppa_event <- list.files(path = suppa_path, pattern = ".ioe", full.names = TRUE, include.dirs = TRUE, no.. = TRUE)
   names(files_suppa_event) <- str_extract(files_suppa_event, ".._strict.ioe") %>% gsub("_strict.ioe", "", .)
   all_event_list <- data.frame()
   for (i in files_suppa_event) {
     data <- read.table(i, header = TRUE, sep="\t", strip.white=TRUE, na.strings="NA") %>% dplyr::filter(event_id %in% misreg_event_list$Event) %>% dplyr::select(event_id, alternative_transcripts, gene_id)
     all_event_list <- rbind(all_event_list, data)
   }
   colnames(all_event_list) <- c("Event", "alternative_transcripts", "gene_id")
   all_event_list <- merge(all_event_list, misreg_event_list, by.x = "Event")
   all_event_list2 <- all_event_list %>% separate(alternative_transcripts, as.character(seq(1,24)), ",") %>% pivot_longer(!c(Event, dPSI), names_to = "tx_nb", values_to = "txID") %>% na.omit(.) %>% merge(., ab, by="txID")
   correlation_exp_splicing <- ggplot(all_event_list2, aes(y= abs(dPSI), x= abs(log2FC))) +
     geom_point( color="black", size=0.3) +
     geom_smooth(method=lm , color="#fd9997", fill="#d5e6ff", se=TRUE) + theme_ipsum() + stat_cor() + theme_minimal() + theme(panel.grid.major = element_blank(), text=element_text(size=8),panel.grid.minor = element_blank(), legend.position = "bottom", axis.line=element_line(size = 0.5, colour = "black"), legend.title = element_blank()) + xlab("log2(gene expression fold change)") + ylab("dPSI")
   
  write.csv(all_event_list2, paste(list_dir, "/all_event_list2_cor_dpsi_exp.csv", sep=""))
   
#  gtf <- rtracklayer::import(paste(flair_path,"/transcriptome_productivity.gtf", sep="")) %>% sortSeqlevels %>% sort %>% dplyr::as_tibble()  %>% filter(transcript_id %in%  selection[selection$time==20,]$txID) 
#  gtf$dummy<-selection$dummy_var[match(gtf$transcript_id, selection$txID)]
    
 #  gene_structure <- function(gene_name,...){
 #    gene_annotation_from_gtf <- gtf %>% 
 #      filter(gene_name == !!gene_name) %>%  
#       dplyr::select(seqnames, start, end, strand, type, gene_id, transcript_id, gene_name, dummy) %>% arrange(dummy)
#     gene_rescaled <- shorten_gaps(
 #      exons = filter(gene_annotation_from_gtf, type == "exon"), 
#       introns = to_intron(filter(gene_annotation_from_gtf, type == "exon"), "transcript_id"), 
 #      group_var = "transcript_id")
#     gene_rescaled_exons <- gene_rescaled %>% dplyr::filter(type == "exon") 
#     gene_rescaled_introns <- gene_rescaled %>% dplyr::filter(type == "intron") 
#     str<- ggplot(gene_rescaled_exons, aes(xstart = start, xend = end, y = transcript_id)) +
#       geom_range(aes(fill=dummy,height=0.3),color="black") + scale_fill+ scale_color_manual(values = "black")+
 #      geom_intron(data = gene_rescaled_introns, aes(strand=strand), 
   #                arrow=grid::arrow(ends = "last", length = grid::unit(0.005, "inches")),
  #                 color="black",
  #                 arrow.min.intron.length=2) +
  #     theme_void(base_size = base_size) +
  #     theme(legend.position = "none",...) +
  #     scale_color_brewer() 
  #   return(str)  
  #  }
    
   #x<-plot_grid(DTU_WT20vsKO20_a, gene_structure("Hyi"), DTU_WT20vsKO20_b, gene_structure("Npm1"), rel_widths = c(1,3,1,3),rel_heights = c(1,1,1,1) ,nrow=1)
   #y<-plot_grid(DTU_WT20vsKO20_c, gene_structure("Canx"), DTU_WT20vsKO20_d, gene_structure("Calr"), rel_widths = c(1,3,1,3),rel_heights = c(1,1,1,1) ,nrow=1)
   #z<- plot_grid(x,y, nrow=2)

###################    
#### Work Lies ####
###################
    
library(ggfortify)

df_pca <- data.frame()
    
for (i in psivec_files) {
    event <- str_extract(i, "...psivec") %>% gsub(".psivec", "", .)
    data <- read.table(i, header = TRUE, sep="\t", strip.white=TRUE, na.strings="NA") %>% setNames(c("WT_1", "WT_2", "KO_1", "KO_2")) %>% mutate(type = event)
    df_pca <- rbind(df_pca, data)
}

    
df_pca <- df_pca %>% drop_na()
pca_res <- prcomp(t(df_pca[,-5]), scale. = FALSE)
pca_psi <- autoplot(pca_res, label = TRUE, shape = FALSE, label.size = 4) +
  theme_classic() +
  scale_x_continuous(expand = c(0.2, 0.2)) +
  scale_y_continuous(expand = c(0.2, 0.2))

############################# 
########## Entropy ########## 
############################# 

files <- list.files(path = suppa_path, pattern = ".psi", full.names = TRUE, include.dirs = TRUE, no.. = TRUE)

psi_dataKO <- map(files, function(f) read.delim(f, row.names = 1) %>% rownames_to_column(var = "eventID")) %>%
  bind_rows %>% 
  dplyr::select(c(eventID, contains("20"))) %>%
  pivot_longer(contains("20"), names_to = "time", values_to = "psi") %>%
  mutate(phenotype = ifelse(time %in% c("CT20_1", "CT20_2"), "WT", "KO"),
         time = 20,
         inphase = cos(2*pi*time/24),
         outphase = sin(2*pi*time/24)) %>%
  separate_wider_delim(eventID, names = c("geneID","eventID"), delim = ";")  


entropy_df_KO <- psi_dataKO[,-6:-7] %>%
  group_by(geneID, eventID, phenotype) %>%
  mutate(psi_av = mean(psi)) %>%
  ungroup() %>%
  dplyr::select(-psi) %>% 
  filter(psi_av > 0) %>%
  distinct() %>%
  mutate(psi_x_log2psi = psi_av*log2(psi_av)) %>%
  group_by(geneID, phenotype) %>%
  mutate(entropy=-sum(psi_x_log2psi)) %>%
  ungroup() %>%
  filter(entropy > 0) %>%
  dplyr::select(c(geneID, phenotype, entropy)) %>%
  distinct()

entropy_boxplot <- ggplot(entropy_df_KO, aes(x=phenotype, y=entropy, group=phenotype)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  ylim(0,5) +
  stat_compare_means(method = "t.test")


entropy_boxplot2 <- entropy_df_KO %>%
  filter(geneID %in% all_event_list$gene_id) %>% 
  ggplot(aes(x=phenotype, y=entropy, group=phenotype)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  ylim(0,5) +
  stat_compare_means(method = "t.test") +
  ylab("entropy for mis-spliced genes")


# entropy_scatter <- pivot_wider(entropy_df_KO, names_from = phenotype, values_from = entropy) %>%
#   rename(ensembl_gene_id = geneID) %>%
#   left_join(., gene_symbols, by='ensembl_gene_id') %>%
#   ggplot(aes(x=WT, y=KO, label=mgi_symbol)) +
#   geom_point() +
#   theme_classic() +
#   geom_abline(slope=1, linetype = 'dashed') #+
#   #geom_text(hjust=0, vjust=0)
# entropy_scatter

entropy_candidates <- pivot_wider(entropy_df_KO, names_from = phenotype, values_from = entropy) %>%
  mutate(splicing = ifelse(geneID %in% all_event_list$gene_id, 'mis-spliced', 'ns')) %>%
  ggplot(aes(x=factor(splicing, levels = c('ns', 'mis-spliced')), y=WT)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  ylim(0,5) +
  stat_compare_means(method = "t.test", size=3) +
  ylab("Entropy in WT") +
  theme(axis.title.x = element_blank())


pivot_wider(entropy_df_KO, names_from = phenotype, values_from = entropy) %>% 
  setNames(c('gene_id', 'WT_entropy', 'KO_entropy')) %>% 
  left_join(., distinct(mart[,2:3]), by = 'gene_id') %>%
  write.table(., "../lists/Fig5D.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

######################################### 
########## examples structures ########## 
######################################### 

#Fig5G <- plot_grid(
#  plot_grid(ggplot(df_ko_no_log[df_ko_no_log$time==20,]%>% subset(., mgi_symbol %in%  c("Canx")), aes(fill=dummy_var, y=value, x=sample, col= dummy_var), alpha=0.7) + geom_bar(position="fill", stat="identity", width=0.8) +facet_wrap(~mgi_symbol, scales = "free_y", nrow=3) + theme_bw(base_size = 10) + theme_minimal() + theme(strip.text = element_text(size=7), panel.grid.major = element_blank(),axis.title.y=element_text(size=6), axis.text= element_text(size=6), panel.grid.minor = element_blank(), legend.position = "none", axis.line=element_line(size = 0.5, colour = "black"), axis.title.x=element_blank()) + scale_fill + annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) + scale_y_continuous(breaks = scales::breaks_extended(n=3)) + ylab("Fraction") + scale_x_discrete(labels=c("CT20" = "WT", "KO20" = "KO")) + scale_color, gene_structure("Canx"), align = "h", axis = "bottom", rel_widths = c(1,2)),
#  plot_grid(ggplot(df_ko_no_log[df_ko_no_log$time==20,]%>% subset(., mgi_symbol %in%  c("Mup15")), aes(fill=dummy_var, y=value, x=sample, col= dummy_var), alpha=0.7) + geom_bar(position="fill", stat="identity", width=0.8) +facet_wrap(~mgi_symbol, scales = "free_y", nrow=3) + theme_bw(base_size = 10) + theme_minimal() + theme(strip.text = element_text(size=7), panel.grid.major = element_blank(),axis.title.y=element_text(size=6), axis.text= element_text(size=6), panel.grid.minor = element_blank(), legend.position = "none", axis.line=element_line(size = 0.5, colour = "black"), axis.title.x=element_blank()) + scale_fill + annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) + scale_y_continuous(breaks = scales::breaks_extended(n=3)) + ylab("Fraction") + scale_x_discrete(labels=c("CT20" = "WT", "KO20" = "KO")) + scale_color, gene_structure("Mup15"), align = "h", axis = "bottom", rel_widths = c(1,2)),
#  plot_grid(ggplot(df_ko_no_log[df_ko_no_log$time==20,]%>% subset(., mgi_symbol %in%  c("Glo1")), aes(fill=dummy_var, y=value, x=sample, col= dummy_var), alpha=0.7) + geom_bar(position="fill", stat="identity", width=0.8) +facet_wrap(~mgi_symbol, scales = "free_y", nrow=3) + theme_bw(base_size = 10) + theme_minimal() + theme(strip.text = element_text(size=7), panel.grid.major = element_blank(),axis.title.y=element_text(size=6), axis.text= element_text(size=6), panel.grid.minor = element_blank(), legend.position = "none", axis.line=element_line(size = 0.5, colour = "black"), axis.title.x=element_blank()) + scale_fill + annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) + scale_y_continuous(breaks = scales::breaks_extended(n=3)) + ylab("Fraction") + scale_x_discrete(labels=c("CT20" = "WT", "KO20" = "KO")) + scale_color, gene_structure("Glo1"), align = "h", axis = "bottom", rel_widths = c(1,2)),
#  ncol=1)

############################################################ 
########## Differential AS events mutual influence ########## 
############################################################ 

x<-suppa_data %>% filter(event_id %in% all_event_list2$Event)   # %>% group_by(tx_id) %>% summarise(x= n())  #misregulated events
#from the transcript event file, all the transcripts that undergo misreg > transcripts with all evetns  (or x$tx_id)
y<-suppa_data %>% filter(tx_id %in% x$tx_id) %>% mutate(misreg= ifelse(event_id %in% x$event_id, "yes", "no")) %>% 
  group_by(tx_id, misreg) %>% summarise(count_misreg=n()) #transcripts having misreg or not events
z<- y %>% group_by(tx_id) %>% mutate(count_all= sum(count_misreg)) %>% 
  filter(misreg=="yes") %>% mutate(frac_misreg= count_misreg/count_all) %>% ungroup() %>%
  group_by(frac_misreg) %>% mutate(count= n()) 
plot1<-ggplot(z, aes(x=frac_misreg, y=count_misreg))+ geom_point() + theme_bw() #put size for geom_point based on number of time the frac is coming if want with the size
z_without1<- z %>% filter(count_all !=1) 

Fig5E <- ggplot(z_without1 %>% group_by(count_misreg, frac_misreg) %>% mutate(count=n()) %>% ungroup()) +
  geom_point(aes(x= frac_misreg, y=count_misreg, size = count, alpha=count)) +
  geom_smooth(aes(x= frac_misreg, y=count_misreg), method=lm , color="#fd9997", fill="#d5e6ff", se=TRUE) +
  stat_cor(aes(x= frac_misreg, y=count_misreg), size=3) +
  theme_classic() +
  theme(axis.title = element_text(size=7)) +
  scale_size(range = c(2,8)) +
  ylab("Number of differential \nAS event per transcript") +
  xlab("Fraction of differential AS event per transcript")


############################################################# 
########## Differential AS events mutual influence ########## 
#############################################################
library(reshape2) 
library(magrittr)

colours<- c("ES"="#c0c0c0","ES+2 or more"="#000000", "ES+AF"="#94c4c1", "ES+AL"="#f9c7be", "ES+AS3"="#629183", "ES+AS5"="#ce6a6c", "ES+MX"="#e93237", "ES+RI"="#e79288",
            "AL"="#c0c0c0","AL+2 or more"="#000000", "AL+AF"="#94c4c1", "AL+AS3"="#629183", "AL+AS5"="#ce6a6c", "AL+MX"="#e93237", "AL+RI"="#e79288","AL+ES"="#42899b",
            "AF"="#c0c0c0","AF+2 or more"="#000000", "AF+AS3"="#629183", "AF+AS5"="#ce6a6c", "AF+MX"="#e93237", "AF+RI"="#e79288", "AF+ES"="#42899b",
            "MX"="#c0c0c0","MX+2 or more"="#000000", "MX+AF"="#94c4c1", "MX+AS3"="#629183", "MX+AS5"="#ce6a6c",  "AL+RI"="#e79288", "MX+ES"="#42899b",
            "AS3"="#c0c0c0","AS3+2 or more"="#000000", "AS3+AF"="#94c4c1",  "AS3+AS5"="#ce6a6c", "AS3+MX"="#e93237", "AS3+RI"="#e79288", "AS3+ES"="#42899b",
            "AS5"="#c0c0c0","AS5+2 or more"="#000000", "AS5+AF"="#94c4c1", "AS5+AS3"="#629183","AS5+MX"="#e93237", "AS5+RI"="#e79288","AS5+ES"="#42899b",
            "RI"="#c0c0c0","RI+2 or more"="#000000", "RI+AF"="#94c4c1","RI+AS3"="#629183", "RI+AS5"="#ce6a6c", "RI+MX"="#e93237", "RI+ES"="#42899b",
            "RI+1"="#7a306c","AF+1"="#7a306c","AS3+1"= "#7a306c", "AS5+1"="#7a306c", "MX+1"="#7a306c", "ES+1"="#7a306c","AL+1"="#7a306c"
)
#from the transcript event file, get the misregulated events
misreg_T<-suppa_data %>% filter(event_id %in% all_event_list2$Event)   # %>% group_by(tx_id) %>% summarise(x= n())  #misregulated events

#read the file that makes upset plot in figure1
upset_data_fig1<-read.csv("../data/upset_all.csv")  #for all 
upset_fig4<-upset_input %>% filter(transcript_id %in% misreg_T$tx_id)  #use upset_input %>% filter(transcript_id %in% all_event_list2$txID) if you want the other definition

##############################################################################################
combine_df<- data.frame(matrix(nrow = 0, ncol = 0)) 
chi_sqaure_pvalue<- data.frame(matrix(nrow = 0, ncol = 0)) 

for (x in c("AF", "AL", "AS3","AS5","ES", "MX", "RI") ){
  print(x)
  df <- upset_data_fig1 %>% dplyr::filter(upset_data_fig1[, x]==1)
  column9_values <- character(nrow(df))
  for (i in 1:nrow(df)) {
    cols_to_check <- names(df)[names(df) != "tx_id" & names(df) != x] #Exclude col1 and col8 from the columns to check
    indices <- which(df[i, cols_to_check] == 1)
    if (length(indices) > 0) {
      col_combination <- c(x, cols_to_check[indices])
      column9_values[i] <- paste0(col_combination, collapse = "+")
    } else {
      column9_values[i] <- x
    }
  }
  df$column9 <- column9_values
  column9_values[grep("^[A-Za-z0-9]{2,3}\\+[A-Za-z0-9]{2,3}$", column9_values)] <- paste0(x,"+1") #hash this line if you want for multiple categories 
  column9_values[grep("(\\+[^\\+]+){2,}", column9_values)] <- paste0(x,"+2 or more")
  df_pie <- data.frame(table(column9_values))
  # Calculate percentages
  df_pie$percentage <- (df_pie$Freq / sum(df_pie$Freq)) * 100
  # Sort the data frame by descending order of percentage
  df_pie_1 <- df_pie[order(-df_pie$percentage), ] %>% mutate(dataset="ALL", event_type=x)
  # Create the pie chart
  df <- upset_fig4 %>% filter(upset_fig4[, x]==1)
  
  ###for the KO dataset 
  column9_values <- character(nrow(df))
  for (i in 1:nrow(df)) {
    cols_to_check <- names(df)[names(df) != "tx_id" & names(df) != x] #Exclude col1 and col8 from the columns to check
    indices <- which(df[i, cols_to_check] == 1)
    if (length(indices) > 0) {
      col_combination <- c(x, cols_to_check[indices])
      column9_values[i] <- paste0(col_combination, collapse = "+")
    } else {
      column9_values[i] <- x
    }
  }
  df$column9 <- column9_values
  column9_values[grep("^[A-Za-z0-9]{2,3}\\+[A-Za-z0-9]{2,3}$", column9_values)] <- paste0(x,"+1") #hash this line if you want for multiple categories 
  column9_values[grep("(\\+[^\\+]+){2,}", column9_values)] <- paste0(x,"+2 or more")
  df_pie <- data.frame(table(column9_values))
  # Calculate percentages
  df_pie$percentage <- (df_pie$Freq / sum(df_pie$Freq)) * 100
  # Sort the data frame by descending order of percentage
  df_pie_2 <- df_pie[order(-df_pie$percentage), ] %>% mutate(dataset="KO", event_type=x)
  df_pie<- rbind(df_pie_1, df_pie_2)
  
  #do chisqaure analysis and save pvalue
  test<- df_pie %>% dplyr::select(dataset, column9_values, Freq) %>% dcast(., dataset ~ column9_values)
  rownames(test)<- test$dataset
  pval_df<- data.frame(event=x,p_val=chisq.test(test[,2:4])$p.value) %>% mutate(sig= ifelse(p_val<0.01, "yes", "no"))
  combine_df<- rbind(combine_df, df_pie)
  chi_sqaure_pvalue<- rbind(chi_sqaure_pvalue, pval_df)
}

combine_df <- combine_df %>%
  group_by(event_type, dataset) %>%
  arrange(column9_values) %>%
  mutate(label_y = cumsum(percentage)- 0.5 * percentage) %>%
  ungroup() %>%
  mutate(dataset = ifelse(dataset == "ALL", "Global", ifelse(dataset == "KO", "Diff.splice\nin KO", NA)),
         event_type =  ifelse(event_type == "ES", "SE", event_type))

plot_eventsx_all <- combine_df %>%
  ggplot(aes(x = factor(dataset, levels=c("Global", "Diff.splice\nin KO")), y = percentage, fill = column9_values)) +
  geom_bar(stat = "identity", width = 0.9, linewidth= 0.2, color="black", position = position_stack(reverse = TRUE)) +
  geom_text(aes(y = label_y, label = paste(round(percentage, digits = 1), "%", sep="")), colour = "white", size=2) +
  scale_fill_manual(values=colours) +
  facet_grid(~event_type) +
  theme_classic()  + labs(fill="AS event(s)") +
  theme(legend.key.size = unit(0.2, "cm")) +
  theme(axis.title.x = element_blank(), legend.position = "right",axis.text.x = element_text(angle = 25))
plot_eventsx_all

plot_eventsx_select <- filter(combine_df, event_type %in% c("SE", "RI")) %>%
  ggplot(aes(y = percentage, x ="", fill = column9_values)) +
  geom_bar(stat = "identity") +
  #geom_text(aes(x = label_y, label = paste(round(percentage, digits = 1), "%", sep="")), fontface = "bold",  colour = "white", size=2) +
  scale_fill_manual(values=colours) +
  facet_grid(vars(dataset),vars(event_type) ) +
  theme_void() + 
  coord_polar("y", start=0)+
    theme(axis.text.y = element_text(size=6),axis.title.y = element_blank(), axis.line.y = element_blank())

write.csv(combine_df, paste(list_dir, "/4_combine_df_Pie_chart.csv", sep=""), row.names = FALSE)

##################
##### supp  ######
##################
nb_event_per_gene <- psi_dataKO %>%
  filter(!psi %in% c(0, 1, NaN)) %>%
  dplyr::select(geneID, eventID) %>%
  group_by(geneID)%>%
  mutate(nb_event = n()) %>%
  ungroup %>%
  dplyr::select(-eventID) %>%
  distinct() %>%
  mutate(misreg = ifelse(geneID %in% str_replace(misreg_event_list$Event, ";.*", ""), "misreg", "ns")) %>%
  ggplot(aes(x=misreg, y=nb_event)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(method = "t.test", size=4) +
  ylim(0,50) +
  theme_classic() +
  theme(axis.title.x = element_blank()) +
  ylab("number of event per gene")

# Plot distribution of difference of entropy per gene and see if it’s centered to 0
delta_entropy_distribution <- pivot_wider(entropy_df_KO, names_from = phenotype, values_from = entropy) %>%
  mutate(delta_entropy = KO-WT) %>% 
  ggplot(aes(x=delta_entropy)) +
  geom_histogram(binwidth = 0.01) +
  theme_classic() +
  xlab("Delta entropy per gene (KO - WT)") +
  ylab("Number of genes") +
  xlim(-1,1)

# Run GO from previous delta entropy
delta_entropy_df <- pivot_wider(entropy_df_KO, names_from = phenotype, values_from = entropy) %>%
  mutate(delta_entropy = KO-WT)

ego_delta_entropy_up <- enrichGO(gene = dplyr::filter(delta_entropy_df, delta_entropy > 0.2)$geneID, 
         keyType = 'ENSEMBL', 
         universe = delta_entropy_df$geneID, 
         OrgDb = org.Mm.eg.db, 
         ont = 'All', 
         pAdjustMethod = "BH", 
         pvalueCutoff  = 0.05, 
         qvalueCutoff = 0.05, 
         readable = TRUE)

ego_delta_entropy_down <- enrichGO(gene = dplyr::filter(delta_entropy_df, delta_entropy < -0.2)$geneID, 
                                 keyType = 'ENSEMBL', 
                                 universe = delta_entropy_df$geneID, 
                                 OrgDb = org.Mm.eg.db, 
                                 ont = 'All', 
                                 pAdjustMethod = "BH", 
                                 pvalueCutoff  = 0.05, 
                                 qvalueCutoff = 0.05, 
                                 readable = TRUE)

# correlation delta entropy and nb of isoforms

nb_isoforms_df <- read.csv("../results/expr_tx_list.csv") %>%
  group_by(gene_id) %>%
  mutate(nb_isoforms = n()) %>%
  ungroup() %>%
  dplyr::select(gene_id, nb_isoforms) %>%
  setNames(c('geneID', 'nb_isoforms')) %>%
  distinct() %>%
  left_join(delta_entropy_df, ., by='geneID')

nb_isoforms_df<- merge(nb_isoforms_df, annotation, by.x="geneID", by.y="gene_id", all.x=TRUE)

correlation_deltaentropy_nbiso_gene <- ggplot(nb_isoforms_df, aes(x=nb_isoforms, y=delta_entropy)) +
  geom_point(size=0.1) +
  theme_classic() + 
  stat_cor(method = "pearson") +
  ylab("delta entropy (KO-WT)") +
  xlab("number of isoform per gene") +
  theme(aspect.ratio = 1) + 
  #geom_text(aes(label= ifelse(abs(delta_entropy) > 1 ,mgi_symbol, "" )),size=3) + 
  theme(text= element_text(size=12))

# correlation entropy mis-splicing
entropy_df_KO2 <- filter(entropy_df_KO, phenotype == "WT")
whole_event_list <- data.frame()

for (i in files_suppa) {
  data <- read.table(i, header = TRUE, sep="\t", strip.white=TRUE, na.strings="NA")
  colnames(data) <- c("X", "dPSI", "pval")
  
  data <- data %>% mutate(Event = .$X) %>% dplyr::select(dPSI, Event)
  whole_event_list <- rbind(whole_event_list, data)
}

corr_dpsi_entropy_df <- whole_event_list %>%
  separate(Event, into = c('geneID', "other")) %>%
  dplyr::select(-other) %>%
  mutate(dPSI = abs(dPSI)) %>%
  slice_max(by = geneID, order_by = dPSI, n=1) %>%
  distinct() %>% 
  left_join(., entropy_df_KO2)


x <- corr_dpsi_entropy_df %>%
  ggplot(aes(x= dPSI, y=entropy)) +
  geom_point(color="black", size=0.3) +
  geom_smooth(method=lm , color="#fd9997", fill="#d5e6ff", se=TRUE) + 
  theme_ipsum() + 
  stat_cor() + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom", axis.line=element_line(size = 0.5, colour = "black"), legend.title = element_blank()) +
  ylab("entropy in WT") +
  xlab("max dPSI observed")


####################
##### Work #########
####################
ggsave(paste(fig_dir,"/4_FigS2-2.svg", sep=""), plot=x, width=3, height = 3)
ggsave(paste(fig_dir,"/4_correlation_exp_splicing.svg", sep="") ,plot=correlation_exp_splicing, width=3, height = 3)
ggsave(paste(fig_dir,"/4_delta_entropy_distribution.svg", sep=""), plot=delta_entropy_distribution, width=3, height = 3)
ggsave(paste(fig_dir,"/4_misreg.svg", sep=""), plot=misreg, width=2, height = 2)
ggsave(paste(fig_dir,"/4_scatter.svg", sep=""), plot=scatter, width=8, height = 5)
#ggsave(paste(fig_dir,"/4_DTU_WT_KO.svg", sep=""), plot=DTU_WT20vsKO20, width=5, height = 3)
ggsave(paste(fig_dir,"/4_plot_eventsx_all.svg", sep=""), plot=plot_eventsx_all, width=8, height = 3)
ggsave(paste(fig_dir,"/4_Fig3E.svg", sep=""), plot=Fig5E, width=4, height = 3)
ggsave(paste(fig_dir,"4_FigS3-6.svg", sep=""), plot=entropy_candidates, width=2, height = 3)
ggsave(paste(fig_dir,"/4_pie_SE_RI.svg", sep=""), plot=plot_eventsx_select, width=4, height = 4)
ggsave(paste(fig_dir,"/4_DTU_WTKO.svg", sep=""), plot= dtu_plot_new, width=5, height = 5)
ggsave(paste(fig_dir,"/4_msreg_splice.svg", sep=""), plot= misreg, width=2.5, height = 2.5)
ggsave(paste(fig_dir,"/4_venn_diffsplice_diffExp.svg", sep=""), plot= venn_diffexp_splice, width=1, height = 1)


ggsave(paste(fig_dir,"/4_correlation_deltaentropy_nbissso_gene.svg", sep=""), plot=correlation_deltaentropy_nbiso_gene, width=3, height = 3)

#######
# create the df
# first get the go terms for genes in our list
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
GO_df <- getBM(attributes = c('ensembl_gene_id', 'go_id', 'name_1006', 'definition_1006', 'go_linkage_type', 'namespace_1003'), filters = 'ensembl_gene_id', values = unique(psi_dataKO$geneID), mart = mart) %>%
  setNames(c('ensembl_gene_id', 'GO_term_accession', 'GO_term_name', 'GO_term_definition', 'GO_term_evidence_code', 'GO_domain'))

GO_df<-GO_df %>% filter(GO_term_name != "")

# keep all cellular component go terms that gather at list 50 genes 
nb_go_df <- data.frame()
for (i in unique(GO_df$GO_term_name)){
  nb_go_df <- rbind(nb_go_df, data.frame(GO = i, nb_genes = length(unique(filter(GO_df, GO_term_name == i)$ensembl_gene_id)))) #get the nb of genes involved in the given GO term
}

list_GO_min <- filter(nb_go_df, nb_genes > 50)$GO # list of GO that contains min 50 genes
GO_df_filtered <- data.frame()
for (i in list_GO_min){
  df <- filter(psi_dataKO[,-6:-7], geneID %in% filter(GO_df, GO_term_name == i)$ensembl_gene_id) %>% mutate(GO = i, GO_domain = unique(filter(GO_df, GO_term_name == i)$GO_domain))
  GO_df_filtered <- rbind(GO_df_filtered ,df)
}
GO_df_filtered <- rbind(GO_df_filtered %>% left_join(., nb_go_df, by = 'GO'), mutate(psi_dataKO[,-6:-7], GO = "All", nb_genes = length(unique(psi_dataKO$geneID)), GO_domain = "All")) 

#G0_df_entropy <- GO_df_filtered %>%
#  group_by(geneID, eventID, phenotype) %>%
#  mutate(psi_av = mean(psi)) %>%
#  ungroup() %>%
#  dplyr::select(-psi) %>% 
#  filter(psi_av > 0) %>%
#  distinct() %>%
#  mutate(psi_x_log2psi = psi_av*log2(psi_av)) %>%
#  group_by(GO, phenotype) %>%
#  mutate(entropy=-sum(psi_x_log2psi)) %>%
#  ungroup() %>%
#  group_by(GO,phenotype) %>%
#  mutate(median_GO = median(entropy)) %>%
#  #ungroup() %>% 
#  mutate(relative_entropy = median_GO/nb_genes) %>%
#  dplyr::select(GO_domain,GO, phenotype, relative_entropy, median_GO, nb_genes) %>%
#  distinct() 



### compute entropy and delta entropy for each gene in each GO category
G0_df_entropy_2 <-GO_df_filtered %>%
  group_by(geneID, eventID, phenotype) %>%
  mutate(psi_av = mean(psi)) %>%
  ungroup() %>%
  dplyr::select(-psi) %>% 
  filter(psi_av > 0) %>%
  distinct() %>%
  mutate(psi_x_log2psi = psi_av*log2(psi_av)) %>%
  group_by(GO_domain, GO, geneID, phenotype) %>%
  mutate(entropy=-sum(psi_x_log2psi)) %>%
  ungroup() %>% 
  dplyr::select(-time, -psi_av, -psi_x_log2psi) %>% 
  pivot_wider(names_from = phenotype, values_from = entropy) %>% 
  mutate(del_entropy = WT-KO) %>% 
  ungroup()

### Extract entropy for each gene (to use later to permutation p-value calculation)
G0_gene_entropy <- G0_df_entropy_2 %>% 
  distinct(geneID, .keep_all = TRUE) %>% 
  dplyr::select(geneID, del_entropy)

### Calculate the median entropy per GO
G0_df_del_entropy <-  G0_df_entropy_2 %>% 
  group_by(GO_domain, GO, nb_genes) %>% 
  summarize(WT = median(WT),
            KO = median(KO),
            del_entropy = median(del_entropy))

### Shuffle gene level entropy values to computate permuation median GO values and from that p-values
pvalue <- rep_len(0, length(G0_df_del_entropy$del_entropy))

ITER <- 1999
for (i in seq(ITER)){
  G0_gene_entropy_shuffle <- G0_gene_entropy %>% 
    mutate(del_entropy = sample(del_entropy))
  G0_df_del_entropy_back <- G0_df_entropy_2 %>% 
    dplyr::select(geneID, GO, GO_domain) %>% 
    left_join(G0_gene_entropy_shuffle, by=join_by(geneID)) %>% 
    group_by(GO_domain, GO) %>% 
    summarize(del_entropy = median(del_entropy), .groups = "drop")
  pvalue <- pvalue + (G0_df_del_entropy_back$del_entropy >= G0_df_del_entropy$del_entropy)
}

G0_df_del_entropy$pvalue = (pvalue+1)/(ITER+1) ###Note, p-value for All category is meaningless.
G0_df_del_entropy$qvalue = p.adjust(G0_df_del_entropy$pvalue, method="fdr")  ### Multiple testing correction

###K###
significant<- G0_df_del_entropy %>% filter(pvalue < 0.05 && qvalue < 0.05)  #filter by a q and p value to get significant
significant<- significant[order(-significant$del_entropy),]  #sort by decreasing delta_entropy
write.csv(G0_df_del_entropy, paste(list_dir, "/del_entropy_GO_WTKO.csv", sep=""), row.names = FALSE)

#remove manually the ones we do not want
dont_want<- c("photoreceptor inner segment", "sarcoplasmic reticulum", "dendritic shaft",
              "skeletal system development", "neuron projection morphogenesis", "axoneme", 
              "multicellular organism growth", "negative regulation of neuron differentiation", 
              "defense response to Gram-positive bacterium", "protein K48-linked ubiquitination" )
significant<- significant[1:30,] 
significant<- significant %>% filter(!(GO %in% dont_want)) %>% 
  mutate(GO = str_replace(GO, "activating", "activating\n"),
         str_replace(GO, "coupled", "coupled\n"),
         )
                                                                  
GO_df_plot<-ggplot(significant, aes(x = -log10(pvalue), y = GO,size = del_entropy))+
  geom_point(alpha = 0.7,aes(color=significant$nb_genes) ) + #geom point for each GO
  scale_colour_gradient(low = "#C24163", high="#3952F5") +
  xlab("-log10(pvalue)") +
  ylab("GO term") +
  geom_text(data = significant,  
            mapping = aes(x =  -log10(pvalue),
            y = GO,
            label = GO), 
            vjust = -1.2, #to bring them above the dots
            hjust = 0.2, #to move them a bit to left
            inherit.aes = FALSE,
            size=3
           ) + 
  scale_size(range = c(3,10)) +
  theme_bw() + 
  theme(axis.title.y=element_blank(), #Remove background panels 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        panel.grid= element_blank()) +labs(color="nb of genes", size="delta entropy")

ggsave(paste(fig_dir,"/4_GO_df_plot.svg", sep=""), plot=GO_df_plot, width=5, height = 5)

#interesting <- c("All", "insulin receptor signaling pathway", "cellular response to insulin stimulus", "cell cycle","cholesterol metabolic process", "glucose homeostasis")
#G0_df_entropy_table <- dplyr::select(G0_df_entropy, c(GO_domain, GO,relative_entropy, phenotype, nb_genes)) %>% distinct() %>% setNames(c("GO_domain", "GO","relative_entropy","phenotype","no_genes"))
#delta_entropy_df <- pivot_wider(G0_df_entropy_table, names_from = phenotype, values_from = relative_entropy) %>%
#  mutate(delta_entropy = KO-WT, col= ifelse(GO %in% interesting, GO , "Others"))
#
#write.table(delta_entropy_df, "../delta_df_entropy_relative.tsv", row.names = FALSE)
#p_relative <- ggplot(delta_entropy_df, aes(x= delta_entropy)) + 
#  geom_density(aes(y=..count.., fill= "#800204", col="#6B2425" ),show_guide=FALSE, alpha=0.2) +
#  #geom_density(aes(fill= col),bins=150) +
#  xlab("delta Relative entropy") +
#  geom_text(aes(label=ifelse(abs(delta_entropy) >= 0.03, GO,  ""), y=10,angle = 90))+
#  geom_rug(alpha = 1/2, aes(y = 0), size=0.5)+
#  theme_bw() + 
#  theme(panel.background=element_blank(), panel.grid= element_blank()) 

#ggsave(paste(fig_dir,"/4_entropy_different_WTKO.svg", sep=""), plot=p_relative, width=4, height = 4)

################################################################################
#G0_df_entropy_table <- dplyr::select(G0_df_entropy, c(GO_domain, GO,median_GO, phenotype, nb_genes)) %>% distinct() %>% setNames(c("GO_domain", "GO","median_GO","phenotype","no_genes"))
#delta_entropy_df <- pivot_wider(G0_df_entropy_table, 
#                                names_from = phenotype, 
#                                values_from = median_GO) %>%
#                    mutate(delta_entropy = KO-WT, 
#                          col= ifelse(GO %in% interesting, GO , "Others"))

#p_median <- ggplot(delta_entropy_df, aes(x= delta_entropy)) + 
#  #geom_density(aes(y=..count.., fill= col), alpha=0.2)
#  geom_density(aes(fill= col),bins=150) +
#  xlab("delta Median_GO entropy") + theme_bw() +
#  geom_text()
#  geom_rug(alpha = 1/2)
#range(delta_entropy_df$delta_entropy) 
################################################################################

#gene_GO_insulin<-GO_df %>% filter(GO_term_name== "cellular response to insulin stimulus") %>% select(ensembl_gene_id) %>% unique()
#mis_Reg_list<- misreg_event_list %>% mutate(gene= gsub(";.*$","", misreg_event_list$Event))
#common<-intersect(gene_GO_insulin$ensembl_gene_id,mis_Reg_list$gene)
#small_pathway<- mis_Reg_list %>% filter(gene %in% common)
#write.table(small_pathway, "../cellular_resopnse_genes.tsv", row.names= FALSE)
#########################
#load clip-data
#srsf6<- read.csv("../data/ENCORI_mm10_RBP-RNA_Srsf6_all.csv")
#srsf1<- read.csv("../data/ENCORI_mm10_RBP-RNA_Srsf1_all.csv")
#srsf1_common<- srsf1 %>% filter(geneID %in% intersect(insulin_genes$ensembl_gene_id, srsf1$geneID)) %>% select(geneName)
#srsf6_common<- srsf6 %>% filter(geneID %in% intersect(insulin_genes$ensembl_gene_id, srsf6$geneID)) %>% select(geneName)
#intersect(srsf1_common$geneName, srsf6_common$geneName)

#entropy_df_KO <- psi_dataKO[,-6:-7] %>%
#  group_by(geneID, eventID, phenotype) %>%
#  mutate(psi_av = mean(psi)) %>%
#  ungroup() %>%
#  dplyr::select(-psi) %>% 
#  filter(psi_av > 0) %>%
#  distinct() %>%
#  mutate(psi_x_log2psi = psi_av*log2(psi_av)) %>%
#  group_by(geneID, phenotype) %>%
#  mutate(entropy=-sum(psi_x_log2psi)) %>%
#  ungroup() %>%
  #filter(entropy > 0) %>%
 # dplyr::select(c(geneID, phenotype, entropy)) %>%
 # distinct()
#insulin_genes<- GO_df_filtered %>% filter(GO=="insulin receptor signaling pathway")
#insulin_genes_delta<- pivot_wider(entropy_df_KO %>% filter(geneID %in% insulin_genes$geneID), names_from = phenotype, values_from = entropy) %>%
#                        mutate(delta_entropy = KO-WT) 
#plot<- insulin_genes_delta %>% ggplot(aes(x=delta_entropy)) +
#  geom_histogram(binwidth = 0.03) +
#  theme_classic() +
#  xlab("Delta entropy per gene (KO - WT)") +
#  ylab("Number of genes")

#ggsave(paste(fig_dir,"/gene_wise_entropy_insulin.svg", sep=""), plot=plot, width=5, height = 5)
#write.csv(insulin_genes_delta, "../extra/insuling_genes_delta.csv", row.names = FALSE)