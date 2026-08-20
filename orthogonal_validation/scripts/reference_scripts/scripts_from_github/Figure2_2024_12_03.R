#setwd and load packages
setwd("C:/Users/kmamgain/prettycoolPores_new/scripts")

library(magrittr)
library(tximport)
library(DESeq2)
library(stageR)
library(ggnewscale)
library(DEXSeq)
library(viridis)
library(GenomicFeatures)
library(DRIMSeq)
library(scales)
library(reshape2)
library(patchwork)
library(ggmosaic)
library(tidyverse) 
library(lemon)
library(ggthemes)
library(cowplot)
library(stringr)
library(ggplot2)
library(extrafont)
library(ggtranscript)
library(svglite)

fdr_cutoff <- 0.05
base_size<-8
retrieve_from_biomart <- FALSE
tx_names_with_version <- TRUE
class_names <- c("Annotated", "Novel-Annotated", "Novel-Novel")
expressed_t<-read.csv("../results/expr_tx_list.csv")

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

##-------------Functions-----------------------------------------------##
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

########

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

my_pretty_plot = function(data, theme_scale, te_families, filter, plotted_col) {
  data$famSub = te_families[match(data$te_name,
                                  te_families$te_name), "family"]
  tmp = str_split(data$famSub, "/", simplify = TRUE)
  data$fam = tmp[, 1]
  data$subfam = tmp[, 2]
  data[is.na(data$famSub), "famSub"] = "merged"
  data[is.na(data$subfam), "subfam"] = "merged"
  data[is.na(data$fam), "fam"] = "merged"
  data[data$subfam == "", "subfam"] = data[data$subfam == "", "famSub"]
  
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
  statuses = statuses[order(statuses)]
  
  data$matching_transcript =
    factor(x = data$matching_transcript,
           levels = unique(
             data[order(data$fam,
                        data[, plotted_col],
                        data$trcov),
                  "matching_transcript"]))
  
  widths = sapply(statuses, function(x) {
    nrow(data[data$status == x, ])
  })  / nrow(data)
  
  theme_plot = theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     panel.grid = element_blank(),
                     text = element_text(size = 15),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank())
  
  discrete_factor = 1 / 2
  data$disc_trcov = as.factor(round(data$trcov / discrete_factor)
                              * discrete_factor)
  colors = c("#ffc9fc", "red", "#ac00ac")
  names(colors) = c(0, 0.5, 1)
  
  labels = c("Low", "Half-half", "High")
  names(labels) = c(0, 0.5, 1)
  labels = as_labeller(labels)
  
  main_plots = list()
  titles = list()
  for (x in statuses) {
    tmp_plot = ggplot(data = data[data$status == x, ],
                      aes(x = matching_transcript,
                          y = fam, fill = disc_trcov)) +
      geom_tile() +
      scale_fill_discrete(type = colors,
                          labels = labels,
                          drop = FALSE) +
      xlab(NULL) + ylab(NULL) + theme_poster +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid = element_blank(),
            legend.position = "right",
            text = element_text(size = 15)) +
      scale_y_discrete(na.translate = FALSE, drop = FALSE) +
      guides(fill = guide_legend(title = "Transcript TE content")) +
      ggtitle(x)
    if (x == statuses[1]) {
      y_axis = get_y_axis(tmp_plot)
    }
    tmp_plot = tmp_plot + theme(axis.text.y = element_blank())
    if (x == statuses[length(statuses)]) {
      legend = get_legend(tmp_plot)
    }
    tmp_plot = tmp_plot + theme(legend.position = "none")
    titles[[which(statuses == x)]] = get_title(tmp_plot)
    tmp_plot = tmp_plot + ggtitle(NULL)
    main_plots[[which(statuses == x)]] = tmp_plot
  }
  
  lower_plots = list()
  for (x in statuses) {
    plotted_col = ensym(plotted_col)
    tmp_plot = ggplot(data = data[data$status == x, ],
                      aes(x = matching_transcript,
                          y = 1, fill = !!plotted_col)) +
      geom_tile() +
      theme_scale +
      xlab(NULL) + ylab(NULL) + theme_poster +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid = element_blank(),
            legend.position = "bottom",
            text = element_text(size = 15),
            axis.text.y = element_blank())
    if (x == statuses[length(statuses)]) {
      legend2 = get_legend(tmp_plot)
    }
    tmp_plot = tmp_plot + theme(legend.position = "none")
    lower_plots[[which(statuses == x)]] = tmp_plot
  }
  
  grid3 = plot_grid(plotlist = titles, nrow = 1, rel_widths = widths)
  grid1 = plot_grid(plotlist = main_plots, nrow = 1, rel_widths = widths)
  grid2 = plot_grid(plotlist = lower_plots, nrow = 1, rel_widths = widths)
  final_grid = plot_grid(NULL, grid3, NULL,
                         y_axis, grid1, legend,
                         NULL, grid2, legend2,
                         rel_heights = c(1, 9, 1),
                         rel_widths = c(1, 10, 3),
                         ncol = 3, nrow = 3)
  return(final_grid)
}

###--------------------------------------------------------###

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


files <- list.files(path = salmon_path, pattern = "quant_CT", 
                    full.names = TRUE, include.dirs = FALSE, no.. = TRUE, recursive = TRUE)
names(files) <- str_extract(files, "([CT|KO]\\d+_\\d)(\\.sf)$", group = 1) %>% gsub("T","CT",.) 

txi_nanopore_wt <- tximport(files, type = "salmon", txOut = TRUE, countsFromAbundance = "no", ignoreTxVersion = TRUE, 
                            importer = function(x) read.delim(x))

exp_design <- data.frame(batch = factor(str_extract(colnames(txi_nanopore_wt$abundance), "(?<=_)\\d")),
                         time = as.numeric(str_extract(colnames(txi_nanopore_wt$abundance), "\\d+")),
                         sample_id = colnames(txi_nanopore_wt$abundance)) %>%
  mutate(inphase = cos(2*pi*time/24),
         outphase = sin(2*pi*time/24))
rownames(exp_design) <- exp_design$sample_id


counts <- magrittr::set_colnames(txi_nanopore_wt$counts, exp_design$sample_id)
if (tx_names_with_version){
  rownames(counts) <- str_extract(rownames(counts), "([A-Za-z0-9-]+)")
}
mode(counts) <- "integer"
##---------Circadian analysis of transcriptome--------------------##

circadian_dds <- DESeqDataSetFromMatrix(countData = counts,
                                        colData = exp_design,
                                        design= ~ batch + inphase + outphase)

rowData(circadian_dds) <- annotation[rownames(circadian_dds),]

keep <- rownames(counts) %in% expressed_t$transcript_id

circadian_dds <- circadian_dds[keep, ]

circadian_dds <- DESeq(circadian_dds, reduced = ~batch, test = "LRT", fitType = "local")

res <- DESeq2::results(circadian_dds, independentFiltering = TRUE, alpha = fdr_cutoff, tidy = TRUE)

res$groupID <- rowData(circadian_dds)$gene_id

res <- cbind(res, coef(circadian_dds))

pScreen <- perGeneQValueHack(res)
pConfirmation <- matrix(res$pvalue, ncol = 1)

dimnames(pConfirmation) <- list(res$row, "transcript")

stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=TRUE, tx2gene=annotation[c("transcript_id", "gene_id")])

stageRObj <- stageWiseAdjustment(stageRObj, method="dte", alpha=fdr_cutoff, allowNA=TRUE)

suppressWarnings({
  circadian_analysis <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=FALSE)
})


anno_bed<-read.table(paste(flair_path,"/productivity.txt", sep=""), sep="\t") %>% 
  mutate(transcript_id= .$V4, productivity_flair= .$V13) %>% left_join(., annotation) %>% 
  mutate(productivity_flair = fct_relevel(productivity_flair, c("PRO","PTC","NGO","NST")),
         productivity_flair = fct_recode(productivity_flair, 
                                         !!!c("productive" = "PRO",
                                              "premature termination codon" = "PTC",
                                              "no start codon" = "NGO",
                                              "has start but no stop codon" = "NST"))) %>% filter(transcript_id %in% expressed_t$transcript_id) %>% dplyr::select(c(1:13,16))
colnames(anno_bed) <- c("chrom", "chromStart", "chromEnd", "tx_name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "exons_nb", "exons_size", "exons_start",  "productivity","gene_name")

circadian_analysis %<>% cbind(coef(circadian_dds)) %>%
  cbind(rowData(circadian_dds)[,1:6]) %>%
  dplyr::select(-geneID, -txID) %>%
  filter(gene < fdr_cutoff) %>%
  mutate(phase = 12/pi * (atan2(outphase, inphase) %% (2*pi)),
         amp = sqrt(inphase^2 + outphase^2)) %>%
  relocate(transcript_id, gene_id) %>%
  group_by(gene_id) %>%
  mutate(n_tx = n(),
         tx_name = sub("\\.\\d+","", transcript_id)) %>%
  left_join(anno_bed[c("tx_name","productivity")])

mosaic_annotation <- ggplot(data = circadian_analysis %>%
                    mutate(rhythmic = ifelse(transcript< fdr_cutoff, "Rhythmic", "Non-rhythmic"), 
                           annot = factor(ifelse(grepl("^ENS", transcript_id), "Annotated", "Novel"), 
                                          levels = c("Annotated", "Novel")))) + 
  geom_mosaic(aes(x=product(rhythmic, annot), fill=annot), offset = 0.05) + 
  scale_fill_manual(values = c("#66c2a5","grey20")) + xlab("") + ylab("") + theme_classic(base_size = base_size) +
  theme(legend.position = "none", aspect.ratio = 1)  + theme(text = element_text(size=base_size+2)) #+ guides(x=guide_axis(n.dodge = 2))

mosaic_prod <- ggplot(data = circadian_analysis %>% mutate(rhythmic = ifelse(transcript< fdr_cutoff, "Rhythmic", "Non-Rhythmic"), 
                                                           prod = factor(ifelse(grepl("PRO", productivity), "Productive", "Non-productive"), levels = c("Productive", "Non-productive")))) + 
  geom_mosaic(aes(x=product(rhythmic, prod), fill=prod), offset = 0.05) + 
  scale_fill_colorblind() + xlab("") + ylab("") + theme_classic(base_size = base_size) +
  theme(legend.position = "none", aspect.ratio = 1, axis.title.x= element_blank()) + theme(text = element_text(size=base_size+2))  #+ guides(x=guide_axis(n.dodge = 2))


xtabs(~ (transcript<fdr_cutoff) + (grepl("^ENS", transcript_id)), data=circadian_analysis, n_tx>1) %>% chisq.test()

xtabs(~ (transcript<fdr_cutoff) + (grepl("PRO", productivity)), data=circadian_analysis, n_tx>1) %>% chisq.test()

abundance = varianceStabilizingTransformation(circadian_dds, fitType = "local")
abundance_no_batch_effect <- limma::removeBatchEffect(assay(abundance),
                                                      batch = colData(circadian_dds)$batch,
                                                      design = model.matrix(~inphase + outphase, data = colData(circadian_dds)))
tx2gene=annotation[c("transcript_id", "gene_id")]

phase_distribution1 <- ggplot(circadian_analysis %>% filter(n_tx>1 & transcript < fdr_cutoff) %>% 
                               mutate(fill = ifelse(productivity=="PRO", "Productive", "Non-productive"),
                                      fill = factor(fill, levels=c("Productive", "Non-productive"))) %>%
                               filter(!is.na(fill)) %>% filter(fill=="Productive")) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 12,
            fill = '#FAF8F8')+
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 12, ymax = 24,
            fill = '#D6D6D6', alpha = 0.05)+
  geom_histogram(aes(y=phase, fill=fill), binwidth = 1, center=0.5, color="white", position="dodge") +
  coord_polar(theta="y", start=0) + 
  scale_y_continuous(breaks = seq(0,24,2), limits = c(0,24), expand = expansion()) + 
  scale_fill_colorblind(breaks=c("Productive", "Non-productive"), name = "") + 
  theme_void(base_size = base_size) +
  geom_vline(aes(xintercept = x), data.frame(x = seq(3)*12), color = "#928E8E", linewidth=0.5) +
  geom_vline(aes(xintercept = 48)) +
  geom_hline(aes(yintercept = y), data.frame(y = seq(0,24,2)), color="#928E8E", linewidth=0.5) +
  theme(axis.text.x = element_text(face="plain", size=rel(0.7)), strip.text  = element_blank(), legend.position = "none", aspect.ratio = 1) + 
  scale_x_continuous(expand = expansion(), limits = c(0,52))
  
phase_distribution2 <- ggplot(circadian_analysis %>% filter(n_tx>1 & transcript < fdr_cutoff) %>% 
                                 mutate(fill = ifelse(productivity=="PRO", "Productive", "Non-productive"),
                                        fill = factor(fill, levels=c("Productive", "Non-productive"))) %>%
                                 filter(!is.na(fill)) %>% filter(fill=="Non-productive")) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 12,
            fill = '#FAF8F8')+
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 12, ymax = 24,
            fill = '#D6D6D6', alpha = 0.05)+
  geom_histogram(aes(y=phase), binwidth = 1, center=0.5, fill="#E69F00", color="white", position="dodge") +
  coord_polar(theta="y", start=0) + 
  scale_y_continuous(breaks = seq(0,24,2), limits = c(0,24), expand = expansion()) + 
  theme_void(base_size = base_size) +
  geom_vline(aes(xintercept = x), data.frame(x = seq(1)*12), color = "#928E8E", linewidth=0.5) +
  geom_vline(aes(xintercept = 15)) +
  geom_hline(aes(yintercept = y), data.frame(y = seq(0,24,2)), color="#928E8E", linewidth=0.5) +
  theme(axis.text.x = element_text(face="plain", size=rel(0.7)), strip.text  = element_blank(), legend.position = "none", aspect.ratio = 1) + 
  scale_x_continuous(expand = expansion(), limits = c(0,16))
write.table(circadian_analysis, paste(list_dir, "/circadian_analysis.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE)

##---------------Heat maps of cycling transcripts-------------------------##
phase_sorted_tx <- circadian_analysis %>% 
  subset(gene < fdr_cutoff) %>%
  group_by(gene_id) %>%
  mutate(med_phase = median(phase),
         n_tx = n()) %>%
  filter(n_tx > 1) %>%
  arrange(desc(med_phase), desc(phase)) %>% 
  pull(transcript_id)

circadian_transcripts<-  circadian_analysis %>% 
  subset(transcript < fdr_cutoff) %>%
  group_by(gene_id) 

circadian_transcripts_gene_name<- merge(circadian_transcripts, annotation, by="transcript_id", all.x=TRUE, all.y=FALSE)
write.table(circadian_transcripts_gene_name, paste(list_dir, "/circadian_transcripts_gene_name.tsv", sep=""), quote = FALSE, sep = "\t", row.names = FALSE)
write.csv(circadian_transcripts, paste(list_dir, "/circadian_transcripts.csv", sep=""), row.names = FALSE)

df_1 <- abundance_no_batch_effect %>% 
  as.data.frame %>% rownames_to_column(var = "transcript_id")  %>%
  filter(transcript_id %in% phase_sorted_tx) %>%
  pivot_longer(-transcript_id, names_to = c("sample","rep"), names_sep="_", values_to = "tx_exp") %>%
  separate(sample, c("cond","time"), convert = TRUE, sep = 2) %>%
  group_by(transcript_id) %>%
  dplyr::select(-cond) %>%
  left_join(tx2gene) %>%
  group_by(gene_id, time, rep) %>%
  mutate(gene_exp = mean(tx_exp)) %>%
  group_by(gene_id) %>%
  mutate(gene_exp = scale(gene_exp)) %>%
  group_by(transcript_id) %>%
  mutate(tx_exp = scale(tx_exp)) %>%
  mutate(transcript_id = factor(transcript_id, levels = phase_sorted_tx, ordered = TRUE),
         time = factor(time, levels = seq(0, 20, 4), ordered = TRUE)) %>%
  arrange(transcript_id)

n_transcripts_rhy <- length(phase_sorted_tx)

LD_bars <- data.frame(xmin = c(rep(c(0.5, 3.5),2)),
                      xmax = c(rep(c(3.5, 6.5),2)),
                      col = c("grey50", "black"),
                      rep = c(1, 1, 2, 2))

n_genes <- circadian_analysis %>% group_by(gene_id) %>% mutate(n_tx = n()) %>% filter(n_tx>1 & gene<fdr_cutoff) %>% pull(gene_id) %>% unique %>% length()
n_tx <- circadian_analysis %>% group_by(gene_id) %>% mutate(n_tx = n()) %>% filter(n_tx>1 & transcript<fdr_cutoff) %>% pull(transcript_id) %>% length()

heatmap_tx <- ggplot(df_1) + 
  geom_tile(aes(y=transcript_id, x=time, fill=tx_exp), width=1) +
  facet_wrap(~rep, nrow = 1) + 
  scale_fill_gradient2_tableau(name="Exp/Med", limits = c(-2.5,2.5), palette = "Green-Blue Diverging")  +
  new_scale_fill() +
  geom_rect(aes(xmin=xmin, xmax=xmax, ymin=n_transcripts_rhy, ymax=n_transcripts_rhy + 35, fill = col), data = LD_bars,inherit.aes=FALSE) +
  scale_fill_manual(values = c("black","grey80"), guide="none") +
  ylab(paste(n_tx, "cycling tx","from", n_genes,  sep= " ")) + 
  xlab("circadian time") + theme_minimal(base_size = base_size) +
  scale_x_discrete(expand=expansion()) + scale_y_discrete(expand=expansion()) +
  theme_legend +
  theme(axis.text.y = element_blank(), axis.line = element_blank(),
        axis.ticks.y = element_blank(), panel.spacing.x = unit(0.2, "lines"), 
        legend.position = "right", panel.grid = element_blank(), panel.grid.major = element_blank(), legend.text = element_text(size=6), axis.text.x=element_blank(), axis.title.x=element_blank())

heatmap_gene <- ggplot(df_1 %>% group_by(gene_id, time, rep) %>% filter(row_number()==1)) + 
  geom_tile(aes(y=transcript_id, x=time, fill=gene_exp), width=1) +
  facet_wrap(~rep, nrow=1) + 
  scale_fill_gradient2_tableau(name="Exp/Med", limits = c(-2.5,2.5), palette = "Green-Blue Diverging")  +
  scale_x_discrete(expand=expansion()) + scale_y_discrete(expand=expansion()) +
  new_scale_fill() +
  geom_rect(aes(xmin=xmin, xmax=xmax, ymin=n_genes, ymax=n_genes + 35, fill = col), data = LD_bars) +
  scale_fill_manual(values = c("black","grey80")) +
  ylab(sprintf("%d genes with AS and\n >= 1 rhythmic isoform", n_genes)) + 
  xlab("circadian time") + theme_minimal(base_size = base_size) +
  theme(axis.text.y = element_blank(), axis.line = element_blank(), aspect.ratio = 
          4*n_genes/n_transcripts_rhy,
        axis.ticks.y = element_blank(), panel.spacing.x = unit(0.2, "lines"), 
        legend.position = "none", panel.grid = element_blank(), panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "black", colour="black", color="black"), legend.title=element_blank())


##-----------Differential transcript circadian usage analysis----------------##
cts <- data.frame(gene_id = rowData(circadian_dds)$gene_id,
                  feature_id = rowData(circadian_dds)$transcript_id,
                  counts(circadian_dds), check.names = FALSE)
d <- DRIMSeq::dmDSdata(counts = cts  %>% filter(!is.na(feature_id)), samples = exp_design)
n <- 12
n.small <- 4
d <- DRIMSeq::dmFilter(d,
                       min_samps_feature_expr=n.small, min_feature_expr=10,
                       min_samps_feature_prop=n.small, min_feature_prop=0.1,
                       min_samps_gene_expr=n-n.small, min_gene_expr=10)

dxd <- DEXSeqDataSet(countData = as.matrix(counts(d)[,-c(1,2)]),
                     sampleData = DRIMSeq::samples(d),
                     design =~sample_id + exon + batch:exon + inphase*exon + outphase*exon,
                     featureID = counts(d)$feature_id,
                     groupID = counts(d)$gene_id %>% sub(':','loc',.))
BPPARAM <-  BiocParallel::MulticoreParam(2)
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, quiet=TRUE, fitType = 'local', BPPARAM=BPPARAM)
dxd <- testForDEU(dxd, reducedModel=~sample_id + exon + batch:exon + inphase + outphase, BPPARAM=BPPARAM)
annotation$gene_id<-  sub(':','loc',annotation$gene_id) #why not above?
tx2gene<- annotation[c("transcript_id", "gene_id")]

dxr <- DEXSeqResults(dxd, independentFiltering=TRUE)
qval <- DEXSeq::perGeneQValue(dxr)

pScreen <- qval
pConfirmation <- matrix(dxr$pvalue, ncol = 1)
dimnames(pConfirmation) <- list(dxr$featureID, "transcript")

stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=TRUE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.1, allowNA=TRUE)
suppressWarnings({
  dtu_padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                 onlySignificantGenes=TRUE)
})
dtu_padj %<>%
  merge(tx2gene, by.x = "txID", by.y = "transcript_id") %>%
  arrange(gene) 

write_rds(dtu_padj, "../results/dtu_padj.rds", compress = "bz2")

dtu_padj %>% dplyr::select(txID, gene, geneID, transcript) %>% 
  filter(gene<fdr_cutoff) %>%
  write_excel_csv(file = "../lists/DTU_genes.csv")

clock_genes <- c("Arntl", "Per1", "Dbp", "Rorc", "Clock", "Hlf", "Bhlhe40", "Nr1d2", "Cry1", "Cry2", "Nr1d1", "Per2", "Bhlhe40", "Arntl2", "Nampt", "Ciart")
dtu_genes<- dtu_padj %>% dplyr::select(txID, gene, geneID, transcript) %>% filter(gene<fdr_cutoff) %>% pull(geneID) %>%unique()

##--------Examples of genes with multiple isoforms in-phase and out-of-phase--------##
colnames(anno_bed) <- c("chrom", "chromStart", "chromEnd", "transcript_id", "score", "strand", "thickStart", "thickEnd", "itemRgb", "exons_nb", "exons_size", "exons_start",  "productivity","gene_name")

df_3 <- abundance_no_batch_effect %>% 
  as.data.frame %>% rownames_to_column(var = "transcript_id") %>%
  filter(transcript_id %in% rownames(circadian_dds)) %>%
  pivot_longer(-transcript_id, names_to = c("sample","rep"), names_sep="_") %>%
  separate(sample, c("cond","time"), convert = TRUE, sep = 2) %>%
  merge(tx2gene, by = "transcript_id") %>% merge(anno_bed, by="transcript_id") %>% inner_join(.,annotation, by= "transcript_id") %>% mutate(productivity = factor(ifelse(grepl("PRO", productivity), "productive", "non-productive"), levels = c("productive", "non-productive"))) %>% 
  group_by(mgi_symbol) %>%
  mutate(tx_seq = as.numeric(factor(transcript_id)))
##------------ WT timecourse-------------------##
clock_tx <- rowData(circadian_dds) %>% as.data.frame %>%
  filter(mgi_symbol %in% clock_genes) %$%
  transcript_id
time_course_plots_WT <- ggplot(data = filter(df_3, transcript_id %in% clock_tx)) + 
  geom_point(aes(x=time, y=value, color = productivity), size=0.4) + 
  geom_smooth(aes(x=time, y=value, group= transcript_id, color = productivity), size = 0.7, formula = "y ~ x", method = "loess", se=FALSE) +
  facet_wrap(~mgi_symbol, scales = "free_y", nrow=3) + 
  theme(legend.position = "right") + scale_x_continuous(breaks = seq(0, 20, 4)) + 
  xlab("circadian time") +   scale_y_continuous(breaks = scales::extended_breaks(n=4)) + 
  ylab(bquote("log"[2]~"isoform expression")) + theme_classic(base_size = base_size) + 
  theme(panel.grid.major = element_blank(), 
        strip.text = element_text(size=8), panel.grid.minor = element_blank(), 
        legend.position = "top", legend.title=element_blank(), 
        aspect.ratio = 0.8, legend.text = element_text(size=8), strip.background = element_rect(fill="white",color="white")) + 
  scale_color_colorblind() #+ ggtitle( "Expression timecourse for clock genes")

ggsave(paste(fig_dir, "/2_WT_timecourse_circgenes.svg", sep=""), plot = time_course_plots_WT, height=5, width=5)

##------------Figure Splicing DTU WT timecourse-------------------##
dtu_tx <- dtu_padj %>% merge(annotation,by.x="txID", by.y="transcript_id") %>% 
  filter(geneID %in% dtu_genes)  %>% filter(mgi_symbol != "Gm8883")  %$%
  txID

#to order the facet, make the order here and use it as level later
#neworder<- c("Prkd3","Ugt2b38","Pdcd4","Insig2","Mup21","Clpx","Pxmp4","Slc1a2","Ccdc127","Rpp21","Pla2g12a","Fgg","Atp5g1","Tmem256","0610031O16Rik")
neworder<- c("Prkd3","Ugt2b38","Pdcd4","Mup21","Pxmp4","Slc1a2","Ccdc127","Rpp21","Pla2g12a","Fgg","Atp5g1","Tmem256", "Insig2", "Mocs2", "Slc11a2") #to make it smaller

time_course_plots_DTU <- ggplot(data = filter(df_3, transcript_id %in% dtu_tx & mgi_symbol %in% neworder) %>% mutate(mgi_symbol=factor(mgi_symbol,levels=neworder)) ) +
  geom_point(aes(x=time, y=value, color = productivity), size=0.4) +
  geom_smooth(aes(x=time, y=value, group= transcript_id, color = productivity), size = 0.7, formula = "y ~ x", method = "loess", se=FALSE) +
  facet_wrap(~mgi_symbol, scales = "free_y",  nrow=3 )+
  scale_y_continuous(breaks = scales::extended_breaks(n=4)) + 
  scale_x_continuous(breaks = seq(0, 20, 4)) + xlab("circadian time") +
  ylab(bquote("log"[2]~"isoform expression")) + theme_classic(base_size = base_size) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "none", strip.text = element_text(size=11), strip.background = element_rect(fill="white",color="white"),
        legend.title=element_blank(),plot.background = element_rect(fill= "white"), aspect.ratio = 0.85) + 
  scale_color_colorblind() #+ ggtitle("Differential transcript circadian usage")
ggsave(paste(fig_dir, "/2_DTU_timecourse.svg", sep=""), plot = time_course_plots_DTU, height=3.5, width=5)

time_course_plots_DTU <- ggplot(data = filter(df_3, transcript_id %in% dtu_tx ) %>% mutate(mgi_symbol=mgi_symbol) ) +
  geom_point(aes(x=time, y=value, color = productivity), size=0.4) +
  geom_smooth(aes(x=time, y=value, group= transcript_id, color = productivity), size = 0.7, formula = "y ~ x", method = "loess", se=FALSE) +
  facet_wrap(~mgi_symbol, scales = "free_y",  nrow=3 )+
  scale_y_continuous(breaks = scales::extended_breaks(n=4)) + 
  scale_x_continuous(breaks = seq(0, 20, 4)) + xlab("circadian time") +
  ylab(bquote("log"[2]~"isoform expression")) + theme_classic(base_size = base_size) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "none", strip.text = element_text(size=11), strip.background = element_rect(fill="white",color="white"),
        legend.title=element_blank(),plot.background = element_rect(fill= "white"), aspect.ratio = 0.85) + 
  scale_color_colorblind() #+ ggtitle("Differential transcript circadian usage")

ggsave(paste(fig_dir, "/2_DTU_timecourse_all.svg", sep=""), plot = time_course_plots_DTU, height=10, width=10)
                      
manualcolors<-c("#80ce9f", '#ba56c5', '#c9503b', '#404f51', '#c8a450', '#8dcd49', '#c75985', '#a4b0ba', '#757ac5', '#537137', '#6e402e', '#5a355f', "black")

polar_clock_tx <- filter(circadian_analysis, transcript_id %in% clock_tx) %>%
  left_join(., mart[,2:3], by='gene_id') %>%
  mutate(productivity = ifelse(productivity == "PRO", "PRO", "ncRNA")) %>%
  ggplot(aes(y=phase, x=amp, col=mgi_symbol, shape = productivity))  +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 12, fill = '#FAF8F8', color="white")+
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 12, ymax = 24, fill = '#D6D6D6', alpha = 0.05, color="white")+
  geom_vline(aes(xintercept = x), data.frame(x = seq(4)*1), color = "#928E8E", linewidth=0.6) +
  geom_vline(aes(xintercept = 4), linewidth=0.6) +
  geom_hline(aes(yintercept = y), data.frame(y = seq(0,24,2)), color="#928E8E", linewidth=0.6)  +
  geom_point(size=2) + 
  scale_shape_manual(values=c(1, 16)) +
  coord_polar(theta="y", start=0) + 
  scale_y_continuous(breaks = seq(0,24,2), limits = c(0,24), expand = expansion()) + 
  theme_void(base_size = base_size*2) +
  theme(legend.title = element_blank(), axis.ticks = element_blank(), text = element_text(size = base_size+2), axis.text.y = element_blank(), axis.text.x = element_text(face="plain", size=rel(1)), strip.text  = element_blank(), legend.position = "right", aspect.ratio = 1, legend.text=element_text(size=12), legend.key.height = unit(0.3, 'cm')) +
  scale_color_manual(values = manualcolors) +
  guides(color = guide_legend(override.aes = list(size = 2, linetype = 0, fill=NA))) +
  labs(color="Gene") + 
  scale_x_continuous(expand = expansion(), limits = c(0,4.5)) 

############################################################################################################################################
############################################################################################################################################
data <- read.table("../data/data_for_final_plot.tsv", sep = "", header = TRUE) %>% 
  filter(tecov > 0.8)
te_families <- read.table("../data/mm39_families.tsv") %>% setNames(c("te_name", "family", "count", "mean size"))#


data$famSub = te_families[match(data$te_name,
                                te_families$te_name), "family"]
tmp = str_split(data$famSub, "/", simplify = TRUE)
data$fam = tmp[, 1]
data$subfam = tmp[, 2]
data[is.na(data$famSub), "famSub"] = "merged"
data[is.na(data$subfam), "subfam"] = "merged"
data[is.na(data$fam), "fam"] = "merged"
data[data$subfam == "", "subfam"] = data[data$subfam == "", "famSub"]

data <- dplyr::select(data, c(matching_transcript, fam, subfam)) %>%
  setNames(c('transcript_id', 'fam', 'subfam'))

plot_df <- dplyr::select(circadian_analysis, c(gene_id, transcript_id, transcript, phase)) %>%
  filter(transcript < fdr_cutoff) %>%
  left_join(., data, by='transcript_id')

polar_TE_fun <- function(subfam_name){
  col_fill <- ifelse(filter(plot_df, subfam == subfam_name) %>% .[1,5] == "SINE", "orange",
                     ifelse(filter(plot_df, subfam == subfam_name) %>% .[1,5] == "LINE", "firebrick4",
                            ifelse(filter(plot_df, subfam == subfam_name) %>% .[1,5] == "LTR", "darkblue",
                                   ifelse(filter(plot_df, subfam == subfam_name) %>% .[1,5] == "DNA", "forestgreen","grey"))))
  max_val <- filter(plot_df, subfam == subfam_name) %>% mutate(phase = round(phase)) %>% group_by(phase) %>% summarise(n = n()) %>% .$n %>% max()
  ggplot(filter(plot_df, subfam == subfam_name) %>% mutate(phase = round(phase)) ) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 12, fill = '#FAF8F8')+
    geom_rect(xmin = -Inf, xmax = Inf, ymin = 12, ymax = 24, fill = '#D6D6D6')+
    geom_vline(aes(xintercept = x), data.frame(x = seq(3)*(max_val/3)), color = "#928E8E", linewidth=0.3) +
    geom_vline(aes(xintercept = max_val+(max_val/5))) +
    geom_hline(aes(yintercept = y), data.frame(y = seq(0,24,2)), color="#928E8E", linewidth=0.3) +
    geom_histogram(aes(y=phase), binwidth = 1, center=0.5, fill=rep(col_fill, 24), color="white", position="dodge") +
     coord_polar(theta="y", start=0) +  
    scale_y_continuous(breaks = seq(0,24,2), limits = c(0,24), expand = expansion()) + 
    theme_void(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, size=base_size-1), plot.caption = element_text(hjust = 0.5),axis.text.x = element_text(face="plain", size=rel(0.7)), strip.text  = element_blank(), legend.position = "none", aspect.ratio = 1) + 
    scale_x_continuous(expand = expansion(), limits = c(0,max_val+(max_val/3))) + 
    labs(title = subfam_name, caption = paste("n=", nrow(filter(plot_df, subfam == subfam_name)), sep=""))
}

polar_TE <- plot_grid(polar_TE_fun("Alu"), polar_TE_fun("B2"), polar_TE_fun("B4"), polar_TE_fun("ID"), polar_TE_fun("ERVL-MaLR"), polar_TE_fun("ERVL"), polar_TE_fun("ERVK"), polar_TE_fun("ERV1"), polar_TE_fun("L2"), polar_TE_fun("CR1"), polar_TE_fun("L1"), polar_TE_fun("RTE-BovB"), polar_TE_fun("hAT-Charlie"), polar_TE_fun("hAT-Ac"), nrow=4)
polar_TE <- plot_grid(polar_TE_fun("Alu"), polar_TE_fun("B2"), polar_TE_fun("B4"), polar_TE_fun("ID"), polar_TE_fun("ERVL-MaLR"), polar_TE_fun("ERVL"), polar_TE_fun("ERVK"), polar_TE_fun("ERV1"), polar_TE_fun("L2"), polar_TE_fun("CR1"), polar_TE_fun("L1"), polar_TE_fun("RTE-BovB"), nrow=3)

ggsave(paste(fig_dir, "/2_heatmap_tx.svg", sep=""), plot=heatmap_tx, width=2.5, height = 4)
ggsave(paste(fig_dir, "/2_phasedistribution1.svg", sep=""), plot=phase_distribution1, width=2.5, height = 2.5)
ggsave(paste(fig_dir, "/2_phasedistribution2.svg", sep=""), plot=phase_distribution2, width=2.5, height = 2.5)

ggsave(paste(fig_dir, "/2_polar_clock_tx.svg", sep=""), plot=polar_clock_tx, width=5.5, height = 5)
ggsave(paste(fig_dir, "/2_timecourseplots_DTU.svg", sep=""), plot=time_course_plots_DTU, width=5, height = 4)
ggsave(paste(fig_dir, "/2_polarTE.svg", sep=""), plot=polar_TE, width=4, height = 4)
ggsave(paste(fig_dir, "/2_mosaic_annotation.svg", sep=""), plot=mosaic_annotation, width=3, height = 3)
ggsave(paste(fig_dir, "/2_mosaic_product.svg", sep=""), plot=mosaic_prod, width=3, height = 3)

##############################################
# old/extra part not in the figures anymore###
##############################################

#neworder<- c("Prkd3","Ugt2b38","Pdcd4","Insig2","Pxmp4","Clpx", "Mup21", "Slc1a2","Ccdc127","Rpp21","Pla2g12a","0610031O16Rik","Atp5g1","Tmem256","Fgg")
#time_course_plots_DTU <- ggplot(data = filter(df_3, transcript_id %in% dtu_tx) %>% mutate(mgi_symbol=factor(mgi_symbol,levels=neworder)) ) +
#  geom_point(aes(x=time, y=value, color = productivity), size=0.4) +
#  geom_smooth(aes(x=time, y=value, group= transcript_id, color = productivity), size = 0.7, formula = "y ~ x", method = "loess", se=FALSE) +
#  facet_wrap(~mgi_symbol, scales = "free_y",  nrow=3 )+
#  scale_y_continuous(breaks = scales::extended_breaks(n=4)) + 
#  scale_x_continuous(breaks = seq(0, 20, 4)) + xlab("circadian time") +
#  ylab(bquote("log"[2]~"isoform expression")) + theme_classic(base_size = base_size) + 
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#        legend.position = "none", strip.text = element_text(size=8), strip.background = element_rect(fill="white",color="white"),
#        legend.title=element_blank(), aspect.ratio = 0.85) + 
#  scale_color_colorblind()

### TE part, file + code from fabien, K made some changes only for aesthetics and layout. 
# load("../data/less_good_than_pringles_but_still.Rdata")
# theme_polar = theme(panel.border = element_blank(), plot.background = element_blank(), panel.background = element_blank(), plot.title = element_text(hjust = 0.5))
# colors = c("#BE55FF", "#2FCC56", "#CC2F9D", "#AD8165")
# names(colors) = keep_fam
###

#list_grid2 = lapply(keep_fam, function(te_fam) {
#  subfdf = plot_df[plot_df$fam == te_fam, ]
#   subfdf$subfam = as.factor(subfdf$subfam)
#   list_subgrid = lapply(unique(subfdf$subfam), function(subfam) {
#     max_y = max(aggregate(subfdf[subfdf$fam == te_fam & subfdf$subfam == subfam, "x"], by = list(subfdf[subfdf$fam == te_fam & subfdf$subfam == subfam, c("phase")]), FUN = sum)$x)
#     total_transcripts = sum(subfdf[subfdf$fam == te_fam & subfdf$subfam == subfam, "x"])
#     tmp_plot = ggplot(data = subfdf[subfdf$fam == te_fam & subfdf$subfam == subfam, ], aes(x = phase, y = x, fill = fam)) + geom_col(width = 0.5, position = "stack") + ylab(NULL) + xlab(NULL) + ggtitle(paste(subfam, " (", total_transcripts, ")", sep = "")) + geom_hline(yintercept = max_y) + scale_fill_manual(values = colors) + scale_x_continuous(breaks = c(0, 6, 12, 18)) + scale_y_continuous(expand = expansion(mult = c(0, .1))) + theme_polar + theme(axis.ticks = element_blank(), legend.position = "None", text = element_text(size = 7), axis.text.y = element_blank()) + coord_polar()
#   })   + ggthemes::scale_fill_tableau(breaks = keep_fam, name = "", drop = F)
#        print(te_fam)
#        print(list_subgrid)
#   if (te_fam == "DNA") {
#     subsubgrid = plot_grid(plotlist = list_subgrid, ncol = 1) + theme(plot.background = element_rect(color = "black", linetype = "solid", linewidth = 0.5))
#     max_y = max(aggregate(leftover[leftover$fam %in% c("Low_complexity", "Simple_repeat"), "x"], by = list(leftover[leftover$fam %in% c("Low_complexity", "Simple_repeat"), c("phase")]), FUN = sum)$x)
#     total_transcripts = sum(leftover[leftover$fam %in% c("Low_complexity", "Simple_repeat"), "x"])
#     added_plot = ggplot(data = leftover[leftover$fam %in% c("Low_complexity", "Simple_repeat"), ], aes(x = phase, y = x)) + geom_col(width = 0.5, fill = "black", position = "stack") + ggtitle(paste("Simple repeats", " (", total_transcripts, ")", sep = "")) + geom_hline(yintercept = max_y)  + theme_polar + theme(axis.ticks = element_blank(),legend.position = "None", text = element_text(size = 7), axis.text.y = element_blank()) + coord_polar()#
#    return(plot_grid(subsubgrid, added_plot, rel_heights = c(3, 1), ncol = 1) + draw_text("DNA", hjust = 2.4, vjust = 16, size= 7, fontface= "bold")+ panel_border(color = "black", size = 0.5, linetype = 1, remove = FALSE))
#   } else {
#     return(plot_grid(plotlist = list_subgrid, ncol = 1) + theme(plot.background = element_rect(color = "black", linetype = "solid", linewidth = 0.5)))
#   }
#        plot_grid(plotlist =  list_subgrid, nrow = 1)#
# })
# grid2 = plot_grid(plotlist = list_grid2, align = c("hv"), ncol = 2, nrow = 2, labels = keep_fam[1:3], label_size= 8, label_x = 0, label_y = 0,hjust = -0.5, vjust = -0.5) #-25 for middle ish
# final_grid<- grid2
# layout<- c(area(1,1,3,2), area(1,3,3,4), area(1,5,3,5),  area(4,1,6,5), area(7,1,9,5), area(1,6,9,9))

#list_grid2 = lapply(keep_fam, function(te_fam){
#  subfdf = plot_df[plot_df$fam == te_fam,]
#  subfdf$subfam = as.factor(subfdf$subfam)
#  list_subgrid = lapply(unique(subfdf$subfam), function(subfam) {
#    max_y = max(aggregate(subfdf[subfdf$fam == te_fam & subfdf$subfam == subfam, "x"], by = list(subfdf[subfdf$fam == te_fam & subfdf$subfam == subfam, c("phase")]), FUN = sum)$x)
#    total_transcripts = sum(subfdf[subfdf$fam == te_fam & subfdf$subfam == subfam, "x"])
#    tmp_plot = ggplot(data = subfdf[subfdf$fam == te_fam & subfdf$subfam == subfam,], aes(x = phase, y = x, fill = subfam)) + geom_col(width = 0.5, position = "stack") + theme(legend.position = "None", text = element_text(size = 8),axis.ticks = element_blank(), plot.title = element_text(size=8)) + ylab(NULL) + xlab(NULL) + ggtitle(paste(subfam, " (", total_transcripts, ")", sep = "")) + ggthemes::scale_fill_colorblind(breaks = unique(subfdf$subfam), name = "", drop = F) + theme_polar + geom_hline(yintercept = max_y)  + coord_polar() + theme(axis.text.y = element_blank())
#     		tmp_plot = ggdraw(tmp_plot + theme(axis.text.y = element_blank())) + draw_grob(plot_y, x = 0.725, y = 1.05, hjust = 1, vjust = 1, scale = 0.6)
#     		tmp_plot = ggdraw(tmp_plot + theme(axis.text.y = element_blank())) + draw_text(max_y, x = 0.5, y = 0.5)
##  })
#legend1 <- get_legend(heatmap_gene + guides(color = guide_legend(nrow = 2)) +theme(legend.position = "bottom"))
#legend2 <- get_legend(mosaic_annotation +guides(color = guide_legend(nrow = 2))+theme(legend.position = "bottom"))
#legend3 <- get_legend(mosaic_prod +guides(color = guide_legend(nrow = 2)) + theme(legend.position = "bottom"))
#plot_legends<- plot_grid(legend2, legend3, nrow=1)
#fig2<-plot_grid(plot_legends,f2,ncol=1,rel_heights = c(0.05, 1))#

#----------Plotting Isoforms Structures-------------------##
####
#gtf <- rtracklayer::import("../results/FLAIR_2023-02-22/transcriptome_ext.gtf") %>% sortSeqlevels %>% sort %>% dplyr::as_tibble() %>% filter(transcript_id %in%  anno_bed$transcript_id) %>% left_join(.,(anno_bed %>% dplyr::select(c("transcript_id", "productivity"))), by="transcript_id")
#gene_structure <- function(gene_name,...){
#  gene_annotation_from_gtf <- gtf %>% 
#    filter(gene_id == !!gene_name) %>%  
#    dplyr::select(seqnames, start, end, strand, type, gene_id, transcript_id, source,productivity)
#  gene_rescaled <- shorten_gaps(
#    exons = filter(gene_annotation_from_gtf, type == "exon"), 
#    introns = to_intron(filter(gene_annotation_from_gtf, type == "exon"), "transcript_id"), 
#    group_var = "transcript_id"
#  )
#  gene_rescaled_exons <- gene_rescaled %>% dplyr::filter(type == "exon") 
#  gene_rescaled_introns <- gene_rescaled %>% dplyr::filter(type == "intron") 
#  ggplot(gene_rescaled_exons, aes(xstart = start, xend = end, y = transcript_id)) +
#    geom_range(aes(fill=(productivity!="PRO")), height=0.5) +
#    geom_intron(data = gene_rescaled_introns, aes(strand=strand), 
#                arrow=grid::arrow(ends = "last", length = grid::unit(0.075, "inches")),
#                color="gray20",
#                arrow.min.intron.length=2) +
#    theme_void(base_size = base_size) +
#   theme(legend.position = "none", aspect.ratio = 0.2,...) +
#    scale_fill_colorblind() + scale_color_brewer()
#}


# Full TE dynamics

#list_full_TE_circ <- c("16loc84745000", "4loc78513000", "4loc49376000", "19loc3084000", "6loc18847000", "11loc86811000", "ENSMUSG00000099291", "5loc16729000", "9loc56223000")

#labels <- c("ERV3","ERV2","L1","ERV1","ERV1","ERV1","snRNA","ERV2","rRNA")
#names(labels) <- list_full_TE_circ#
#
#full_TE_timecourse <- ggplot(data = filter(df_3, gene_id %in% list_full_TE_circ)) + 
# geom_point(aes(x=time, y=value, color = "black"), size=0.3) + 
#  geom_smooth(aes(x=time, y=value, group= transcript_id, color = productivity), size = 0.5, formula = "y ~ x", method = "loess", se=FALSE) + 
#  scale_x_continuous(breaks = seq(0, 20, 4)) + xlab("circadian time") + 
#  scale_y_continuous(breaks = scales::breaks_extended(n=4)) + 
#  ylab(bquote("log"[2]~"isoform expression")) + theme_classic(base_size = base_size) + 
#  theme(panel.grid.major = element_blank(), strip.background = element_blank(),
#        panel.grid.minor = element_blank(), legend.position = "none", 
#        legend.title=element_blank(), aspect.ratio = 0.75) + 
#  scale_color_manual(values=c("black", "black")) + 
#  facet_rep_wrap(~ gene_id, scales = "free_y", nrow = 2, labeller=labeller(gene_id = labels))



#rescaled exons/intron for plotting isoform structures#
#gene_of_interest <- c("ENSMUSG00000000876")
#transcript_list <- c("ENSMUST00000000896", "ENSMUST00000109703")
#gene_annotation_from_gtf <- gtf %>% filter(!is.na(gene_id), gene_id == gene_of_interest) %>%  dplyr::select(seqnames, start, end, strand, type, gene_id, transcript_id)
#gene_exons <- filter(gene_annotation_from_gtf, type == "exon")
#genes_rescaled <- shorten_gaps(exons = gene_exons, introns = to_intron(gene_exons, "transcript_id"), group_var = "transcript_id")
#rescaled_exons <- genes_rescaled %>% dplyr::filter(type == "exon")
#rescaled_introns <- genes_rescaled %>% dplyr::filter(type == "intron")
#Pxmp4_rescaled_structure <-  ggplot(filter(rescaled_exons, transcript_id %in% transcript_list), aes(xstart = start, xend = end, y = transcript_id, fill = transcript_id)) + geom_range(height = 0.25) + geom_intron(data = filter(rescaled_introns, transcript_id %in% transcript_list), aes(strand = strand), arrow.min.intron.length = 500) + scale_fill_manual(values = c("#d0ecb6", "#d5e6ff")) + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", text=element_blank())
#svg("../Manuscript_Draft/Fig2/Pxmp4.svg", width=3, height=3)
#Pxmp4_rescaled_structure
#dev.off()
#to inverse the 2 isforms : Csnk1d_rescaled_structure + scale_y_discrete(limits = rev)
# jtk_analysis <- meta2d("NA", outdir = "../results/metaout/", "csv", colData(dds)$time, cycMethod = "JTK",
#                        analysisStrategy = "auto", outputFile = FALSE, minper = 24, maxper = 24,
#                        inDF = abundance_no_batch_effect %>% as.data.frame %>% rownames_to_column(var="txID"))
# 
# res_jtk <- jtk_analysis$meta %>% 
#               dplyr::rename(pvalue = JTK_pvalue, padj = JTK_BH.Q, transcript_id = CycID) %>% 
#               left_join(tx2gene) %>% dplyr::rename(groupID=gene_id)
# 
# pScreen <- perGeneQValueHack(res_jtk)
# pConfirmation <- matrix(res_jtk$pvalue, ncol = 1)
# 
# dimnames(pConfirmation) <- list(res_jtk$transcript_id, "transcript")
# 
# stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
#                       pScreenAdjusted=TRUE, tx2gene=tx2gene[c("transcript_id","gene_id")])
# 
# stageRObj <- stageWiseAdjustment(stageRObj, method="dte", alpha=fdr_cutoff, allowNA=TRUE)
# 
# suppressWarnings({
#   circadian_analysis <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=FALSE)
# })


##------------ Lies-------------------##

# manualcolors<-c('magenta','forestgreen', 'red2', 'orange', 'cornflowerblue', 
#                        'black', 'darkolivegreen4', 'indianred1', 'tan4', 'darkblue', 
#                        'mediumorchid1','firebrick4',  'yellowgreen', 'lightsalmon', 'tan3',
#                        "tan1",'darkgray', 'wheat4', '#DDAD4B', 'chartreuse', 
#                        'seagreen1', 'moccasin', 'mediumvioletred', 'seagreen','cadetblue1',
#                        "darkolivegreen1" ,"tan2" ,   "tomato3" , "#7CE3D8","gainsboro")

############################################################################################################################################
############################################################################################################################################
# f2_LC <- plot_grid(plot_grid(heatmap_tx, phase_distribution, plot_grid(mosaic_annotation, mosaic_prod, ncol=1), polar_clock_tx, nrow=1, rel_widths = c(1.5,1,0.9,2)), plot_grid(time_course_plots_DTU, plot_grid(NULL, plot_grid(Pxmp4_structure, NULL, rel_widths = c(6,1)), NULL, ncol=1, rel_heights = c(1,1.2,3)), nrow=1, rel_widths = c(2,1)), polar_TE, ncol=1, rel_heights = c(1.3,1.7,0.8))
# ggsave("../results/figure2bis.pdf", plot=f2_LC, width=7, height = 6)

############################################################################################################################################
############################################################################################################################################

#gtf <- rtracklayer::import(gtf_path) %>% sortSeqlevels %>% sort %>% dplyr::as_tibble() %>% filter(transcript_id %in%  anno_bed$transcript_id) %>% left_join(.,(anno_bed %>% dplyr::select(c("transcript_id", "productivity"))), by="transcript_id")

# Pxmp4 : ENSMUST00000000896 & ENSMUST00000109703 - ENSMUSG00000000876
#gene_of_interest <- c("ENSMUSG00000000876")
#transcript_list <- c("ENSMUST00000000896", "ENSMUST00000109703")
#gene_annotation_from_gtf <- gtf %>% filter(!is.na(gene_id), gene_id == gene_of_interest & transcript_id %in% transcript_list) %>%  dplyr::select(seqnames, start, end, strand, type, gene_id, transcript_id)
#gene_exons <- filter(gene_annotation_from_gtf, type == "exon" & transcript_id %in% transcript_list)

#Pxmp4_structure <-  ggplot(filter(gene_exons, transcript_id %in% transcript_list), aes(xstart = start, xend = end, y = transcript_id, fill = transcript_id)) + 
#  geom_rangeX(data = gene_exons, height = 0.4, fill="grey90") + 
#  geom_rangeX(data = filter(gene_annotation_from_gtf, type == "CDS"), height = 0.4, fill="black") +
#  geom_intron(data = to_intron(gene_exons, "transcript_id"), aes(strand = strand), arrow.min.intron.length = 1400) + 
#  theme_minimal() + 
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", axis.text=element_blank(), axis.title = element_blank()) + 
#  scale_y_discrete(limits = rev) + 
#  ggtitle("Pxmp4") 

# also put Fgg structure and Rpp21

# Fgg : ENSMUST00000048486 & ENSMUST00000194175 - ENSMUSG00000033860
#gene_of_interest <- c("ENSMUSG00000033860")
#transcript_list <- c("ENSMUST00000048486", "ENSMUST00000194175")
#gene_annotation_from_gtf <- gtf %>% filter(!is.na(gene_id), gene_id == gene_of_interest & transcript_id %in% transcript_list) %>%  dplyr::select(seqnames, start, end, strand, type, gene_id, transcript_id)
#gene_exons <- filter(gene_annotation_from_gtf, type == "exon" & transcript_id %in% transcript_list)

#Fgg_structure <-  ggplot(filter(gene_exons, transcript_id %in% transcript_list), aes(xstart = start, xend = end, y = transcript_id, fill = transcript_id)) + 
#  geom_rangeX(data = gene_exons, height = 0.4, fill="grey90") + 
#  geom_rangeX(data = filter(gene_annotation_from_gtf, type == "CDS"), height = 0.4, fill="black") +
#  geom_intron(data = to_intron(gene_exons, "transcript_id"), aes(strand = strand), arrow.min.intron.length = 1000) + 
#  theme_minimal() + 
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", axis.text=element_blank(), axis.title = element_blank()) + 
#  scale_y_discrete(limits = rev) + 
#  ggtitle("Fgg")
#Fgg_structure

# Rpp21 : ENSMUST00000173486 & ENSMUST00000025319 - ENSMUSG00000024446
#gene_of_interest <- c("ENSMUSG00000024446")
#transcript_list <- c("ENSMUST00000173486", "ENSMUST00000025319")
#gene_annotation_from_gtf <- gtf %>% filter(!is.na(gene_id), gene_id == gene_of_interest & transcript_id %in% transcript_list) %>%  dplyr::select(seqnames, start, end, strand, type, gene_id, transcript_id)
#gene_exons <- filter(gene_annotation_from_gtf, type == "exon" & transcript_id %in% transcript_list)

#Rpp21_structure <-  ggplot(filter(gene_exons, transcript_id %in% transcript_list), aes(xstart = start, xend = end, y = transcript_id, fill = transcript_id)) + 
#  geom_rangeX(data = gene_exons, height = 0.4, fill="grey90") + 
#  geom_rangeX(data = filter(gene_annotation_from_gtf, type == "CDS"), height = 0.4, fill="black") +
#  geom_intron(data = to_intron(gene_exons, "transcript_id"), aes(strand = strand), arrow.min.intron.length = 1000) + 
#  theme_minimal() + 
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", axis.text=element_blank(), axis.title = element_blank()) + 
#  scale_y_discrete(limits = rev) + 
#  ggtitle("Rpp21")
#Rpp21_structure 



##------------ Lies-------------------##

#df_1 <- abundance_no_batch_effect %>% 
#  as.data.frame %>% rownames_to_column(var = "transcript_id")  %>%
#  filter(transcript_id %in% circadian_transcripts$transcript_id) %>%
#  pivot_longer(-transcript_id, names_to = c("sample","rep"), names_sep="_", values_to = "tx_exp") %>%
#  separate(sample, c("cond","time"), convert = TRUE, sep = 2) %>%
#  group_by(transcript_id) %>%
#  dplyr::select(-cond) %>%
#  left_join(tx2gene) %>%
#  group_by(gene_id, time, rep) %>%
#  mutate(gene_exp = mean(tx_exp)) %>%
#  group_by(gene_id) %>%
#  mutate(gene_exp = scale(gene_exp)) %>%
#  group_by(transcript_id) %>%
#  mutate(tx_exp = scale(tx_exp)) %>%
#  mutate(transcript_id = factor(transcript_id, levels = phase_sorted_tx, ordered = TRUE),
#         time = factor(time, levels = seq(0, 20, 4), ordered = TRUE)) %>%
#  arrange(transcript_id)
  
#heatmap_tx2 <- ggplot(df_1) + 
#  geom_tile(aes(y=transcript_id, x=time, fill=tx_exp), width=1) +
#  facet_wrap(~rep, nrow = 1) + 
#  scale_fill_gradient2_tableau(name="Exp/Med", limits = c(-2.5,2.5), palette = "Classic Red-White-Green")  +
#  new_scale_fill() +
#  geom_rect(aes(xmin=xmin, xmax=xmax, ymin=n_transcripts_rhy, ymax=n_transcripts_rhy + 35, fill = col), data = LD_bars,inherit.aes=FALSE) +
#  scale_fill_manual(values = c("black","grey80"), guide="none") +
#  ylab("") + 
#  xlab("circadian time") + theme_minimal(base_size = base_size) +
#  scale_x_discrete(expand=expansion()) + scale_y_discrete(expand=expansion()) +
#  theme_legend +
#  theme(axis.text.y = element_blank(), axis.line = element_blank(),
#        axis.ticks.y = element_blank(), panel.spacing.x = unit(0.2, "lines"), 
#        legend.position = "right", panel.grid = element_blank(), panel.grid.major = element_blank(), legend.text = element_text(size=6), axis.text.x=element_blank(), axis.title.x=element_blank()
#        )

##layout <- c(area(1,1), area(1,2,2,2), area(2,1), area(1,3,2,3), area(3,1,3,2), area(3,3))
#layout<- c(area(1,1,1,3), area(1,3,3,4), area(2,1), area(1,3,2,3), area(3,1,3,2), area(3,3))
#layout<- c(area(1,1,1,2), area(1,3,3,4), area(2,1,3,1))#, area(2,2,3,2))
#for heat,apgene==area(2,5,,5),
#+heatmap_gene

#jpeg("../Figures/qPCR_validation_plot.jpg", width = 1800, height = 900, units = "px")
#qPCR_val_DTU
#dev.off()

#------clock-isoform-structures----#
#cry1<-gene_structure("ENSMUSG00000020038") #cry1:twotranscripts
#per1<-gene_structure("ENSMUSG00000020893") #per1
#Per2<-gene_structure("ENSMUSG00000055866")#per2
#Ciart<-gene_structure("ENSMUSG00000038550") #ciart
#dbp<-gene_structure("ENSMUSG00000059824") #dbp
#Hlf<-gene_structure("ENSMUSG00000003949") #Hlf
#arntl<-gene_structure("ENSMUSG00000015522") #arntl
#cry2<-gene_structure("ENSMUSG00000068742") #cry2
#clock<-gene_structure("ENSMUSG00000029238") #clock
#nr1d1<-gene_structure("ENSMUSG00000020889") #nr1d1
#nr1d2<-gene_structure("ENSMUSG00000021775") #nr1d2
#rorc<-gene_structure("ENSMUSG00000028150") #rorc
#noct<-gene_structure("ENSMUSG00000023087") #noct
#bhlhe40<-gene_structure("ENSMUSG00000030103") #bhlhe40
#nampt<-gene_structure("ENSMUSG00000020572") #nampt
#ggsave("../extra/cry1.png", plot=cry1, width=8, height = 2)
#ggsave("../extra/per1.png", plot=per1, width=8, height = 2)
#ggsave("../extra/Per2.png", plot=Per2, width=8, height = 2)
#ggsave("../extra/Ciart.png", plot=Ciart, width=8, height = 2)
#ggsave("../extra/dbp.png", plot=dbp, width=8, height = 2)
#ggsave("../extra/Hlf.png", plot=Hlf, width=8, height = 2)
#ggsave("../extra/arntl.png", plot=arntl, width=8, height = 2)
#ggsave("../extra/cry2.png", plot=cry2, width=8, height = 2)
#ggsave("../extra/clock.png", plot=clock, width=8, height = 2)
#ggsave("../extra/nr1d1.png", plot=nr1d1, width=8, height = 2)
#ggsave("../extra/nr1d2.png", plot=nr1d2, width=8, height = 2)
#ggsave("../extra/rorc.png", plot=rorc, width=8, height = 2)
#ggsave("../extra/noct.png", plot=noct, width=8, height = 2)
#ggsave("../extra/bhlhe40.png", plot=bhlhe40, width=8, height = 2)
#ggsave("../extra/nampt.png", plot=nampt, width=8, height = 2)
#dtu_structures<- sapply(dtu_genes, gene_structure, simplify = FALSE, USE.NAMES = TRUE)
#clock_structures<- sapply((annotation %>%  filter(mgi_symbol %in% clock_genes) %>% dplyr::pull(gene_id) %>% unique()), gene_structure, simplify = FALSE, USE.NAMES = TRUE)
#plot_grid(dtu_structures$ENSMUSG00000003721,dtu_structures$ENSMUSG00000021578,dtu_structures$ENSMUSG00000015357,dtu_structures$ENSMUSG00000024070, nrow=2, ncol=2)
#same for other genes
# Fgg: ENSMUST00000048486 & ENSMUST00000194175 - ENSMUSG00000033860
# Pxmp4 : ENSMUST00000000896 & ENSMUST00000109703 - ENSMUSG00000000876
# 0610031O16Rik : ENSMUST00000185046 & ENSMUST00000185122 - ENSMUSG00000099146
# Rrbp1 : ENSMUST00000016072 & ENSMUST00000037875 - ENSMUSG00000027422
### Gathering
# qPCR validation of new clock gene isoforms
# heatmap_plot
# phase_distribution
# time_course_plots_WT
# time_course_plots_DTU
# qPCR_val_DTU
# iso structures: Fgg_structure Rik_structure Pxmp4_structure Rrbp1_structure