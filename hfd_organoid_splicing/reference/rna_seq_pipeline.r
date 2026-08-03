### RNA Seq exploratory analysis on organoide

library(DESeq2)
library(tximeta)
library(readr)
library(pheatmap)
library(RColorBrewer)
library(org.Mm.eg.db)
library(clusterProfiler)
library(forcats)
library(cowplot)
library(RColorBrewer)
library(apeglm)
library(limma)
library(ashr)


### Load samples

working_dir <- getwd()
coldata <- read.csv(paste0(working_dir, "/data/metafile.csv"))
coldata$files <- file.path(working_dir, "results/salmon",
    coldata$samples, "quant.sf")
coldata$names <- coldata$samples

se <- tximeta(coldata)
gse <- summarizeToGene(se)


### Filter genes with very low counts accross all samples
dds <- DESeqDataSet(gse, design = ~ conditions)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

dds$conditions <- relevel(factor(dds$conditions), ref = "WT_SD")
dds$conditions <- factor(dds$conditions, 
    levels = c("WT_SD", "W8_WHFD", "W18_WHFD", "W24_WHFD", "W42_WHFD"))

vsd <- vst(dds, blind = TRUE)
pcaplot <- plotPCA(vsd, intgroup = "conditions") +
    ggplot2::geom_point(size = 8) +
    ggplot2::theme_minimal() + 
    ggplot2::theme(
        plot.title = ggplot2::element_text(size = 24),
        axis.title.x = ggplot2::element_text(size = 20),
        axis.text.x = ggplot2::element_text(size = 20),
        axis.title.y = ggplot2::element_text(size = 20),
        axis.text.y = ggplot2::element_text(size = 20),
        legend.text = ggplot2::element_text(size = 20),
        legend.title = ggplot2::element_text(size = 22)
    ) +
    ggplot2::ggtitle("PCA")
ggplot2::ggsave(file.path(working_dir, "results/differential_exprs/pcaplot.png"),
    pcaplot, bg = "white", width = 12, height = 12)

############################## SAFE mRNA MATRIX FOR COFE
vst_mat <- assay(vsd)
gene_symbols <- mcols(dds)$symbol
rownames(vst_mat) <- make.unique(ifelse(is.na(gene_symbols), 
                                         rownames(vst_mat), 
                                         gene_symbols))

head(vst_mat)
col_names <- c("WT_SD_1", "WT_SD_2", "WT_SD_3",
    "W8_WHFD_1", "W8_WHFD_2", "W8_WHFD_3", "W18_WHFD_1", "W18_WHFD_2", "W18_WHFD_3",
    "W24_WHFD_1", "W24_WHFD_2", "W24_WHFD_3", "W42_WHFD_1", "W42_WHFD_2", "W42_WHFD_3")
colnames(vst_mat) <- col_names

col_names <- colnames(vst_mat)
conditions <- str_match(col_names, "^(.+)_\\d+$")[, 2]
metadata_pred <- data.frame(
  sample_id = col_names,
  time = NA,              
  condition = conditions,
  train = FALSE,          
  row.names = col_names
)

ddds <- estimateSizeFactors(dds)
normalized_counts <- counts(ddds, normalized = TRUE)
query_counts_log <- log2(normalized_counts + 1)
colnames(query_counts_log) <- col_names
rownames(query_counts_log) <- make.unique(ifelse(is.na(gene_symbols), 
    rownames(query_counts_log), gene_symbols))
head(query_counts_log)

write.csv(query_counts_log, 
          file.path(working_dir, "COFE/rna_matrix.csv"), row.names = TRUE)
write.csv(metadata_pred, 
          file.path(working_dir, "COFE/rna_metadata.csv"), 
          row.names = FALSE)


##############################################################

################# Differential Expression #######################
### map genes names to ensembl names
mcols(dds)$symbol <- mapIds(org.Mm.eg.db,
    keys = rownames(dds),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first")

### Statical approaches 
# Deseq2 uses negative binomial distribution to estimate dispersion.
# Wald test fits a statistical model between two groups to estimate the logFC.
# find if it's statistically != 0, finds a signal to noise by dividing logFC by
# standard error to get zscore then finds the pvalue of that gene. Basically, It
# test the specific comparision : does group A differ from group B?

# Likelihood ration test fits two models to data and not each gene like Wald test.
# I mean takes a gene and fit accross all conditions. reduce the model to the conditions
# to try to find if there are conditions that explains the data/gene very well than the 
# bigger model. if not then that gene does not respond to those conditions or treatment.
### To keep it short, it ask the questions after accounting/removing for conditions,
### does the data still mater?  it's tessting the effect of what is removed from the model.
design(dds) <- ~ conditions

#### Likelihood Ratio TEst 
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~1)
res_lrt <- results(dds_lrt)
### the output is genes that vary 
# significantly in at least one condition compared to the others. But it 
# doesn't tell you which condition is different because any pattern is possible.


##### Pairwise comparison with Wald test
dds <- DESeq(dds)
resultsNames(dds)

str(dds, 2)
dds$conditions
resultsNames(dds)
# extract each pairwise comparision
# group_list <- list(
#     "res_8w" = "conditions_W8_WHFD_vs_WT_SD",
#     "res_18" = "conditions_W18_WHFD_vs_WT_SD",
#     "res_24w" = "conditions_W24_WHFD_vs_WT_SD",
#     "res_42w" = "conditions_W42_WHFD_vs_WT_SD"
# )

# res <- lapply(names(group_list), function(names) {
#     r <- lfcShrink(dds, coef=group_list[[names]], type="apeglm")
#     r$symbol <- mcols(dds)$symbol
#     return(r)
# }) 
# names(res) <- names(group_list)
# sig <- lapply(res, function(res_group) {
#     subset(res_group, padj < 0.05 & abs(log2FoldChange) > 1)
# })

###### Add for custom pairwise comparision
res_18vs8  <- lfcShrink(dds, 
    contrast = c("conditions", "W18_WHFD", "W8_WHFD"),
    type     = "ashr")
res_24vs18 <- lfcShrink(dds,
    contrast = c("conditions", "W24_WHFD", "W18_WHFD"),
    type     = "ashr")
res_42vs24 <- lfcShrink(dds,
    contrast = c("conditions", "W42_WHFD", "W24_WHFD"),
    type     = "ashr")
#### to ease workflow, keep the same settings for non custom as above
res_8w  <- lfcShrink(dds, coef = "conditions_W8_WHFD_vs_WT_SD",  type = "apeglm")
res_18w <- lfcShrink(dds, coef = "conditions_W18_WHFD_vs_WT_SD", type = "apeglm")
res_24w <- lfcShrink(dds, coef = "conditions_W24_WHFD_vs_WT_SD", type = "apeglm")
res_42w <- lfcShrink(dds, coef = "conditions_W42_WHFD_vs_WT_SD", type = "apeglm")

# ### With limma, modify the above method to not only include the wt vs others comparision but also
# ### make a pairwise comparision between non wildtype experiments
# condition <- factor(colData(dds)$conditions,
#     levels = c("WT_SD", "W8_WHFD", "W18_WHFD", "W24_WHFD", "W42_WHFD"))
# design <- model.matrix(~ 0 + condition)
# colnames(design) <- levels(condition)
# mat <- assay(vsd)
# # head(mat)
# # colnames(vsd)
# # colnames(dds)
# dim(mat)
# dim(design)
# fit <- lmFit(mat, design)
# fit$design
# contrast_matrix <- makeContrasts(
#     W8_vs_WT   = `W8_WHFD`  - WT_SD,
#     W18_vs_WT  = `W18_WHFD` - WT_SD,
#     W24_vs_WT  = `W24_WHFD` - WT_SD,
#     W42_vs_WT  = `W42_WHFD` - WT_SD,
#     W18_vs_W8  = `W18_WHFD` - `W8_WHFD`,
#     W24_vs_W18 = `W24_WHFD` - `W18_WHFD`,
#     W42_vs_W24 = `W42_WHFD` - `W24_WHFD`,
#     levels     = design
# )
# colnames(design)
# rownames(contrast_matrix)
# fit2 <- contrasts.fit(fit, contrast_matrix)
# fit2 <- eBayes(fit2)


# # check if results make sense
# res <- lapply(colnames(contrast_matrix), function(coef) {
#     topTable(fit2, coef = coef, number = Inf, adjust.method = "BH")
# })
# names(res) <- colnames(contrast_matrix)

# # check significant genes per comparison
# sapply(res, function(r) {
#     nrow(subset(r, adj.P.Val < 0.05 & abs(logFC) > 1))
# })
# str(res, 1)
# names(res)
# head(res$W8_vs_WT)


# ### add gene symbol
# res <- lapply(res, function(r) {
#     r$symbol <- mapIds(org.Mm.eg.db,
#                        keys     = rownames(r),
#                        column   = "SYMBOL",
#                        keytype  = "ENSEMBL",
#                        multiVals = "first")
#     # rename logFC to log2FoldChange for EnhancedVolcano
#     colnames(r)[colnames(r) == "logFC"] <- "log2FoldChange"
#     colnames(r)[colnames(r) == "adj.P.Val"] <- "padj"
#     return(r)
# })
### After trying limma, I choose to proceed with Deseq2 because it has a more stable    
### variation of dge based on logfoldchange significance and genes abundance

###### Visualizations
res <- list(
    "res_8w"    = res_8w,
    "res_18"    = res_18w,
    "res_24w"   = res_24w,
    "res_42w"   = res_42w,
    "res_18vs8"  = res_18vs8,
    "res_24vs18" = res_24vs18,
    "res_42vs24" = res_42vs24
)
head(res[[1]])
### add gene symbol
res <- lapply(res, function(r) {
    r$symbol <- mapIds(org.Mm.eg.db,
        keys     = rownames(r),
        column   = "SYMBOL",
        keytype  = "ENSEMBL",
        multiVals = "first")
    return(r)
})

# Volcano plot
sig <- lapply(res, function(res_group) {
    subset(res_group, padj < 0.05 & abs(log2FoldChange) > 1)
})

lapply(names(res), function(names) {
    r <- as.data.frame(res[[names]])
    r$gene <- r$symbol

    # pick top 10 genes to label
    top <- r[!is.na(r$padj), ]
    top <- top[order(top$padj), ]
    top <- head(top$gene, 10)

    res_volcan <- EnhancedVolcano::EnhancedVolcano(res[[names]],
        lab = "",#res[[names]]$symbol,
        selectLab = top,
        x = "log2FoldChange",
        y = "padj",
        pCutoff = 0.05,
        FCcutoff = 1,
        pointSize = 2, 
        labSize = 8,
        drawConnectors = TRUE,
        max.overlaps = Inf,
        title = paste("Volcano Plot :", names)) +#,
    #     subtitle = paste0(nrow(sig[[names]]), " significant genes (padj<0.05, |LFC|>1)")) +
    ggplot2::theme(
        axis.title = ggplot2::element_text(size = 24),
        axis.text.x = ggplot2::element_text(size = 24),
        axis.text.y = ggplot2::element_text(size = 24),
        plot.title = ggplot2::element_text(size = 28),
        plot.subtitle = ggplot2::element_text(size = 24),
        text = ggplot2::element_text(size = 24)
    )
    ggplot2::ggsave(file.path(working_dir, paste0("results/differential_exprs/limma/no_labels/volcan_",
        names, ".png")),
        res_volcan, bg = "white", width = 12, height = 12)
})
## uncomment and change save directory for plots with 0 labels.

### Make MA plots
# lapply(names(res), function(names) {
#     png(file.path(working_dir, 
#         paste0("results/differential_exprs/limma/ma_", names, ".png")), 
#         width = 800, height = 800)
#     plotMA(res[[names]], ylim = c(-10, 10), main = paste("MA Plot :", group_list[[names]]))
#     dev.off()
# })
head(res[[1]])
lapply(names(res), function(names) {
    png(file.path(working_dir,
        paste0("results/differential_exprs/limma/ma_", names, ".png")),
        width = 800, height = 800)
    
    par(mar = c(5, 6, 4, 2))
    
    plot(x    = log2(res[[names]]$baseMean + 1),
        y    = res[[names]]$log2FoldChange,
        ylim = c(-4, 4),
        xlab = "log2 mean expression",
        ylab = "log2 fold change",
        main = paste("MA Plot :", names),
        pch  = 20,
        cex  = 0.7,
        cex.main = 3,
        cex.axis = 3,
        cex.lab = 3,
         col  = ifelse(res[[names]]$padj < 0.05, "blue", "grey"))
    
    abline(h = 0, col = "black", lwd = 2)
    dev.off()
})

### make a pheatmap/GO plot from pairwise statistics
go_enrichment <- function(cluster_df, clust_no) {
    genes_smb <- cluster_df[cluster_df$cluster == clust_no, "genes"]
    id <- bitr(genes_smb, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
    genes <- id$ENTREZID
    print(genes)
    enrich <- enrichGO(gene = genes,
        OrgDb = org.Mm.eg.db,  
        ont = "BP",
        minGSSize = 1,
        pvalueCutoff = 0.05,
        qvalueCutoff = 1,
        readable = TRUE)

    pp <- enrich@result %>%
        filter(p.adjust < 0.05) %>%
        dplyr::slice_min(p.adjust, n = 15) %>%
        dplyr::mutate(Description = fct_reorder(Description, Count))
    return(pp)
}

keggenrich_list <- function(cluster_df, clust_no) {
    genes_smb <- cluster_df[cluster_df$cluster == clust_no, "genes"]
    id <- bitr(genes_smb, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    genes <- id$ENTREZID
    print(paste("Working on KEGG genes:", length(genes)))
    if (length(genes) <= 1) {return(NULL)}
    
    # You already converted to ENTREZID above, don't do it again
    ekegg <- enrichKEGG(
        gene = genes,
        organism = "mmu",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
    )
    if (is.null(ekegg) || nrow(ekegg) == 0) {return(NULL)}
    
    top_ekegg <- ekegg@result %>%
        dplyr::filter(p.adjust < 0.05) %>%
        dplyr::slice_min(p.adjust, n = 15)
    print(top_ekegg$Description)
    
    if (nrow(top_ekegg) == 0) {
        message(paste0("  No significant KEGG terms — skipping"))
        return(NULL)
    }
    return(top_ekegg)
}
print(group_list)
print(names(group_list))
print(names(res))
str(res, 1)
group_list <- list(
    "res_8w"     = "conditions_W8_WHFD_vs_WT_SD",
    "res_18"     = "conditions_W18_WHFD_vs_WT_SD",
    "res_24w"    = "conditions_W24_WHFD_vs_WT_SD",
    "res_42w"    = "conditions_W42_WHFD_vs_WT_SD",
    "res_18vs8"  = "conditions_W18_WHFD_vs_W8_WHFD",
    "res_24vs18" = "conditions_W24_WHFD_vs_W18_WHFD",
    "res_42vs24" = "conditions_W42_WHFD_vs_W24_WHFD"
)
lapply(names(res), function(names) {
    # rows
    siggenes <- subset(res[[names]], padj < 0.05 & abs(log2FoldChange) > 1)
    siggenes <- siggenes[order(siggenes$padj), ]
    # cols
    # cond <- gsub("conditions_(.*)_vs_WT_SD", "\\1", group_list[[names]])
    # keep_samples <- colData(vsd)$conditions %in% c("WT_SD", cond)
    # detect if comparison is vs WT or progressive
    if (grepl("vs_WT_SD", group_list[[names]])) {
        # vs WT comparison
        cond <- gsub("conditions_(.*)_vs_WT_SD", "\\1", group_list[[names]])
        keep_samples <- colData(vsd)$conditions %in% c("WT_SD", cond)
    } else {
        # progressive pairwise
        cond1 <- gsub("conditions_(.*)_vs_(.*)", "\\1", group_list[[names]])
        cond2 <- gsub("conditions_(.*)_vs_(.*)", "\\2", group_list[[names]])
        keep_samples <- colData(vsd)$conditions %in% c(cond1, cond2)
    }
    # subset rows and cols
    mat <- assay(vsd)[rownames(siggenes), keep_samples]
    print(head(mat))

    gene_symbols <- mapIds(org.Mm.eg.db,
        keys = rownames(mat),
        column = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first")
    rownames(mat) <- make.unique(ifelse(is.na(gene_symbols), rownames(mat), gene_symbols))

    # get cluster anotation before anootating
    pheatoplot <- pheatmap(mat,
        scale = "row",
        clustering_method = "ward.D2",
        cutree_rows = 2,
        cluster_cols = TRUE,
        silent = FALSE)

    # Draw second pheatmap and annotate cluster columns
    clusters <- cutree(pheatoplot$tree_row, k = 2)
    cluster_labels <- paste("Cluster", clusters)
    cluster_anno <- data.frame(Cluster = factor(cluster_labels, 
        levels = paste("Cluster", 1:2)))
    rownames(cluster_anno) <- names(clusters)
    cluster_colors <- list(Cluster = setNames(
        brewer.pal(2, "Set1"), 
        paste("Cluster", 1:2)))
    par(mar = c(5, 7, 4, 2), mgp = c(4, 1, 0))
    pheatoplot <- pheatmap(mat,
        scale = "row",
        annotation_col = as.data.frame(colData(vsd)[keep_samples, 
            "conditions", drop = FALSE]),
        annotation_row = cluster_anno,
        annotation_colors = cluster_colors,
        cluster_cols = TRUE,
        show_rownames = FALSE,
        show_colnames = FALSE,
        dendrogram = "none",
        fontsize = 24,
        fontsize_col = 20,
        fontsize_row = 20,
        cellwidth = 25,
        clustering_method = "complete",
        #clustering_distance_rows = "canberra",
        main = paste("Pairwise significant DE genes\n", group_list[[names]],
            "\n",nrow(mat), "genes,", "padj < 0.05, |LFC| > 1"))
    filename = file.path(working_dir,
        paste0("results/differential_exprs/sgnft_pheatoplot_", names, ".png"))
    ggplot2::ggsave(filename, pheatoplot, bg = "white", width = 11.5, height = 15)

    ############ GO enrichment
    ## Get genes and plot functional annotations along side.
    clusters <- cutree(pheatoplot$tree_row, k = 2)
    cluster_df <- data.frame(genes = names(clusters), cluster = clusters)
    write.csv(cluster_df,
        file.path(working_dir, paste0("results/differential_exprs/limma/pw_genes_", names, ".txt")))

    # Plot go enrichment for both clusters
    plot_list <- lapply(seq(2), function(cn) {
        go_res <- go_enrichment(cluster_df, cn)
        mapped_genes <- length(cluster_df[cluster_df$cluster == cn, "genes"])
        pp <- ggplot2::ggplot(go_res, ggplot2::aes(x = Count, y = Description, fill = p.adjust)) +
            ggplot2::geom_col() +
            ggplot2::scale_fill_gradient(low = "firebrick3", high = "steelblue") +
            ggplot2::labs(
                title = paste("GO - Cluster", cn, ": ", group_list[[names]], "\nmapped:", mapped_genes),
                x     = "Gene Count",
                y     = NULL,
                fill  = "adj. p-value") +
            ggplot2::theme_bw() +
            ggplot2::theme(
                plot.margin = ggplot2::margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
                plot.title  = ggplot2::element_text(size = 24, face = "bold"),
                axis.text.y = ggplot2::element_text(size = 24),
                axis.text.x = ggplot2::element_text(size = 24),
                axis.title = ggplot2::element_text(size = 24),
                legend.title = ggplot2::element_text(size = 24),
                legend.text = ggplot2::element_text(size = 20)

            )
        ggplot2::ggsave(paste0(working_dir, "/results/differential_exprs/pw_go_", names, "_", cn, ".png"), 
            pp, bg = "white", width = 20, height = 10)
        return (pp)
    })

    plot_list_kegg <- lapply(seq(2), function(cn) {
        mapped_genes <- length(cluster_df[cluster_df$cluster == cn, "genes"])
        top_ekegg <- keggenrich_list(cluster_df, cn)
        top_ekegg <- top_ekegg %>%
            dplyr::arrange(desc(Count)) %>%
            dplyr::mutate(Description = factor(Description, levels = rev(Description)))

        pkegg <- ggplot2::ggplot(top_ekegg, ggplot2::aes(x = Count, y = Description, fill = p.adjust)) +
            ggplot2::geom_col(width = 0.6) +
            ggplot2::scale_fill_gradient(low = "firebrick3", high = "steelblue", name = "adj. p-value") +
            ggplot2::labs(
                title = paste("KEGG : ", group_list[[names]], "\nCluster", cn, " : mapped - ", mapped_genes),
                x = "Gene Count",
                y = NULL
            ) +
            ggplot2::theme_bw() +
            ggplot2::theme(
                plot.margin = ggplot2::margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
                plot.title  = ggplot2::element_text(size = 24, face = "bold"),
                axis.text.y = ggplot2::element_text(size = 24),
                axis.text.x = ggplot2::element_text(size = 24),
                axis.title = ggplot2::element_text(size = 24),
                legend.title = ggplot2::element_text(size = 24),
                legend.text = ggplot2::element_text(size = 20)
            )
        ggplot2::ggsave(paste0(working_dir, "/results/differential_exprs/pw_kegg_", names, "_", cn, ".png"), 
            pkegg, bg = "white", width = 18, height = 10)
    return(pkegg)
    })

    len <- length(plot_list)
    tt <- plot_grid(
        plotlist = plot_list[1:len],
        ncol = 1,
        align = "v",
        axis = "tb"
    )

    len <- length(plot_list_kegg)
    kk <- plot_grid(
        plotlist = plot_list_kegg[1:len],
        ncol = 1,
        align = "v",
        axis = "tb"
    )

    pheato_ego <- plot_grid(
        pheatoplot$gtable,
        tt,
        kk,
        nrow = 1,
        align = "v",
        axis = "l",
        rel_widths = c(0.2, 0.4, 0.4)
    )
    ggplot2::ggsave(paste0(working_dir, "/results/limma/differential_exprs/limma/pw_go_", names, ".png"), 
            pheato_ego, bg = "white", width = 32, height = 22)
})



######################################## LIKELIHOOD RATIO TEST

# Make pheatmap for Likelihood ratio stats test method
top_gene_n <- 150
sgnfct_genes <- subset(res_lrt, padj < 0.05 & abs(log2FoldChange) > 1.5)
sgnfct_genes <- sgnfct_genes[order(sgnfct_genes$padj), ]
mat_vsd <- assay(vsd)[head(rownames(sgnfct_genes), top_gene_n), ]

gsymbols <- mapIds(org.Mm.eg.db,
    keys = rownames(mat_vsd),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first")
rownames(mat_vsd) <- ifelse(is.na(gsymbols), rownames(mat_vsd), gsymbols)
# Draw first pheatmap to get clusters for annotation
bb1 <- pheatmap(mat_vsd,
    scale = "row",
    clustering_method = "average",
    clustering_distance_rows = "euclidean",
    cutree_rows = 9,
    cluster_cols = FALSE,
    silent = TRUE)

# Draw second pheatmap and annotate cluster columns
clusters <- cutree(bb1$tree_row, k = 9)
cluster_labels <- paste("Cluster", clusters)
cluster_anno <- data.frame(Cluster = factor(cluster_labels, 
    levels = paste("Cluster", 1:9)))
rownames(cluster_anno) <- names(clusters)
library(RColorBrewer)
cluster_colors <- list(Cluster = setNames(
    brewer.pal(9, "Set1"), 
    paste("Cluster", 1:9)))

bb1 <- pheatmap(mat_vsd,
    scale = "row",
    annotation_col = as.data.frame(colData(vsd)[, "conditions", drop = FALSE]),
    annotation_row = cluster_anno,
    annotation_colors = cluster_colors,
    show_colnames = FALSE,
    fontsize = 12,
    clustering_method = "average",
    clustering_distance_rows = "euclidean",
    cutree_rows = 9,
    cluster_cols = FALSE,
    main = paste("LRT: Genes changing across HFD",
        "\n", nrow(mat_vsd), "set1 top genes,", "padj < 0.05, |LFC| > 1.5"))
ggplot2::ggsave(file.path(working_dir, "results/differential_exprs/limma/pheatomap_lrt1.png"),
    bb1, bg = "white", width = 12, height = 24)

# test if clusters match by mapping manually gene counts
clusters <- cutree(bb1$tree_row, k = 9)
cluster_df <- data.frame(genes = names(clusters), cluster = clusters)
cluster_df[cluster_df$cluster == 2, "genes"]

############ GO enrichment
## Get genes and plot functional annotations along side.
clusters <- cutree(bb1$tree_row, k = 9)
cluster_df <- data.frame(genes = names(clusters), cluster = clusters)

plot_list <- lapply(seq(9), function(cn) {
    go_res <- go_enrichment(cluster_df, cn)
    mapped_genes <- length(cluster_df[cluster_df$cluster == cn, "genes"])
    pp <- ggplot2::ggplot(go_res, ggplot2::aes(x = Count, y = Description, fill = p.adjust)) +
        ggplot2::geom_col() +
        ggplot2::scale_fill_gradient(low = "firebrick3", high = "steelblue") +
        ggplot2::labs(
            title = paste("Cluster", cn, "\nmapped:", mapped_genes),
            x     = "Gene Count",
            y     = NULL,
            fill  = "adj. p-value") +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.title  = ggplot2::element_text(size = 16, face = "bold"),
            axis.text.y = ggplot2::element_text(size = 16)
        )
    ggplot2::ggsave(paste0(working_dir, "/results/differential_exprs/limma/go_", cn, ".png"), 
        pp, bg = "white", width = 12, height = 12)
    return (pp)
})

### Combine go and pheatmap plot
len <- length(plot_list)
tt <- plot_grid(
    plotlist = plot_list[1:len],
    ncol = 2,
    align = "v",
    axis = "tb"
)

pheato_ego <- plot_grid(
    bb1$gtable,
    tt,
    nrow = 1,
    align = "v",
    axis = "l",
    rel_widths = c(0.2, 0.6)
)
ggplot2::ggsave(paste0(working_dir, "/results/differential_exprs/ttgo_", ".png"), 
        pheato_ego, bg = "white", width = 34, height = 22)




############# READ REF CIRCADIAN TIME AND TRANFER TO COFE
ref_rna_seq <- read.delim(file.path(working_dir, "data/GSE117134_casava_gene_expression.tsv"))
head(ref_rna_seq)
rownames(ref_rna_seq) <- ref_rna_seq$Gene
ref_counts <- ref_rna_seq[, -1]
ref_counts <- log2(ref_counts + 1)
head(ref_counts)
colnames(ref_counts)
ref_counts_wt <- ref_counts[, 1:18]
head(ref_counts_wt)
dim(ref_counts_wt)

sample_names <- colnames(ref_counts_wt)
library(stringr)
zt_hours <- as.numeric(str_match(sample_names, "ZT(\\d+)")[, 2])
conditions <- str_match(sample_names, "^([A-Za-z0-9_\\.]+)_ZT")[, 2]
metadata <- data.frame(
  sample_id = sample_names,
  time = zt_hours,
  condition = conditions,
  train = TRUE, # Marking this reference data as training
  row.names = sample_names
)
head(metadata)

write.csv(ref_counts_wt, "data/ref_expression_raw.csv", row.names = TRUE)
write.csv(metadata, "data/ref_metadata.csv", row.names = FALSE)




######################## TIME TELLER #########################
library(TimeTeller)
help(train_model)


cn <- colnames(ref_counts_wt)
time      <- as.numeric(sub(".*ZT0*(\\d+).*", "\\1", colnames(ref_counts_wt))) 
replicate <- sub(".*_BR(\\d+)$", "\\1", cn)                
group     <- rep("WT_liver", length(cn))                   

panda_data$genes_used
clock_genes <- intersect(panda_data$genes_used, rownames(ref_counts_wt))
## missing gene : "Rorb" 

range(panda_data$expr_mat)
range(ref_counts_wt)


tt_model <- train_model(
    exp_matrix = ref_counts_wt, 
    genes = clock_genes,
    group_1 = group, 
    time = time, 
    replicate = replicate, 
    mat_normalised = TRUE,
    cores = 24, 
    method = "intergene",
    num_PC = 2,
    log_thresh = -5)

?test_model

########## Get test sample data for projection
library(edgeR)
counts <- counts(dds, normalized = FALSE) 
counts_logcpm <- cpm(counts, log = TRUE, prior.count = 2)

gene_symbols <- mcols(dds)$symbol
rownames(counts_logcpm) <- make.unique(ifelse(is.na(gene_symbols),
    rownames(counts_logcpm), gene_symbols))
colnames(counts_logcpm) <- col_names
head(counts_logcpm)
dim(counts_logcpm)

range(ref_counts_wt) 
range(vst_mat)

cn <- colnames(counts_logcpm)
test_group <- sub("_[0-9]+$", "", cn)
test_replicate <- sub(".*_([0-9]+)$", "\\1", cn)

confirm_test_genes <- intersect(clock_genes, rownames(counts_logcpm))
rownames(counts_logcpm)[rownames(counts_logcpm) == "Bmal1"] <- "Arntl"

test_time <- rep(NA, length(cn))

tt_result <- test_model(
    object = tt_model,
    exp_matrix = counts_logcpm,
    test_group_1 = test_group,
    test_time = test_time,
    test_replicate = test_replicate, 
    mat_normalised = TRUE,
    cores = 24, 
    log_thresh = -5)

str(tt_result, max.level = 2)    
names(tt_result) 

tt_result$Test_Data$Results_df       
str(tt_result$Test_Data$Thetas_Test)   
