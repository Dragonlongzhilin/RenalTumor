#' @description: this script idnetify TF target gene and construct the TF network

library(openxlsx)
library(circlize)
library(dplyr)
library(tibble)
library(data.table)
library(Signac)
library(Seurat)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVARmotifs)
data("human_pwms_v2")
library(ggplot2)
library(ggseqlogo)
library(ggpubr)
library(openxlsx)

setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC")
scATAC.data <- readRDS("scATAC.data.pro.rds")
DefaultAssay(scATAC.data) <- "Peaks"
scRNA.data <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/data.merge.pro.rds")
DefaultAssay(scRNA.data) <- "RNA"

celltype <- "Tumor"
Tumor.TF.motifs <- readRDS("5.Motif/Analysis/tumor.specific.TFs.rds")
# scRNA.DEGs <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/2.Cluster/AnnotateCellType/cellType.sig.pos.DEGs.rds")
scRNA.all.DEGs <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/2.Cluster/AnnotateCellType/cellType.all.DEGs.rds")
DEGs <- scRNA.all.DEGs %>% filter(cluster == "Tumor" & p_val_adj < 0.05 & abs(avg_log2FC) > 0.25) %>% .$gene # 1451

# motif information
motif.info <- data.frame(originName = names(scATAC.data@assays$Peaks@motifs@motif.names), TF = unlist(scATAC.data@assays$Peaks@motifs@motif.names))
rownames(motif.info) <- NULL
#                           originName     TF
#   ENSG00000008196_LINE2_TFAP2B_D_N1 TFAP2B
#      ENSG00000008197_LINE6_TFAP2D_D TFAP2D
#   ENSG00000087510_LINE7_TFAP2C_D_N3 TFAP2C
#     ENSG00000116819_LINE19_TFAP2E_I TFAP2E
#  ENSG00000137203_LINE29_TFAP2A_D_N3 TFAP2A
#  ENSG00000116017_LINE44_ARID3A_D_N1 ARID3A

# subset by celltype
scATAC.sub <- subset(scATAC.data, AnnotatedcellType == celltype)
scRNA.sub <- subset(scRNA.data, cellType_low == celltype)

scATAC.sub <- FindTopFeatures(scATAC.sub, min.cutoff=ncol(scATAC.sub)*0.05) # present in what % of cells
scRNA.sub <- FindVariableFeatures(scRNA.sub, nfeatures = 3000)

# use chipseeker annotate each peak in scATAC-seq
peak.info <- readRDS("4.Peak/peak.annotation.simple.ChIPseeker.rds")

PWMatrixToProbMatrix <- function(x){
        if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
        m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
        m <- t(t(m)/colSums(m))
        m
}

cCREs <- readRDS(paste0("6.Co-Accessible/", celltype,"_peak_gene_correlation.rds")) # 4368 cCREs
#gene-cCRE link --- enhancer
cCREs <- cCREs %>% subset(distance_bp/1000>=10 & peak2_type != 'Promoter')
cCREs$target_gene <- cCREs$peak1_nearestGene
cCREs$peak_gene <- paste0(cCREs$Peak2, '_', cCREs$peak1_nearestGene)
cCREs$Peak1 <- gsub("_", "-", cCREs$Peak1)
cCREs$Peak2 <- gsub("_", "-", cCREs$Peak2)

################################################################################
# Find target genes for each TFs
################################################################################
# settings for featureplot
order_values <- TRUE
reduction <- 'umap'

# for each motif, find genes:
use_variable_genes <- T # require the variable gene

res <- sapply(Tumor.TF.motifs$Name, function(x){
    motif_name <- x
    motif_ID <- motif.info$originName[which(motif.info$TF==x)] # e.g., ENSG00000008196_LINE2_TFAP2B_D_N1

    #### get list of promoter and enhancer targets for these TFs
    ## target genes with promoter type
    DefaultAssay(scATAC.sub) <- "Peaks"
    motif_accessible <- Motifs(scATAC.sub)@data[, motif_ID]
    motif_accessible <- names(motif_accessible)[motif_accessible > 0]
    motif_accessible <- motif_accessible[motif_accessible %in% VariableFeatures(scATAC.sub)]

    # which of these peaks are at promoters?
    motif_accessible_promoters <- motif_accessible[motif_accessible %in% peak.info$peaks[which(peak.info$originType == "Promoter (<=1kb)")]]
    # which genes are associated with these promoters?
    motif_promoter_target_genes <- peak.info$SYMBOL[match(motif_accessible_promoters, peak.info$peaks)]
    motif_promoter_target_genes <- as.character(motif_promoter_target_genes)

    # variable genes only? which of these genes are highly expressed in snRNA-seq?
    if(use_variable_genes){
      motif_promoter_target_genes <- motif_promoter_target_genes[motif_promoter_target_genes %in% VariableFeatures(scRNA.sub)] %>% as.character %>% unique
    } else{
      motif_promoter_target_genes<- motif_promoter_target_genes %>% as.character %>% unique
    }

    ## target genes with enhancer type
    motif_accessible_enhancers <- motif_accessible[motif_accessible %in% cCREs$Peak2]
    motif_cCREs <- subset(cCREs, Peak2 %in% motif_accessible_enhancers)
    motif_enhancer_target_genes <- motif_cCREs$target_gene

    # variable genes only? which of these genes are highly expressed in snRNA-seq?
    if(use_variable_genes){
      motif_enhancer_target_genes <- motif_enhancer_target_genes[motif_enhancer_target_genes %in% VariableFeatures(scRNA.sub)] %>% as.character %>% unique
    }else{
      motif_enhancer_target_genes<- motif_enhancer_target_genes %>% as.character %>% unique
    }

    # check if there are promoter targets:
    if(length(motif_promoter_target_genes) > 0){
      promoter_edge_df <- data.frame(
        from = motif_name,
        to = as.character(motif_promoter_target_genes),
        type = 'promoter'
      )
    } else{promoter_edge_df <- data.frame()}

    # check if there are enhancer targets:
    if(length(motif_enhancer_target_genes) > 0){
      enhancer_edge_df <- data.frame(
        from = motif_name,
        to = as.character(motif_enhancer_target_genes),
        type = 'enhancer'
      )
    } else{enhancer_edge_df <- data.frame()}

    target_genes <- rbind(promoter_edge_df, enhancer_edge_df)
    write.xlsx(target_genes, file = paste0("7.TF.analysis/targets/", motif_name, "_targetGenes.xlsx"), overwrite = T)

    #### Compute module score for these target genes
    gene_list <- list(unique(target_genes$to))
    names(gene_list) <- paste0(motif_name, '_targets')

    # add in whole data
    scRNA.data <- AddModuleScore(
        scRNA.data,
        features=gene_list,
        pool = rownames(scRNA.data),
        name=paste0(motif_name, '_targets')
    )

    #### plot module score feature plot:
    # plot promoter target gene module score for this TF:
    p1 <- FeaturePlot(scRNA.data, features = paste0(motif_name, '_targets1'), order = order_values, reduction = reduction) +
    scale_color_gradient2(low=scales::muted('blue'), mid='white', high=scales::muted('red')) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "in"))  +
    ggtitle(paste0(motif_name, ' target score'))
    pdf(paste0('7.TF.analysis/targets/', celltype, '_', motif_name, '_targets.pdf'), width=5, height=4, useDingbats=FALSE)
    print(p1)

    # plot chromVAR deviation for this TF:
    DefaultAssay(scATAC.data) <- "chromvar"
    p2 <- FeaturePlot(scATAC.data, features = gsub("_", "-", motif_ID), order = order_values, reduction = reduction) +
    scale_color_gradient2(
        low=rgb(32, 67, 37, maxColorValue=255), mid='white', high=rgb(58, 22, 72, maxColorValue=255)) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "in")) +
        ggtitle(paste0(motif_name, ' motif'))
    print(p2)

    #### cluster violin plot for target expression modules
    signature.scores <- scRNA.data@meta.data
    signature.scores$cellType_low <- factor(signature.scores$cellType_low)
    p3 <- ggviolin(signature.scores, y=paste0(motif_name, '_targets1'), x='cellType_low', pt.size=0, fill = 'cellType_low', add = "mean_sd", error.plot = "crossbar") +
    xlab('') + ylab(paste0(motif_name, ' target score')) + ggtitle('') +
    theme(plot.margin = unit(c(0, 0, 0, 0.1), "in"), axis.title.y=element_text(face='bold'), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    NoLegend() + stat_compare_means(ref.group = "Tumor", label = "p.signif")
    print(p3)

    p4 <- VlnPlot(scATAC.data, assay='chromvar', features=gsub("_", "-", motif_ID), group.by='AnnotatedcellType', pt.size=0) +
    stat_compare_means(label='p.signif', ref.group = "Tumor", label.y = 12)  +
    xlab('') + ylab(paste0(motif_name, ' deviation')) + ggtitle('') +
    theme(plot.margin = unit(c(0, 0, 0, 0.1), "in"), axis.title.y=element_text(face='bold'), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    NoLegend() + stat_compare_means(ref.group = "Tumor", label = "p.signif")
    print(p4)

    #### motif logo
    PPM <- PWMatrixToProbMatrix(human_pwms_v2[[motif_ID]])
    p5 <- ggseqlogo(PPM) + ggtitle(motif_name) + theme(plot.title = element_text(hjust = 0.5))
    print(p5)    
    dev.off()
})

################################################################################
# Construct TF nets
################################################################################
####---- 1.build the vertex and edge

library(igraph)
library(RColorBrewer)

interest.TFs <- c("OTP", "ISL1", "VENTX", "HOXC5")
TFs.info <- motif.info[match(interest.TFs, motif.info$TF),]

# for each motif, find target genes:
use_variable_genes <- T # require the variable gene
edge_df.list <- data.frame()
vertex_df.list <- data.frame()
for(i in 1:nrow(TFs.info)){
  motif_name <- TFs.info[i,2]
  motif_ID <- TFs.info[i,1]

  # get list of promoter and enhancer targets of these TFs
  DefaultAssay(scATAC.sub) <- "Peaks"
  motif_accessible <- Motifs(scATAC.sub)@data[,motif_ID]
  motif_accessible <- names(motif_accessible)[motif_accessible > 0]
  motif_accessible <- motif_accessible[motif_accessible %in% VariableFeatures(scATAC.sub)]  
  motif_accessible_promoters <- motif_accessible[motif_accessible %in% peak.info$peaks[which(peak.info$originType == "Promoter (<=1kb)")]]
  motif_promoter_target_genes <- peak.info$SYMBOL[match(motif_accessible_promoters, peak.info$peaks)]
  motif_promoter_target_genes <- as.character(motif_promoter_target_genes)

  # variable genes only?
  if(use_variable_genes){
    motif_promoter_target_genes <- motif_promoter_target_genes[motif_promoter_target_genes %in% VariableFeatures(scRNA.sub)] %>% as.character %>% unique
  } else{
    motif_promoter_target_genes<- motif_promoter_target_genes %>% as.character %>% unique
  }

  # enhancer target genes
  motif_accessible_enhancers <- motif_accessible[motif_accessible %in% cCREs$Peak2]
  motif_cCREs <- subset(cCREs, Peak2 %in% motif_accessible_enhancers)
  motif_enhancer_target_genes <- motif_cCREs$target_gene

  # variable genes only?
  if(use_variable_genes){
    motif_enhancer_target_genes <- motif_enhancer_target_genes[motif_enhancer_target_genes %in% VariableFeatures(scRNA.sub)] %>% as.character %>% unique
  }else{
    motif_enhancer_target_genes<- motif_enhancer_target_genes %>% as.character %>% unique
  }

  vertex_df <- data.frame(name = c(motif_name, as.character(unique(c(motif_enhancer_target_genes, motif_promoter_target_genes)))))

  # check if there are promoter targets:
  if(length(motif_promoter_target_genes) > 0){
    promoter_edge_df <- data.frame(
      from = motif_name,
      to = as.character(motif_promoter_target_genes),
      type = 'promoter'
    )
  } else{promoter_edge_df <- data.frame()}

  # check if there are enhancer targets:
  if(length(motif_enhancer_target_genes) > 0){
    enhancer_edge_df <- data.frame(
      from = motif_name,
      to = as.character(motif_enhancer_target_genes),
      type = 'enhancer'
    )
  } else{enhancer_edge_df <- data.frame()}

  edge_df <- rbind(promoter_edge_df, enhancer_edge_df)

  edge_df.list <- rbind(edge_df.list, edge_df)
  vertex_df.list <- rbind(vertex_df.list, vertex_df)
}

vertex_df.list <- data.frame(name=na.omit(as.character(unique(vertex_df.list$name))))
vertex_df.list$name <- as.character(vertex_df.list$name)
edge_df.list <- na.omit(edge_df.list)

####---- 2.visual settings for network
cellType.DEGs <- scRNA.all.DEGs[which(scRNA.all.DEGs$cluster == celltype),]

# color vertices based on Diagnosis DEGs:
up_genes <- cellType.DEGs %>% filter(p_val_adj < 0.05 & avg_log2FC > 1) %>% .$gene
#down_genes <- cellType.DEGs %>% filter(p_val_adj < 0.05 & avg_log2FC < -1) %>% .$gene # not include these genes

# remove labels if gene is not DE, or not a TF:
de_targets <- as.character(vertex_df.list$name[vertex_df.list$name %in% unique(c(up_genes))])
vertex_df.list$label <- ifelse(vertex_df.list$name %in% de_targets, vertex_df.list$name, '')
vertex_df.list$label <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), vertex_df.list$name, vertex_df.list$label)

# set node color based on DEGs:
vertex_df.list$color <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), '#1E90FF', "#FFFFFF")
vertex_df.list$color <- ifelse(vertex_df.list$name %in% up_genes, "#E87D72",  vertex_df.list$color)
#vertex_df.list$color <- ifelse(vertex_df.list$name %in% down_genes, '#7BD0D3', vertex_df.list$color)
vertex_df.list$color <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), '#1E90FF',vertex_df.list$color)

# italics font for genes:
vertex_df.list$font <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), 2, 1)

# set size to larger if the gene is a TF:
vertex_df.list$size <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), 10, 2)

other_tfs <- setdiff(motif.info$TF, interest.TFs)

vertex_df.list$size <- ifelse((vertex_df.list$name %in% de_targets | vertex_df.list$name %in% other_tfs), 5, 2)
vertex_df.list$size <- ifelse(vertex_df.list$name %in% interest.TFs, 10, vertex_df.list$size)

# graph all nodes
enhancer_color <- '#FFC125'
promoter_color <- '#00CED1'

g <- igraph::graph_from_data_frame(edge_df.list, directed=TRUE, vertices = vertex_df.list)
l <- layout_nicely(g)

edge_colors <- ifelse(E(g)$type == 'promoter', promoter_color, enhancer_color)

pdf(paste0('7.TF.analysis/network/', celltype, '_TF_interaction_graph.pdf'), width=10, height=10, useDingbats=FALSE)
plot(
  g, layout=l,
  vertex.size=vertex_df.list$size,
  edge.color=edge_colors,
  edge.alpha=0.5,
  vertex.color=vertex_df.list$color,
  vertex.label=vertex_df.list$label, vertex.label.family='Helvetica', vertex.label.font=vertex_df.list$font,
  vertex.label.color = 'black',
  vertex.frame.color = "grey",
  vertex.alpha = 0.8,
  edge.arrow.size=0.25
)
dev.off()

###########################Comparative analysis of the target genes of these 4 TFs###################################
setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/7.TF.analysis/targets")
library(VennDiagram)
# 2022.1.1
TFs <- c("HOXC5", "OTP", "ISL1", "VENTX")
#### 1. target gene overlap
target.genes.list <- lapply(TFs, function(x){
  targets <- read.xlsx(paste0("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/7.TF.analysis/targets/", x, "_All.targetGenes.xlsx"))
  targets <- unique(targets$to)
  return(targets)
})
names(target.genes.list) <- TFs
venn.plot <- venn.diagram(
  x = list(
    HOXC5 = target.genes.list$HOXC5,
    OTP = target.genes.list$OTP,
    ISL1 = target.genes.list$ISL1,
    VENTX = target.genes.list$VENTX
  ),
  filename = "four.targets.png",
  col = "black",
  lty = "dotted", #边框线型改为"dotted"虚线
  lwd = 3, # 边框线的宽度
  fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
  alpha = 0.50,
  label.col = c("orange", "white", "darkorchid4", "white", "white", "white",
                "white", "white", "darkblue", "white",
                "white", "white", "white", "darkgreen", "white"),
  cex = 2.0,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
  cat.cex = 1.8,
  cat.fontface = "bold",
  cat.fontfamily = "serif"
)

#### enhancer
target.genes.list <- lapply(TFs, function(x){
  targets <- read.xlsx(paste0("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/7.TF.analysis/targets/", x, "_All.targetGenes.xlsx"))
  targets <- targets[which(targets$type == "enhancer"),]
  targets <- unique(targets$to)
  return(targets)
})
names(target.genes.list) <- TFs
venn.plot <- venn.diagram(
  x = list(
    HOXC5 = target.genes.list$HOXC5,
    OTP = target.genes.list$OTP,
    ISL1 = target.genes.list$ISL1,
    VENTX = target.genes.list$VENTX
  ),
  filename = "four.targets.enhancer.png",
  col = "black",
  lty = "dotted", #边框线型改为"dotted"虚线
  lwd = 3, # 边框线的宽度
  fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
  alpha = 0.50,
  label.col = c("orange", "white", "darkorchid4", "white", "white", "white",
                "white", "white", "darkblue", "white",
                "white", "white", "white", "darkgreen", "white"),
  cex = 2.0,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
  cat.cex = 1.8,
  cat.fontface = "bold",
  cat.fontfamily = "serif"
)

#### pathway enrichment of target genes
source("/home/longzhilin/Analysis_Code/PathwayEnrichment/clusterProfiler.enricher.R")
TFs <- c("HOXC5", "OTP", "ISL1", "VENTX")
pdf("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/7.TF.analysis/pathwayEnrichment/pathwayEnrichment.pdf")
enrich.result <- lapply(TFs, function(x){
  cat("Analysis ", x, "\n")
  targets <- read.xlsx(paste0("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/7.TF.analysis/targets/", x, "_targetGenes.xlsx"))
  cat("target gene ", length(unique(targets$to)), "\n")
  res <- cluterProfiler.enricher(gene = unique(targets$to), geneType = "SYMBOL", db.type = "MsigDB", saveDir = "/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/7.TF.analysis/pathwayEnrichment",
                          title = x, qvalueCutoff = 0.2, pvalueCutoff = 0.2, showCategory = 10)
  return(res)
})
dev.off()

####cluster profiler enrichment result display
library(ggthemes)
#top10
HOXC5.enrich <- read.csv("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/7.TF.analysis/pathwayEnrichment/HOXC5_enricherResult.csv", header = T)
VENTX.enrich <- read.csv("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/7.TF.analysis/pathwayEnrichment/VENTX_enricherResult.csv", header = T)
ISL1.enrich <- read.csv("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/7.TF.analysis/pathwayEnrichment/ISL1_enricherResult.csv", header = T)
OTP.enrich <- read.csv("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/7.TF.analysis/pathwayEnrichment/OTP_enricherResult.csv", header = T)
enrich.list <- list(HOXC5=HOXC5.enrich, VENTX=VENTX.enrich, ISL1=ISL1.enrich, OTP=OTP.enrich)
enrich <- list(HOXC5 = HOXC5.enrich[c(5,13,15,26,30,47:49,91,108,132), c("ID")],
               VENTX = VENTX.enrich[c(1:3,5:6,12,14:16,37), c("ID")],
               ISL1 = ISL1.enrich[c(1,12:16,18,34,113), c("ID")],
               OTP = OTP.enrich[c(1,3:4,6:8,10,15,20:21,40), c("ID")])
enrich.path <- Reduce(function(x,y) c(x,y), enrich)
enrich.path <- unique(enrich.path)
enrich.path <- enrich.path[-c(1,2,3,5,7,8,9,10,12,13,14,15,18,22,23,25,26,27,28,29)]
enrich.list.select <- lapply(enrich.list, function(x){
  idx <- which(x$ID %in% enrich.path)
  return(x[idx, c("ID", "Count", "pvalue", "p.adjust")])
})
enrich.res <- Reduce(function(x,y) rbind(x,y), enrich.list.select)
a <- sapply(enrich.list.select, function(x){nrow(x)})
enrich.res$Type <- c(rep("HOXC5", as.numeric(a[1])), rep("VENTX", as.numeric(a[2])), rep("ISL1", as.numeric(a[3])), rep("OTP", as.numeric(a[4])))
enrich.res$ID <- factor(enrich.res$ID, levels = unique(enrich.res$ID))
enrich.res$pvalue <- -log10(enrich.res$pvalue)
#处理名字
enrich.res$ID <- gsub("_", " ", enrich.res$ID)
paths <- unique(enrich.res$ID)
enrich.res$ID <- factor(enrich.res$ID, levels = paths[c(10,9,6,7,8,5,2,3,4,1)])

#source(file = "/home/longzhilin/Analysis_Code/Visualization/ggplot.bubble.plot.R")
pdf("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/7.TF.analysis/pathwayEnrichment/ggplot.bubble.pathway.pdf", width = 6,height = 4)
p1 <- ggplot(enrich.res, 
            aes(x = Type, 
                y = ID, 
                size = Count,
                fill = pvalue))
p2 <- p1 + geom_point(shape = 21, alpha = 0.7) + theme_few()
p2 <- p2 + scale_fill_gradientn(colours = c("#ffe11a", "#f9ac1d", "#e92925"), breaks = c(1.82, 3, 6), labels = c(0.015, 0.001, 0.000001)) + scale_size_continuous(range = c(3,7), breaks = c(4,6,8,12), labels = c(4,6,8,12))
p2 <- p2 + theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1, size = 8), axis.text.y = element_text(size = 8)) + xlab("") + ylab("")
print(p2)
dev.off()

## FXYD2和CRYAB
plot.cellType <- c("CD4+ T cell", "Treg", "CD8+ T cell", "NK/NKT cell", "Proliferative CD8+ T cell", "B cell", "Macrophage", "Monocyte", "Dendritic cell", "Neutrophil", "Mast cell", "Endothelium (VCAM1+)", "Endothelium (VCAM1-)", "Mesangial cell", "Tumor")
scRNA.data$cellType_low <- factor(scRNA.data$cellType_low, levels = plot.cellType)
pdf("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/7.TF.analysis/FXYD2_CRYAB.pdf")
VlnPlot(scRNA.data, group.by = "cellType_low", features = c("FXYD2", "CRYAB"), ncol = 1, pt.size = 0) & xlab("")
dev.off()

## replot the target score
tumor.TFs <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/5.Motif/Analysis/tumor.specific.TFs.rds")
four.TFs <- c("HOXC5", "VENTX", "ISL1", "OTP")

target.genes.list <- lapply(four.TFs, function(x){
  targets <- read.xlsx(paste0("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/7.TF.analysis/targets/", x, "_All.targetGenes.xlsx"))
  targets <- unique(targets$to)
  return(targets)
})
names(target.genes.list) <- four.TFs
target.genes.list$AllTFs <- tumor.TFs$Name
TF.names <- names(target.genes.list)

for(i in TF.names){
  scRNA.data <- AddModuleScore(scRNA.data, features = target.genes.list[i], name = i, pool = rownames(scRNA.data))
}

#scRNA.data <- AddModuleScore(scRNA.data, features = target.genes.list, name = TF.names, pool = rownames(scRNA.data))
colnames(scRNA.data@meta.data)[(ncol(scRNA.data@meta.data)-4):ncol(scRNA.data@meta.data)] <- paste0(TF.names, "_targets")
scRNA.data$cellType_low <- factor(scRNA.data$cellType_low, levels = c(setdiff(names(table(scRNA.data$cellType_low)), "Tumor"), "Tumor"))

signature.scores <- scRNA.data@meta.data
pdf("7.TF.analysis/targets/targetScores.pdf", height = 4)
p <- lapply(TF.names, function(x){
    p <- ggviolin(signature.scores, y=paste0(x, '_targets'), x='cellType_low', pt.size=0, fill = 'cellType_low', add = "mean_sd", error.plot = "crossbar") +
    xlab('') + ylab(paste0(x, ' target score')) + ggtitle('') +
    theme(plot.margin = unit(c(0, 0, 0, 0.1), "in"), axis.title.y=element_text(face='bold'), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    NoLegend() + stat_compare_means(ref.group = "Tumor", label = "p.signif")
    print(p)
    return(p)
})
dev.off()