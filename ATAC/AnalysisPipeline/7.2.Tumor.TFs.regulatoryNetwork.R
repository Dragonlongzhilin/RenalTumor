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
scRNA.DEGs <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/2.Cluster/AnnotateCellType/cellType.sig.pos.DEGs.rds")

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
scATAC.sub <- subset(scATAC.data, cellType_low == celltype)
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

################################################################################
# Find promoters & genes with accessible TF binding sites
################################################################################
# settings for featureplot
order_values <- TRUE
reduction <- 'umap'

res <- sapply(Tumor.TF.motifs$Name, function(x){
    motif_name <- x
    motif_ID <- motif.info$originName[match(x, motif.info$TF)] # e.g., ENSG00000008196_LINE2_TFAP2B_D_N1

    # get all regions with motif binding site:
    DefaultAssay(scATAC.sub) <- "Peaks"
    motif_accessible <- Motifs(scATAC.sub)@data[, motif_ID]
    motif_accessible <- names(motif_accessible)[motif_accessible > 0]

    # subset this list by top features
    motif_accessible <- motif_accessible[motif_accessible %in% VariableFeatures(scATAC.sub)]

    # which of these peaks are at promoters?
    motif_accessible_promoters <- motif_accessible[motif_accessible %in% peak.info$peaks[which(peak.info$originType == "Promoter (<=1kb)")]]
    # which genes are associated with these promoters?
    motif_target_genes <- peak.info$SYMBOL[match(motif_accessible_promoters, peak.info$peaks)]

    # optional:
    # which of these genes are highly expressed in snRNA-seq?
    motif_target_genes.RNA <- motif_target_genes[motif_target_genes %in% VariableFeatures(scRNA.sub)] %>% as.character %>% unique

    # remove genes that are not in the seurat obj
    motif_target_genes <- as.character(motif_target_genes)
    motif_target_genes <- motif_target_genes[motif_target_genes %in% rownames(scRNA.sub)]

    # output 
    target.genelists <- list(data.frame(motif_target_genes.RNA), data.frame(motif_target_genes))
    write.xlsx(target.genelists, file = paste0("7.TF.analysis/targets/", x, "_targetGenes.xlsx"), sheetName = c("overlap with hvgs in scRNA-seq", "all target genes"))
    
    ################################################################################
    # Compute module score for these target genes
    ################################################################################
    gene_list <- list(motif_target_genes)
    names(gene_list) <- paste0(motif_name, '_targets')

    # add in whole data
    scRNA.data <- AddModuleScore(
        scRNA.data,
        features=gene_list,
        pool = rownames(scRNA.data),
        name=paste0(motif_name, '_targets')
    )

    ################################################################################
    # plot module score feature plot:
    ################################################################################

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

    ################################################################################
    # cluster violin plot for target expression modules
    ################################################################################
    p3 <- VlnPlot(scRNA.data, features=paste0(motif_name, '_targets1'), group.by='cellType_low', pt.size=0) +
    stat_compare_means(method='wilcox.test', label='p.signif', label.y=0.1) +
    geom_hline(yintercept = 0, linetype='dashed') +
    xlab('') + ylab(paste0(motif_name, ' targets')) + ggtitle('') +
    theme(plot.margin = unit(c(0, 0, 0, 0.1), "in"), axis.title.y=element_text(face='bold')) +
    NoLegend()
    print(p3)

    p4 <- VlnPlot(scATAC.data, assay='chromvar', features=gsub("_", "-", motif_ID), group.by='cellType_low', pt.size=0) +
    stat_compare_means(method='wilcox.test', label='p.signif', label.y=5)  +
    xlab('') + ylab(paste0(motif_name, ' deviation')) + ggtitle('') +
    theme(plot.margin = unit(c(0, 0, 0, 0.1), "in"), axis.title.y=element_text(face='bold')) +
    NoLegend() + geom_hline(yintercept = 0, linetype='dashed')
    print(p4)

    # motif logo
    PPM <- PWMatrixToProbMatrix(human_pwms_v2[[motif_ID]])
    p5 <- ggseqlogo(PPM) + ggtitle(motif_name) + theme(plot.title = element_text(hjust = 0.5))
    print(p5)    
    dev.off()
})

##enrichR 功能富集分析
source(file = "/home/longzhilin/Analysis_Code/PathwayEnrichment/pathwayEnrichment.Enrichr.R")
interest.TFs <- c("OTP", "ISL1", "VENTX", "HOXC5")
targetGenes <- lapply(interest.TFs, function(x){
  genes <- read.xlsx(paste0("7.TF.analysis/targets/", x, "_targetGenes.xlsx"), sheet = 2)
  return(as.character(genes[,1]))
})
names(targetGenes) <- interest.TFs
dbs <- c("MSigDB_Hallmark_2020", "KEGG_2019_Human", "GO_Biological_Process_2018", "BioCarta_2016", "Reactome_2016")
all.enrichr.res <- pathwayEnrichment.Enrichr(geneList = targetGenes, dbs = dbs, plot = F, orderBy = "P.value")
combined_output <- data.frame()
for(db in dbs){
  all.enrichr.res[[db]]$db <- db
  all.enrichr.res[[db]]$NMF_module <- cur_mod
  all.enrichr.res[[db]] <-all.enrichr.res[[db]] %>% subset(P.value <= 0.01)
  print(dim(all.enrichr.res[[db]]))
  combined_output <- rbind(combined_output, cur_results[[db]])
}

################################################################################
# Construct TF nets
################################################################################
library(igraph)
library(RColorBrewer)
cCREs <- readRDS("6.Co-Accessible/Tumor_peak_gene_correlation.rds") # 4733 cCREs
cCREs <- cCREs %>% subset(distance_bp/1000>=10 & peak2_type != 'Promoter')
cCREs$target_gene <- cCREs$peak1_nearestGene
cCREs$peak_gene <- paste0(cCREs$Peak2, '_', cCREs$peak1_nearestGene)

# get links for this cell type
#cCREs <- subset(top_link_df, celltype %in% cur_celltype)
cCREs$Peak1 <- gsub("_", "-", cCREs$Peak1)
cCREs$Peak2 <- gsub("_", "-", cCREs$Peak2)

interest.TFs <- c("OTP", "ISL1", "VENTX", "HOXC5")
#interest.TFs <- c("HNF1A", "HNF1B")
TFs.info <- motif.info[match(interest.TFs, motif.info$TF),]

# for each motif, find genes:
use_variable_genes <- T # require the variable gene
motif_list <- list()
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
  motif_target_genes <- peak.info$SYMBOL[match(motif_accessible_promoters, peak.info$peaks)]
  motif_target_genes <- as.character(motif_target_genes)

  # variable genes only?
  if(use_variable_genes){
    motif_target_genes <- motif_target_genes[motif_target_genes %in% VariableFeatures(scRNA.sub)] %>% as.character %>% unique
  } else{
    motif_target_genes<- motif_target_genes %>% as.character %>% unique
  }

  # enhancer target genes
  motif_accessible_enhancers <- motif_accessible[motif_accessible %in% cCREs$Peak2]
  motif_cCREs <- subset(cCREs, Peak2 %in% motif_accessible_enhancers)
  motif_enhancers_target_genes <- motif_cCREs$target_gene

  # variable genes only?
  if(use_variable_genes){
    motif_enhancers_target_genes <- motif_enhancers_target_genes[motif_enhancers_target_genes %in% VariableFeatures(scRNA.sub)] %>% as.character %>% unique
  }else{
    motif_enhancers_target_genes<- motif_enhancers_target_genes %>% as.character %>% unique
  }

  vertex_df <- data.frame(name = c(motif_name, as.character(unique(c(motif_enhancers_target_genes, motif_target_genes)))))

  # check if there are promoter targets:
  if(length(motif_target_genes) > 0){
    promoter_edge_df <- data.frame(
      from = motif_name,
      to = as.character(motif_target_genes),
      type = 'promoter'
    )
  } else{promoter_edge_df <- data.frame()}

  # check if there are enhancer targets:
  if(length(motif_enhancers_target_genes) > 0){
    enhancer_edge_df <- data.frame(
      from = motif_name,
      to = as.character(motif_enhancers_target_genes),
      type = 'enhancer'
    )
  } else{enhancer_edge_df <- data.frame()}

  #cur_edge_df <- rbind(cur_promoter_edge_df, cur_enhancer_edge_df, cur_repressors_edge_df)
  edge_df <- rbind(promoter_edge_df, enhancer_edge_df)

  edge_df.list <- rbind(edge_df.list, edge_df)
  vertex_df.list <- rbind(vertex_df.list, vertex_df)
}

vertex_df.list <- data.frame(name=na.omit(as.character(unique(vertex_df.list$name))))
vertex_df.list$name <- as.character(vertex_df.list$name)
edge_df.list <- na.omit(edge_df.list)

################################################################################
# visual settings for network
################################################################################
cellType.DEGs <- scRNA.DEGs[which(scRNA.DEGs$cluster == celltype),]

# color vertices based on Diagnosis DEGs:
up_genes <- cellType.DEGs %>% filter(p_val_adj < 0.05 & avg_log2FC >= 1) %>% .$gene
down_genes <- cellType.DEGs %>% filter(p_val_adj < 0.05 & avg_log2FC < -1) %>% .$gene

# remove labels if gene is not DE, or not a TF:
de_targets <- as.character(vertex_df.list$name[vertex_df.list$name %in% unique(c(up_genes, down_genes))])
vertex_df.list$label <- ifelse(vertex_df.list$name %in% de_targets, vertex_df.list$name, '')
vertex_df.list$label <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), vertex_df.list$name, vertex_df.list$label)
#vertex_df.list$label <- ifelse(vertex_df.list$name %in% gwas_genes, vertex_df.list$name, vertex_df.list$label)

# set node color based on control vs AD DEGs:
vertex_df.list$color <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), '#1E90FF', "#FFFFFF")
vertex_df.list$color <- ifelse(vertex_df.list$name %in% up_genes, "#E87D72",  "#FFFFFF")
vertex_df.list$color <- ifelse(vertex_df.list$name %in% down_genes, '#55BCC2', vertex_df.list$color)
vertex_df.list$color <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), '#1E90FF',vertex_df.list$color)
#vertex_df.list$color <- ifelse(vertex_df.list$name %in% gwas_genes, 'orange', vertex_df.list$color)

# italics font for genes:
vertex_df.list$font <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), 2, 1)

# set size to larger if the gene is a TF:
vertex_df.list$size <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), 10, 2)

other_tfs <- setdiff(motif.info$TF, interest.TFs)

#vertex_df.list$size <- ifelse((vertex_df.list$name %in% de_targets | vertex_df.list$name %in% gwas_genes | vertex_df.list$name %in% other_tfs), 5, 2)
vertex_df.list$size <- ifelse((vertex_df.list$name %in% de_targets | vertex_df.list$name %in% other_tfs), 5, 2)
vertex_df.list$size <- ifelse(vertex_df.list$name %in% interest.TFs, 10, vertex_df.list$size)

################################################################################
# graph all nodes
################################################################################
enhancer_color <- '#FFC125'
promoter_color <- '#00CED1'

# repressor_color <- 'gray'

g <- igraph::graph_from_data_frame(edge_df.list, directed=TRUE, vertices = vertex_df.list)
l <- layout_nicely(g)

edge_colors <- ifelse(E(g)$type == 'promoter', promoter_color, enhancer_color)
# edge_colors <- ifelse(E(g)$type == 'repressor', repressor_color, edge_colors)

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
