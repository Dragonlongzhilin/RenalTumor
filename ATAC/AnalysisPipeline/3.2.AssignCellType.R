#' @description cell type annotation

library(Signac)
library(Seurat)
library(harmony)
library(SeuratWrappers)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86) #---GRCh38 (hg38)
library(ggplot2)
library(patchwork)
library(clustree)
set.seed(101)
library(GenomicRanges)
library(tibble)
library(dplyr)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) #50G

source(file = "/home/longzhilin/Analysis_Code/Plot_colorPaletters.R")
setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC")
source("/home/longzhilin/Analysis_Code/SingleCell/scATAC.Integrate.multipleSample.R")
scATAC.data <- readRDS("scATAC.sub.rds")

#### recluster
DefaultAssay(scATAC.data) <- "ATAC"
scATAC.data <- RunTFIDF(scATAC.data)
scATAC.data <- FindTopFeatures(scATAC.data)
pdf("2.Cluster/cluster.pdf")
scATAC.data <- Harmony.integration.reduceDimension(scATAC.object = scATAC.data, set.resolutions = seq(0.2, 1.2, by = 0.1), groups = "dataset", assay = "ATAC", PC = 15)
ElbowPlot(object = scATAC.data, ndims = 30, reduction = "lsi")
dev.off()
scATAC.data$seurat_clusters <- scATAC.data$ATAC_snn_res.0.4

#### annotated cluster
#marker
DefaultAssay(scATAC.data) <- "ACTIVITY"
cell.type.markers <- read.table(file = "/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/2.Cluster/AnnotateCellType/cellMarker.txt", header = T, stringsAsFactors = F, sep = "\t")
pdf("2.Cluster/cluster.marker.dotplot.pdf")
DotPlot(scATAC.data, cols = c("#1e90ff", "#ff5a36"), features = cell.type.markers$Gene, group.by = "seurat_clusters", dot.scale = 5) + theme(axis.text.x = element_text(size = 6, angle = 90)) + theme(axis.text.y = element_text(size = 6))
DotPlot(scATAC.data, cols = c("#1e90ff", "#ff5a36"), features = cell.type.markers$Gene, group.by = "seurat_clusters", dot.scale = 5) + NoLegend()+ theme(axis.text.x = element_text(size = 6, angle = 90)) + theme(axis.text.y = element_text(size = 6))
ElbowPlot(object = scATAC.data, ndims = 30, reduction = "lsi")
DimPlot(object = scATAC.data, reduction = 'umap',label = TRUE, group.by = "cellType_low") + NoLegend()
DimPlot(object = scATAC.data, reduction = 'tsne',label = TRUE, group.by = "cellType_low") + NoLegend()
dev.off()

cluster.label <- scATAC.data$seurat_clusters
cluster.label <- gsub("^0$", "CD8+ T cell", cluster.label)
cluster.label <- gsub("^1$", "CD4+ T cell", cluster.label)
cluster.label <- gsub("^2$", "NK/NKT cell", cluster.label)
cluster.label <- gsub("^3$", "Macrophage", cluster.label)
cluster.label <- gsub("^4$", "Endothelium (VCAM1-)", cluster.label)
cluster.label <- gsub("^5$", "NK/NKT cell", cluster.label)
cluster.label <- gsub("^6$", "Mesangial cell", cluster.label)
cluster.label <- gsub("^7$", "Tumor", cluster.label)
cluster.label <- gsub("^8$", "Monocyte", cluster.label)
cluster.label <- gsub("^9$", "B cell", cluster.label)
cluster.label <- gsub("^10$", "Endothelium (VCAM1+)", cluster.label)
cluster.label <- gsub("^11$", "Treg", cluster.label)
cluster.label <- gsub("^12$", "Mast cell", cluster.label)
scATAC.data <- AddMetaData(scATAC.data, cluster.label, col.name = "AnnotatedcellType")
saveRDS(scATAC.data, file = "scATAC.data.rds")

## map the umap color between scRNA and scATAC
emb <- scRNA.data@reductions$umap@cell.embeddings
gg <- UMAPPlot(scRNA.data, group.by = "cellType_low")
color.info <- data.frame(colors = ggplot_build(gg)$data[[1]])
Idents(scRNA.data) <- scRNA.data$cellType_low
color.idx <- data.frame(cellType = names(sort(table(scRNA.data$cellType_low))), color = names(sort(table(color.info$colors.colour))))

scATAC.colors <- color.idx[unique(scATAC.data$AnnotatedcellType),]
cell.number <- as.data.frame(sort(table(scATAC.data$AnnotatedcellType), decreasing = T))
idx <- match(cell.number$Var1, color.idx$cellType)
cell.number$colors <- color.idx$color[idx]
colors <- cell.number$colors
names(colors) <- cell.number$Var1
saveRDS(colors, file  = "2.Cluster/AnnotateCellType/cellType.color.rds")
pdf("2.Cluster/AnnotateCellType/cellType.pdf")
DimPlot(object = scATAC.data, cols = colors, reduction = 'tsne',label = TRUE, group.by = "AnnotatedcellType")+NoLegend()
DimPlot(object = scATAC.data, cols = colors, reduction = 'tsne',label = TRUE, group.by = "AnnotatedcellType")
DimPlot(object = scATAC.data, cols = colors, reduction = 'umap',label = TRUE, group.by = "AnnotatedcellType")+NoLegend()
DimPlot(object = scATAC.data, cols = colors, reduction = 'umap',label = TRUE, group.by = "AnnotatedcellType")
dev.off()

####comparative cell type
#define Jaccard Similarity function
jaccard <- function(a, b) {
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return (intersection/union)
}
annotated.info <- scATAC.data@meta.data

scRNA.cellType <- sort(as.character(unique(annotated.info$cellType_low)))
scATAC.cellType <- sort(as.character(unique(annotated.info$AnnotatedcellType)))
scATAC.cluster <- sort(unique(annotated.info$seurat_clusters))

jaccard.matrix.cellType <- sapply(scRNA.cellType, function(x){
    scRNA.idx <- which(annotated.info$cellType_low == x)
    
    partial.matrix <- sapply(scATAC.cellType, function(y){
        scATAC.idx <- which(annotated.info$AnnotatedcellType == y)
        return(jaccard(scRNA.idx, scATAC.idx))
    })
    return(partial.matrix)
})

jaccard.matrix.cluster <- sapply(scRNA.cellType, function(x){
    scRNA.idx <- which(annotated.info$cellType_low == x)
    
    partial.matrix <- sapply(scATAC.cluster, function(y){
        scATAC.idx <- which(annotated.info$seurat_clusters == y)
        return(jaccard(scRNA.idx, scATAC.idx))
    })
    return(partial.matrix)
})
rownames(jaccard.matrix.cluster) <- paste0("Cluster ", scATAC.cluster)
pdf("2.Cluster/AnnotateCellType/jaccard.similarity.pdf")
Heatmap(jaccard.matrix.cellType, name = "Jaccard similarity", rect_gp = gpar(col = "lightgray", lwd = 1), border = "lightgray", cluster_rows = F, cluster_columns = F, col = c("white", "red"), show_row_dend = F, show_column_dend = F, width = unit(10, "cm"), height = unit(8, "cm"))
Heatmap(jaccard.matrix.cluster, name = "Jaccard similarity", rect_gp = gpar(col = "lightgray", lwd = 1), border = "lightgray", cluster_rows = F, cluster_columns = F, col = c("white", "red"), show_row_dend = F, show_column_dend = F, width = unit(10, "cm"), height = unit(10, "cm"))
dev.off()

####--- plot the number ratio of sample in each cell cluster
source(file = "/home/longzhilin/Analysis_Code/code/ratio.plot.R")
pdf("2.Cluster/cluster.number.ratio.pro.pdf", height = 4, width = 7)
ratio.plot(seurat.object = scATAC.data, id.vars1 = "orig.ident", id.vars2 = "seurat_clusters")
dev.off()

pdf("2.Cluster/AnnotateCellType/cellType.number.ratio.pdf", height = 4, width = 7)
ratio.plot(seurat.object = scATAC.data, id.vars1 = "orig.ident", id.vars2 = "AnnotatedcellType")
dev.off()

# cell lineage ratio
cellType.ratio <- as.data.frame(table(scATAC.data$AnnotatedcellType))
cellType.ratio$Type <- rep("Lymphoid cells", nrow(cellType.ratio))
idx <- which(cellType.ratio$Var1 %in% c("Macrophage", "Monocyte", "Mast cell"))
cellType.ratio$Type[idx] <- "Myeloid cells"
idx <- which(cellType.ratio$Var1 %in% c("Endothelium (VCAM1+)", "Endothelium (VCAM1-)", "Mesangial cell"))
cellType.ratio$Type[idx] <- "Other cells"
idx <- which(cellType.ratio$Var1 == "Tumor")
cellType.ratio$Type[idx] <- "Tumor cells"
group.ratio <- tapply(cellType.ratio$Freq, cellType.ratio$Type, sum)
group.ratio <- data.frame(Type = names(group.ratio), Ratio = group.ratio)
labs <- paste0(group.ratio$Type, " (", round(group.ratio$Ratio/sum(group.ratio$Ratio), 4)*100, "%)")
pdf("2.Cluster/AnnotateCellType/lineage.ratio.pdf")
p <- ggpie(group.ratio, "Ratio", label = labs,
      fill = "Type", color = "white", lab.pos = "in",
      palette = "npgs")
print(p)
dev.off()

####---- plot marker expression
plot.cellType <- c("CD4+ T cell", "Treg", "CD8+ T cell", "NK/NKT cell", "B cell", "Macrophage", "Monocyte", "Mast cell", "Endothelium (VCAM1+)", "Endothelium (VCAM1-)", "Mesangial cell", "Tumor")
scATAC.data$AnnotatedcellType <- factor(scATAC.data$AnnotatedcellType, levels = plot.cellType)
DefaultAssay(scATAC.data) <- "ACTIVITY"
cell.type.markers <- read.table(file = "/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/2.Cluster/AnnotateCellType/cellMarker.txt", header = T, stringsAsFactors = F, sep = "\t")
pdf("2.Cluster/AnnotateCellType/cellType.marker.dotplot.pdf", height = 3.5)
DotPlot(scATAC.data, cols = c("#1e90ff", "#F15F30"), features = cell.type.markers$Gene, group.by = "AnnotatedcellType", dot.scale = 3.5) + theme(axis.text.x = element_text(size = 6, angle = 90)) + theme(axis.text.y = element_text(size = 4))
DotPlot(scATAC.data, cols = c("#1e90ff", "#F15F30"), features = cell.type.markers$Gene, group.by = "AnnotatedcellType", dot.scale = 3.5) + NoLegend()+ theme(axis.text.x = element_text(size = 6, angle = 90)) + theme(axis.text.y = element_text(size = 4))
dev.off()

#### plot the number ratio of sample in each cell cell type
scRNA.merge <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/data.merge.pro.rds")
colors <- readRDS("2.Cluster/AnnotateCellType/cellType.color.rds")
ratio.plot(seurat.object = scATAC.data, id.vars1 = "AnnotatedcellType", id.vars2 = "orig.ident", color.len = colors)


#### plot the number ratio of cell types in each sample
names(colors) <- gsub("Endothelium VCAM1-", "Endothelium (VCAM1-)", names(colors))
names(colors) <- gsub("Endothelium VCAM1\\+", "Endothelium (VCAM1+)", names(colors))
library(plotly)
pdf("2.Cluster/AnnotateCellType/cellType.patient.ratio.pdf")
patients <- unique(scATAC.data$orig.ident)
res <- sapply(patients, function(x){
  a <- subset(scATAC.data, subset = orig.ident==x)
  cellType.ratio <- as.data.frame(table(as.character(a$AnnotatedcellType)))
  group.ratio <- cellType.ratio$Freq/sum(cellType.ratio$Freq)
  group.ratio <- data.frame(Type = cellType.ratio$Var1, Ratio = group.ratio)
  group.ratio$Type <- factor(group.ratio$Type, levels = names(colors))
  labs <- paste0(round(group.ratio$Ratio, 4)*100, "%")
  p <- ggpie(group.ratio, "Ratio", label = labs,
        fill = "Type", color = "white", lab.pos = "out",
        palette = "npgs") + scale_fill_manual(values = colors)
  print(p)
})
dev.off()
