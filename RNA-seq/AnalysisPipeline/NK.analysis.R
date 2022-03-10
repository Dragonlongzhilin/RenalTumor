#' @description: analysis NK cell population

library(Seurat)
library(harmony)
library(clustree)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(vegan)
library(openxlsx)
set.seed(101)
library(future)
plan("multiprocess", workers = 5) 
options(future.globals.maxSize = 50000 * 1024^2) # set 50G RAM
source(file = "/home/longzhilin/Analysis_Code/Visualization/colorPalettes.R")
source(file = "/home/longzhilin/Analysis_Code/SingleCell/scRNA.Integrate.multipleSample.R")
source(file = "/home/longzhilin/Analysis_Code/Combined.P.FC.R")
source(file = "/home/longzhilin/Analysis_Code/SingleCell/variableFeatureSelection.R")
source(file = "/home/longzhilin/Analysis_Code/code/ratio.plot.R")
set.resolutions <- seq(0.1, 0.8, by = 0.1)

setwd(dir = "/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA")
data.merge <- readRDS("data.merge.pro.rds")

#### observate the marker expression
NK.markers <- c("CD3D", "CD8A", "KLRD1", "GNLY")
pdf("2.Cluster/AnnotateCellType/NK.marker.pdf")
FeaturePlot(data.merge, features = NK.markers, cols = c("lightgrey", "red"), ncol = 2)
dev.off()

cell.Type <- "NK/NKT cell"
sub.scRNA <- subset(data.merge, subset = cellType_low == cell.Type)
#### plot the sub-NK plot
cell.locations <- as.data.frame(sub.scRNA@reductions$umap@cell.embeddings)
cells <- rownames(cell.locations)[which(cell.locations$UMAP_1 < -4 & cell.locations$UMAP_2 > 0)]
sub.scRNA.NK <- subset(sub.scRNA, cells = cells)
pdf("2.Cluster/AnnotateCellType/NK.marker.pdf")
FeaturePlot(sub.scRNA.NK, features = NK.markers, cols = c("lightgrey", "red"), ncol = 2)
dev.off()


DefaultAssay(sub.scRNA) <- "RNA"
sub.scRNA <- DietSeurat(sub.scRNA, assays = "RNA")
sub.scRNA@meta.data <- sub.scRNA@meta.data[,-grep("SCT_snn", colnames(sub.scRNA@meta.data))]

# observe the batch effect
# sub.scRNA <- SCTransform(sub.scRNA, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
# sub.scRNA <- RunPCA(sub.scRNA, npcs = 30, verbose = FALSE) %>%
# RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
# FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
# FindClusters(resolution = set.resolutions, verbose = FALSE)
# p1 <- DimPlot(sub.scRNA, group.by="orig.ident")
# p2 <- DimPlot(sub.scRNA, group.by="SCT_snn_res.0.3")
# ggarrange(p1, p2)

# split the dataset into a list of two seurat objects (stim and CTRL)
sub.list <- SplitObject(sub.scRNA, split.by = "orig.ident")
sub.list <- lapply(X = sub.list, FUN = function(x) {
    x <- SCTransform(x, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
})
sub.scRNA <- merge(sub.list[[1]], y = sub.list[2:length(sub.list)], project = "Renal_NK")
DefaultAssay(sub.scRNA) <- "SCT"
seurat.features.SCT <- SelectIntegrationFeatures(object.list = sub.list)
VariableFeatures(sub.scRNA) <- seurat.features.SCT
pdf("5.Immune/NK/SCT.Harmony.Integration.PC20.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "SCT", set.resolutions = set.resolutions, PC = 20, npcs = 30, nfeatures = 2000)
dev.off()
sub.scRNA.harmony$seurat_clusters <- sub.scRNA.harmony$SCT_snn_res.0.2
Idents(sub.scRNA.harmony) <- sub.scRNA.harmony$seurat_clusters
DefaultAssay(sub.scRNA.harmony) <- "RNA"
sub.scRNA.harmony <- NormalizeData(sub.scRNA.harmony, verbose = FALSE)
sub.scRNA.harmony <- ScaleData(sub.scRNA.harmony, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mt"), features = rownames(sub.scRNA.harmony))
saveRDS(sub.scRNA.harmony, file = "5.Immune/NK/sub.scRNA.harmony.rds")

#### marker expression
features <- c("CD3D", "CD8A", "CD8B", "CD68", "CD163", "CD1D", "NCAM1", "NKG7", "KLRD1", "KLRB1", "FGFBP2", "FCGR3A", "GNLY")
pdf("5.Immune/NK/marker.expression.pdf")
DotPlot(sub.scRNA.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "seurat_clusters", dot.scale = 3.5) + xlab("") + ylab("") + theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.text.y = element_text(size = 10))
VlnPlot(sub.scRNA.harmony, features = features, group.by="seurat_clusters", same.y.lims=T,flip = T, stack = T) & xlab("") & ylab("Log-normalized expression") & theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.5))
FeaturePlot(sub.scRNA.harmony, features = c("CD3D", "CD8A", "KLRD1", "FGFBP2"), cols = c("lightgrey", "red"), ncol = 2)
dev.off()

#### differential analysis
DefaultAssay(sub.scRNA.harmony) <- "RNA"
cluster.DE <- FindAllMarkers(sub.scRNA.harmony, 
                            group.by = "seurat_clusters", 
                            test.use = "MAST", latent.vars = "orig.ident")
saveRDS(cluster.DE, file = "5.Immune/NK/cellType.DE.rds")

top.genes <- cluster.DE %>% filter(p_val_adj<0.05 & avg_log2FC>0.25) %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
pdf("5.Immune/NK/cluster.topgenes.pdf")
DoHeatmap(sub.scRNA.harmony, features = unique(top.genes$gene), size = 2) + NoLegend()
dev.off()

#### assign cell type
cellType <- sub.scRNA.harmony$seurat_clusters
cellType <- gsub("^0$", "NK cell-1", cellType)
cellType <- gsub("^1$", "NKT cell", cellType)
cellType <- gsub("^2$", "NK cell-2", cellType)
cellType <- gsub("^3$", "Mix", cellType)
sub.scRNA.harmony <- AddMetaData(sub.scRNA.harmony, cellType, "cellType")
pdf("5.Immune/NK/cellType.pdf")
DimPlot(sub.scRNA.harmony, group.by = "cellType", label.size = 5, label = T) + NoLegend()
DotPlot(sub.scRNA.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + xlab("") + ylab("") + theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.text.y = element_text(size = 10))
VlnPlot(sub.scRNA.harmony, features = features, group.by="cellType", same.y.lims=T,flip = T, stack = T) & xlab("") & ylab("Log-normalized expression") & theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.5))
FeaturePlot(sub.scRNA.harmony, features = c("CD3D", "CD8A", "KLRD1", "FGFBP2"), cols = c("lightgrey", "red"), ncol = 2)
dev.off()

features <- c("CD3D", "CD3E", "CD8A", "CD8B", "CD68", "NCAM1", "NKG7", "KLRD1", "KLRB1", "FGFBP2", "FCGR3A", "GNLY")
pdf("5.Immune/NK/cellType.marker.vlnplot.pdf", width = 2.5, height = 5)
VlnPlot(sub.scRNA.harmony, features = features, group.by="cellType", same.y.lims=T,flip = T, stack = T) & xlab("") & ylab("Log-normalized expression") & theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.5))
dev.off()


