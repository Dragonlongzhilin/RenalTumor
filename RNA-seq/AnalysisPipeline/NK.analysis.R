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
