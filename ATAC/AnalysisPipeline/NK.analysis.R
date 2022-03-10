#' @description: analysis NK cells

library(Signac)
library(Seurat)
library(harmony)
library(SeuratWrappers)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86) #---GRCh38 (hg38)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(clustree)
library(vegan)
library(dplyr)
library(tidyverse)
set.seed(101)
library(openxlsx)
library(future)
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 50000 * 1024^2) #50G
library(GenomicRanges)
library(ggpubr)
source(file = "/home/longzhilin/Analysis_Code/Visualization/colorPalettes.R")
source(file = "/home/longzhilin/Analysis_Code/SingleCell/scATAC.Integrate.multipleSample.R")
setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC")
source(file = "/home/longzhilin/Analysis_Code/code/ratio.plot.R")

scATAC <- readRDS("scATAC.data.pro.rds")
DefaultAssay(scATAC) <- "ATAC"
cell.type <- "NK cell"
sub.ATAC <- subset(scATAC, subset = AnnotatedcellType == cell.type)

sub.ATAC <- RunTFIDF(sub.ATAC)
sub.ATAC <- FindTopFeatures(sub.ATAC, min.cutoff = "q1")

pdf("8.Immune/NK/Harmony.integration.PC10.pdf")
sub.scATAC.harmony <- Harmony.integration.reduceDimension(scATAC.object = sub.ATAC, set.resolutions = seq(0.1, 1.2, by = 0.1), group.by.vars = "dataset", assay = "ATAC", PC = 10, npcs = 30)
dev.off()

sub.scATAC.harmony$seurat_clusters <- sub.scATAC.harmony$ATAC_snn_res.0.2
Idents(sub.scATAC.harmony) <- sub.scATAC.harmony$seurat_clusters

#### marker expression
DefaultAssay(sub.scATAC.harmony) <- "ACTIVITY"
features <- c("CD3D", "CD3E", "CD8A", "CD8B", "CD68", "CD163", "NCAM1", "NKG7", "KLRD1", "KLRB1", "FGFBP2", "FCGR3A", "GNLY")
pdf("8.Immune/NK/marker.expression.pdf")
DotPlot(sub.scATAC.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "seurat_clusters", dot.scale = 3.5) + xlab("") + ylab("") + theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.text.y = element_text(size = 10))
VlnPlot(sub.scATAC.harmony, features = features, group.by="seurat_clusters", same.y.lims=T,flip = T, stack = T) & xlab("") & ylab("Log-normalized expression") & theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.5))
FeaturePlot(sub.scATAC.harmony, features = c("CD3D", "CD8A", "KLRB1", "FGFBP2"), cols = c("lightgrey", "red"), ncol = 2)
dev.off()

DefaultAssay(sub.scATAC.harmony) <- "Peaks"
pdf("8.Immune/NK/marker.coverage.pdf")
CoveragePlot(
  object = sub.scATAC.harmony,
  region = c("CD3D", "CD8A", "KLRB1", "FGFBP2"),
  extend.upstream = 3000,
  extend.downstream = 3000
)
dev.off()

#### integrate scRNA and scATAC
source(file = "/home/longzhilin/Analysis_Code/SingleCell/Integrate.scRNA.scATAC.R")
source(file = "/home/longzhilin/Analysis_Code/SingleCell/variableFeatureSelection.R")
sub.scRNA.harmony <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/5.Immune/NK/sub.scRNA.harmony.rds")
DefaultAssay(sub.scRNA.harmony) <- "RNA"
VariableFeatures(sub.scRNA.harmony) <- sub.scRNA.harmony@assays$SCT@var.features
pdf("8.Immune/NK/Integrate.scRNA.scATAC.PC10.pdf")
transfer.anchors <- Integrate.scRNA.scATAC(scATAC.object = sub.scATAC.harmony, scRNA.object = sub.scRNA.harmony,
                                            scATAC.assay = "ACTIVITY", scRNA.assay = "RNA", ref.npcs = 30, feature.assay = "SCT", dims = 20,
                                            nfeatures = 2000, class.label = "cellType", scATAC.reduction = "harmony", PC = 10)
Coembedding <- Coembedding.scATAC.scRNA(transfer.anchors = transfer.anchors$transfer.anchors, 
                                         scATAC.object = transfer.anchors$scATAC.object, 
                                         scATAC.reduction = "harmony", PC = 20, dims = 30, observe.vars = c("dataset", "type"),
                                         scRNA.object = sub.scRNA.harmony, class.label = "cellType")
dev.off()