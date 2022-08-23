#' @description: cluster and dimension reduction

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
library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) #50G
library(GenomicRanges)
library(ggpubr)
source("/home/longzhilin/Analysis_Code/SingleCell/scATAC.Integrate.multipleSample.R")
setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC")

scATAC.merge.pro <- readRDS("scATAC.merge.pro.rds")

#Normalization: Signac performs term frequency-inverse document frequency (TF-IDF) normalization. 
#This is a two-step normalization procedure, that both normalizes across cells to correct for differences in cellular sequencing depth, and across peaks to give higher values to more rare peaks.
scATAC.merge.pro <- RunTFIDF(scATAC.merge.pro)

#Feature selection：though we note that we see very similar results when using only a subset of features (try setting min.cutoff to ‘q75’ to use the top 25% all peaks)
scATAC.merge.pro <- FindTopFeatures(scATAC.merge.pro, min.cutoff = 'q1')

saveRDS(scATAC.merge.pro, file = "scATAC.merge.pro.rds")

#####################1.observe batch effect########################
#Dimension reduction: We next run singular value decomposition (SVD) on the TD-IDF matrix, using the features (peaks) selected above. 
scATAC.merge.pro <- RunSVD(scATAC.merge.pro)
scATAC.merge.pro <- RunUMAP(scATAC.merge.pro, dims = 2:30, reduction = 'lsi')
scATAC.merge.pro <- RunTSNE(scATAC.merge.pro, dims = 2:30, reduction = 'lsi')
pdf("2.Cluster/observe.patient.effect.pdf")
ElbowPlot(object = scATAC.merge.pro, ndims = 30, reduction = "lsi")
DepthCor(scATAC.merge.pro)
DimPlot(scATAC.merge.pro, group.by = 'dataset', pt.size = 0.1, reduction = 'umap')
DimPlot(scATAC.merge.pro, group.by = 'dataset', pt.size = 0.1, reduction = 'tsne')
dev.off()

#####################2.correct batch effect########################
DefaultAssay(scATAC.merge.pro) <- "ATAC"

pdf("2.Cluster/Harmony.integration.PC15.pdf")
Harmony.scATAC.PC15 <- Harmony.integration.reduceDimension(scATAC.object = scATAC.merge.pro, set.resolutions = seq(0.2, 1.2, by = 0.1), groups = "dataset", assay = "ATAC", PC = 15, npcs = 30)
dev.off()
saveRDS(Harmony.scATAC.PC15, file = "Harmony.scATAC.PC15.rds")
#####################3.gene activity########################
#PC 15
gene.activities <- GeneActivity(Harmony.scATAC.PC15)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
Harmony.scATAC.PC15[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
Harmony.scATAC.PC15 <- NormalizeData(
  object = Harmony.scATAC.PC15,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(Harmony.scATAC.PC15$nCount_ATAC)
)
saveRDS(Harmony.scATAC.PC15, file = "Harmony.scATAC.PC15.rds")
