#' @description: proProcess, cluster and remove the batch effect

# load package
library(Seurat)
library(harmony)
library(clustree)
library(ggpubr)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(vegan)
set.seed(101)
library(future)
plan("multiprocess", workers = 10) 
options(future.globals.maxSize = 100000 * 1024^2) # set 50G RAM
setwd(dir = "/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA")
source(file = "/home/longzhilin/Analysis_Code/Visualization/colorPalettes.R")

#### four samples
T1.data <- Read10X(data.dir = "/data/raw_data/lzl/RenalTumor-20200713/cellRanger_result/scRNA-V5.0/T1/outs/filtered_feature_bc_matrix")
colnames(T1.data) <- paste0("T1_", colnames(T1.data))
#36601*8849
T2.data <- Read10X(data.dir = "/data/raw_data/lzl/RenalTumor-20200713/cellRanger_result/scRNA-V5.0/T2/outs/filtered_feature_bc_matrix")
colnames(T2.data) <- paste0("T2_", colnames(T2.data))
#36601*15440
T3.data <- Read10X(data.dir = "/data/raw_data/lzl/RenalTumor-20200713/cellRanger_result/scRNA-V5.0/T3/outs/filtered_feature_bc_matrix")
colnames(T3.data) <- paste0("T3_", colnames(T3.data))
#36601*14644
T4.data <- Read10X(data.dir = "/data/raw_data/lzl/RenalTumor-20200713/cellRanger_result/scRNA-V5.0/T4/outs/filtered_feature_bc_matrix")
colnames(T4.data) <- paste0("T4_", colnames(T4.data))
#36601*11016

#### Preliminary filtration
#min.cell >3 & min.features >200                        
T1 <- CreateSeuratObject(counts = T1.data,
                         project = "T1",
                         min.cells = 3,
                         min.features = 200)
#23137*8823                         
T2 <- CreateSeuratObject(counts = T2.data,
                         project = "T2",
                         min.cells = 3,
                         min.features = 200) 
#23487*15437                         
T3 <- CreateSeuratObject(counts = T3.data,
                         project = "T3",
                         min.cells = 3,
                         min.features = 200)
#24106*14602                         
T4 <- CreateSeuratObject(counts = T4.data,
                         project = "T4",
                         min.cells = 3,
                         min.features = 200) 
#23065*10792

#### Mitochondrial gene ratio
T1[["percent.mt"]] <- PercentageFeatureSet(T1, pattern = "^MT-")
T2[["percent.mt"]] <- PercentageFeatureSet(T2, pattern = "^MT-")
T3[["percent.mt"]] <- PercentageFeatureSet(T3, pattern = "^MT-")
T4[["percent.mt"]] <- PercentageFeatureSet(T4, pattern = "^MT-")

#### Draw a statistical graph of the number of genes/count number/proportion of mitochondrial genes
pdf(file = "1.QualityControl/count.feature.mt.pdf")
VlnPlot(T1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(T1, feature1 = "nCount_RNA", feature2 = "percent.mt")
ploT1 <- FeatureScatter(T1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + ploT1
VlnPlot(T2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(T2, feature1 = "nCount_RNA", feature2 = "percent.mt")
ploT1 <- FeatureScatter(T2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + ploT1
VlnPlot(T3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(T3, feature1 = "nCount_RNA", feature2 = "percent.mt")
ploT1 <- FeatureScatter(T3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + ploT1
VlnPlot(T4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(T4, feature1 = "nCount_RNA", feature2 = "percent.mt")
ploT1 <- FeatureScatter(T4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + ploT1

ggdensity(T1@meta.data, x = "nCount_RNA", title = "T1")
ggdensity(T1@meta.data, x = "nFeature_RNA", title = "T1")
ggdensity(T1@meta.data, x = "percent.mt", title = "T1")
ggdensity(T2@meta.data, x = "nCount_RNA", title = "T2")
ggdensity(T2@meta.data, x = "nFeature_RNA", title = "T2")
ggdensity(T2@meta.data, x = "percent.mt", title = "T2")
ggdensity(T3@meta.data, x = "nCount_RNA", title = "T3")
ggdensity(T3@meta.data, x = "nFeature_RNA", title = "T3")
ggdensity(T3@meta.data, x = "percent.mt", title = "T3")
ggdensity(T4@meta.data, x = "nCount_RNA", title = "T4")
ggdensity(T4@meta.data, x = "nFeature_RNA", title = "T4")
ggdensity(T4@meta.data, x = "percent.mt", title = "T4")
dev.off()

######################### Detect the resolution parameters of each sample cluster. After the parameters are determined, you can block them without executing [test]
set.resolutions <- seq(0.5, 2, by = 0.1)
pdf(file = "1.QualityControl/PCA-test.pdf")

#### T1
T1.pro <- subset(T1, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA > 1000 & nFeature_RNA < 6000) 
T1.pro <- SCTransform(T1.pro, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
T1.pro <- RunPCA(T1.pro, npcs = 100, verbose = FALSE)
ElbowPlot(object = T1.pro, ndims = 100)
T1.pro  <- FindNeighbors(object = T1.pro , dims = 1:50, verbose = FALSE)
T1.pro  <- FindClusters(object = T1.pro , resolution = set.resolutions, verbose = FALSE) 
clustree(T1.pro)
T1.pro  <- RunUMAP(T1.pro , dims = 1:50)
T1.res <- sapply(set.resolutions, function(x){
    p <- DimPlot(object = T1.pro, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
    print(p)
})

#### T2
T2.pro <- subset(T2, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA > 1000 & nFeature_RNA < 6000)
T2.pro <- SCTransform(T2.pro, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
T2.pro <- RunPCA(T2.pro, npcs = 100, verbose = FALSE)
ElbowPlot(object = T2.pro, ndims = 100)
T2.pro  <- FindNeighbors(object = T2.pro , dims = 1:50, verbose = FALSE)
T2.pro  <- FindClusters(object = T2.pro , resolution = set.resolutions, verbose = FALSE) 
clustree(T2.pro)
T2.pro  <- RunUMAP(T2.pro , dims = 1:50)
T2.res <- sapply(set.resolutions, function(x){
    p <- DimPlot(object = T2.pro, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
    print(p)
})

#### T3
T3.pro <- subset(T3, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA > 1000 & nFeature_RNA < 6000) 
T3.pro <- SCTransform(T3.pro, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
T3.pro <- RunPCA(T3.pro, npcs = 100, verbose = FALSE)
ElbowPlot(object = T3.pro, ndims = 100)
T3.pro  <- FindNeighbors(object = T3.pro , dims = 1:50, verbose = FALSE)
T3.pro  <- FindClusters(object = T3.pro , resolution = set.resolutions, verbose = FALSE) 
clustree(T3.pro)
T3.pro  <- RunUMAP(T3.pro , dims = 1:50)
T3.res <- sapply(set.resolutions, function(x){
    p <- DimPlot(object = T3.pro, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
    print(p)
})

#### T4
T4.pro <- subset(T4, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA > 1000 & nFeature_RNA < 6000) 
T4.pro <- SCTransform(T4.pro, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
T4.pro <- RunPCA(T4.pro, npcs = 100, verbose = FALSE)
ElbowPlot(object = T4.pro, ndims = 100)
T4.pro  <- FindNeighbors(object = T4.pro , dims = 1:50, verbose = FALSE)
T4.pro  <- FindClusters(object = T4.pro , resolution = set.resolutions, verbose = FALSE) 
clustree(T4.pro)
T4.pro  <- RunUMAP(T4.pro , dims = 1:50)
T4.res <- sapply(set.resolutions, function(x){
    p <- DimPlot(object = T4.pro, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
    print(p)
})
dev.off()

#### remove doublet
library(DoubletFinder) # Require cleanup of low-quality cells in advance
source(file = "/home/longzhilin/Analysis_Code/SingleCell/doubletDetect.R")
pdf("1.QualityControl/doublet.pdf")

T1.pro1 <- doubletDetect(Seurat.object = T1.pro, PCs = 1:50, doublet.rate = 0.061, annotation = "SCT_snn_res.0.7", sct = T) #7893 ~8000
T2.pro1 <- doubletDetect(Seurat.object = T2.pro, PCs = 1:50, doublet.rate = 0.106, annotation = "SCT_snn_res.1.7", sct = T) #13992 ~14000
T3.pro1 <- doubletDetect(Seurat.object = T3.pro, PCs = 1:50, doublet.rate = 0.091, annotation = "SCT_snn_res.0.7", sct = T) #11973 ~12000
T4.pro1 <- doubletDetect(Seurat.object = T4.pro, PCs = 1:50, doublet.rate = 0.061, annotation = "SCT_snn_res.0.5", sct = T) #8054 ~8000
dev.off()
saveRDS(T1.pro1, file = "T1.pro1.rds")
saveRDS(T2.pro1, file = "T2.pro1.rds")
saveRDS(T3.pro1, file = "T3.pro1.rds")
saveRDS(T4.pro1, file = "T4.pro1.rds")

pdf("1.QualityControl/doublet.cell.pdf")
DimPlot(object = T1.pro1, reduction = 'umap', group.by = "Doublet")
DimPlot(object = T2.pro1, reduction = 'umap', group.by = "Doublet")
DimPlot(object = T3.pro1, reduction = 'umap', group.by = "Doublet")
DimPlot(object = T4.pro1, reduction = 'umap', group.by = "Doublet")
dev.off()

T1.pro2 <- subset(T1.pro1, subset = Doublet == "Singlet") #21247*7449
T2.pro2 <- subset(T2.pro1, subset = Doublet == "Singlet") #21605*12574
T3.pro2 <- subset(T3.pro1, subset = Doublet == "Singlet") #21947*10937
T4.pro2 <- subset(T4.pro1, subset = Doublet == "Singlet") #20459*7605
saveRDS(T1.pro2, file = "T1.pro2.rds")
saveRDS(T2.pro2, file = "T2.pro2.rds")
saveRDS(T3.pro2, file = "T3.pro2.rds")
saveRDS(T4.pro2, file = "T4.pro2.rds")

############################################## merge data and correct the batch effect
DefaultAssay(T1.pro2) <- "RNA"
DefaultAssay(T2.pro2) <- "RNA"
DefaultAssay(T3.pro2) <- "RNA"
DefaultAssay(T4.pro2) <- "RNA"

source(file = "/home/longzhilin/Analysis_Code/SingleCell/variableFeatureSelection.R")
Renal.list <- list(T1 = T1.pro2, T2 = T2.pro2, T3 = T3.pro2, T4 = T4.pro2)
Renal.list.Stardard <- variableFeatureSelection(seurat.lists = Renal.list, method = "Stardard", nfeatures = 3000)
saveRDS(Renal.list.Stardard, file = "Renal.list.Stardard.3000.rds")

Renal.list.SCT <- variableFeatureSelection(seurat.lists = Renal.list, method = "SCT", nfeatures = 3000)
saveRDS(Renal.list.SCT, file = "Renal.list.SCT.3000.rds")

#### 
# assay=SCT
data.merge <- merge(Renal.list.SCT[[1]], y = Renal.list.SCT[2:length(Renal.list.SCT)], project = "Renal")
DefaultAssay(data.merge) <- "SCT"
seurat.features.SCT <- SelectIntegrationFeatures(object.list = Renal.list.SCT, nfeatures = 3000)
VariableFeatures(data.merge) <- seurat.features.SCT
# Remove previous clustering results
index <- match(paste0("SCT_snn_res.", seq(0.5, 2, by=0.1)), colnames(data.merge@meta.data))
data.merge@meta.data <- data.merge@meta.data[,-index]
# assay=RNA
seurat.features.RNA <- SelectIntegrationFeatures(object.list = Renal.list.Stardard, nfeatures = 3000)
DefaultAssay(data.merge) <- "RNA"
VariableFeatures(data.merge) <- seurat.features.RNA
data.merge <- NormalizeData(data.merge, verbose = FALSE)
data.merge <- ScaleData(data.merge, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mt"), features = rownames(data.merge@assays$RNA@data))
DefaultAssay(data.merge) <- "SCT"
saveRDS(data.merge, file = "data.merge.rds")

##############################################2.Evaluation of cellcycle and patient bias
DefaultAssay(data.merge) <- "SCT"
pdf("1.QualityControl/filtered.statistics.pdf")
VlnPlot(object = data.merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident", cols = Palettes$group_pal[1:length(unique(data.merge@meta.data$orig.ident))], pt.size = 0)
FeatureScatter(object = data.merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
FeatureScatter(object = data.merge, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
dev.off()
# Draw the distribution of the number of samples
cell.number <- as.data.frame(table(data.merge$orig.ident))
pdf("1.QualityControl/highQuality.pdf")
ggbarplot(cell.number, x="Var1", y="Freq", fill = "Var1", color = "Var1", palette = Palettes$group_pal[1:length(unique(data.merge@meta.data$orig.ident))],
          sort.by.groups=FALSE, #不按组排序
          label = T, xlab = "", ylab = "Cell Number") + theme(legend.position="none") 
dev.off()

#### Assess the cell cycle effect
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
data.merge <- CellCycleScoring(data.merge, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
data.merge <- RunPCA(data.merge, features = c(s.genes, g2m.genes))
pdf("1.QualityControl/cellCycle.afterMerge.pdf")
DimPlot(data.merge, dims = c(1, 2), reduction = "pca", group.by = "Phase")
DimPlot(data.merge, dims = c(1, 3), reduction = "pca", group.by = "Phase")
DimPlot(data.merge, dims = c(2, 3), reduction = "pca", group.by = "Phase")
dev.off()

set.resolutions <- seq(0.2, 1.2, by = 0.1)
#### First observe whether the clustering effect will depend on the sample
a <- data.merge
a <- RunPCA(a, npcs = 100, verbose = T)
pdf("1.QualityControl/merge.observe.batch.pdf")
ElbowPlot(object = a, ndims = 100)
a <- FindNeighbors(a, dims = 1:50, verbose = T)
a <- FindClusters(object = a, resolution = set.resolutions, verbose = T) 
clustree(a)
a <- RunUMAP(a, dims = 1:50)
DimPlot(object = a, reduction = 'umap',label = TRUE, group.by = "orig.ident")
DimPlot(object = a, reduction = 'umap',label = TRUE, group.by = "Phase")
merge.res <- sapply(set.resolutions, function(x){
    p <- DimPlot(object = a, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
    print(p)
})
dev.off()

##############################################3.correct batch effect
source(file = "/home/longzhilin/Analysis_Code/SingleCell/scRNA.Integrate.multipleSample.R")
#### SCT 3000; PC 50
pdf("2.Cluster/SCT.Harmony.Integration.PC50.feature3000-test.pdf")
data.merge.harmony.PC50.SCT <- Harmony.integration.reduceDimension(seurat.object = data.merge, assay = "SCT", set.resolutions = seq(0.2, 1.2, by = 0.1), PC = 50, nfeatures = 3000, npcs = 50)
dev.off()
saveRDS(data.merge.harmony.PC50.SCT, file = "data.merge.harmony.PC50.SCT.feature3000.rds")

pdf("2.Cluster/SCT.Harmony.Integration.PC40.feature3000.pdf")
data.merge.harmony.PC40.SCT <- Harmony.integration.reduceDimension(seurat.object = data.merge, assay = "SCT", set.resolutions = seq(0.2, 1.2, by = 0.1), PC = 40, nfeatures = 3000, npcs = 50)
dev.off()
saveRDS(data.merge.harmony.PC40.SCT, file = "data.merge.harmony.PC40.SCT.feature3000.rds")

pdf("2.Cluster/SCT.Harmony.Integration.PC30.feature3000.pdf")
data.merge.harmony.PC30.SCT <- Harmony.integration.reduceDimension(seurat.object = data.merge, assay = "SCT", set.resolutions = seq(0.2, 1.2, by = 0.1), PC = 30, nfeatures = 3000, npcs = 50)
dev.off()
saveRDS(data.merge.harmony.PC30.SCT, file = "data.merge.harmony.PC30.SCT.feature3000.rds")

#### FindVariableFeatures 3000; PC 50
pdf("2.Cluster/Stardard.Harmony.Integration.PC50.feature3000.pdf")
data.merge.harmony.PC50.stardard <- Harmony.integration.reduceDimension(seurat.object = data.merge, assay = "RNA", set.resolutions = seq(0.2, 1.2, by = 0.1), PC = 50, nfeatures = 3000, npcs = 50)
dev.off()
saveRDS(data.merge.harmony.PC50.stardard, file = "data.merge.harmony.stardard.PC50.feature3000.rds")

pdf("2.Cluster/Stardard.Harmony.Integration.PC30.feature3000.pdf")
data.merge.harmony.PC30.stardard <- Harmony.integration.reduceDimension(seurat.object = data.merge, assay = "RNA", set.resolutions = seq(0.2, 1.2, by = 0.1), PC = 30, nfeatures = 3000, npcs = 50)
dev.off()
saveRDS(data.merge.harmony.PC30.stardard, file = "data.merge.harmony.stardard.PC30.feature3000.rds")


#############################################4.Determine the best classification model
#Harmony
DefaultAssay(data.merge.harmony.PC50.stardard) <- "RNA"
DefaultAssay(data.merge.harmony.PC50.SCT) <- "RNA"
features <- c("CA9", "PECAM1", "MKI67", "CD8A", "CD4", "IL7R", "KLRD1", "FOXP3", "MS4A1", "IGKC", "CD68", "CD163", "CD14")
pdf("2.Cluster/AnnotateCellType/initial.cluster.harmony.stardard.pdf")
res <- sapply(features, function(x){
    p <- FeaturePlot(data.merge.harmony.PC50.stardard, features = x, cols = c("lightgrey", "red"), reduction = 'umap') 
    print(p)
    p <- FeaturePlot(data.merge.harmony.PC50.stardard, features = x, cols = c("lightgrey", "red"), reduction = 'tsne') 
    print(p)
})
dev.off()
pdf("2.Cluster/AnnotateCellType/initial.cluster.harmony.SCT.pdf")
res <- sapply(features, function(x){
    p <- FeaturePlot(data.merge.harmony.PC50.SCT, features = x, cols = c("lightgrey", "red"), reduction = 'umap') 
    print(p)
    p <- FeaturePlot(data.merge.harmony.PC50.SCT, features = x, cols = c("lightgrey", "red"), reduction = 'tsne') 
    print(p)
})
dev.off()
#Seurat
DefaultAssay(data.merge.seurat.PC50.stardard) <- "RNA"
DefaultAssay(data.merge.seurat.PC50.SCT) <- "RNA"
pdf("2.Cluster/AnnotateCellType/initial.cluster.seurat.stardard.pdf")
res <- sapply(features, function(x){
    p <- FeaturePlot(data.merge.seurat.PC50.stardard, features = x, cols = c("lightgrey", "red"), reduction = 'umap') 
    print(p)
    p <- FeaturePlot(data.merge.seurat.PC50.stardard, features = x, cols = c("lightgrey", "red"), reduction = 'tsne') 
    print(p)
})
dev.off()
pdf("2.Cluster/AnnotateCellType/initial.cluster.seurat.SCT.pdf")
res <- sapply(features, function(x){
    p <- FeaturePlot(data.merge.seurat.PC50.SCT, features = x, cols = c("lightgrey", "red"), reduction = 'umap') 
    print(p)
    p <- FeaturePlot(data.merge.seurat.PC50.SCT, features = x, cols = c("lightgrey", "red"), reduction = 'tsne') 
    print(p)
})
dev.off()