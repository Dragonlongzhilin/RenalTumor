#' @description Integrating scRNA-seq and scATAC-seq data

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
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(future)
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 50000 * 1024^2) #50G

source(file = "/home/longzhilin/Analysis_Code/Plot_colorPaletters.R")
source(file = "/home/longzhilin/Analysis_Code/SingleCell/Integrate.scRNA.scATAC.R")
setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC")
scRNA.merge <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/data.merge.pro.rds")

############################################ integrative model
#PC15
scATAC.merge <- readRDS("Harmony.scATAC.PC15.rds")
scATAC.merge$seurat_clusters <- scATAC.merge$ATAC_snn_res.0.5
Idents(scATAC.merge) <- scATAC.merge$seurat_clusters
DefaultAssay(scATAC.merge) <- "ACTIVITY"

source(file = "/home/longzhilin/Analysis_Code/code/ratio.plot.R")
pdf("2.Cluster/cluster.number.ratio.pdf", height = 4, width = 7)
ratio.plot(seurat.object = scATAC.merge, id.vars1 = "orig.ident", id.vars2 = "seurat_clusters")
dev.off()

pdf("3.Integration/integration.model.pdf")
DimPlot(scRNA.merge, cols = "#00BEC4", group.by = "Doublet") + NoLegend()
DimPlot(scATAC.merge, cols = "#F8766D", group.by = "high.tss") + NoLegend()
dev.off()

# low cell type resolution
pdf("3.Integration/PC15.feature3000.stardard/Integrate.scRNA.scATAC.PC15.low.pdf")
transfer.anchors1 <- Integrate.scRNA.scATAC(scATAC.object = scATAC.merge, scRNA.object = scRNA.merge,
                                            scATAC.assay = "ACTIVITY", ref.npcs = 40, dims = 30, PC = 15, scRNA.assay = "RNA",
                                            nfeatures = 3000, class.label = "cellType_low")
Coembedding1 <- Coembedding.scATAC.scRNA(transfer.anchors = transfer.anchors1$transfer.anchors, scRNA.assay = "RNA",
                                         scATAC.object = transfer.anchors1$scATAC.object, 
                                         scRNA.object = scRNA.merge, PC = 15, dims = 30, class.label = "cellType_low")
dev.off()
saveRDS(transfer.anchors1$transfer.anchors, "3.Integration/PC15.feature3000.stardard/transfer.anchors.feature3000.PC15.low.rds")
saveRDS(Coembedding1$scATAC.object, "3.Integration/PC15.feature3000.stardard/scATAC.object.feature3000.PC15.low.rds")
saveRDS(Coembedding1$coembed, "3.Integration/PC15.feature3000.stardard/Coembedding.feature3000.PC15.low.rds")

############################################ analysis predict score
setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC")
scATAC.data <- readRDS("3.Integration/PC15.feature3000.stardard/scATAC.object.feature3000.PC15.low.rds")

scATAC.data$seurat_clusters <- scATAC.data$ATAC_snn_res.0.5
scATAC.data.high <- subset(scATAC.data, subset = prediction.score.max>=0.5)
pdf("3.Integration/PC15.feature3000.stardard/highQuality.pdf")
DimPlot(scATAC.data.high, group.by = "cellType_low", label = TRUE, repel = TRUE) + NoLegend()
DimPlot(scATAC.data.high, group.by = "cellType_low", reduction = "tsne", label = TRUE, repel = TRUE) + NoLegend()
DimPlot(scATAC.data.high, group.by = "seurat_clusters", label = TRUE, repel = TRUE)+ NoLegend()
DimPlot(scATAC.data.high, group.by = "seurat_clusters", reduction = "tsne", label = TRUE, repel = TRUE) + NoLegend()
dev.off()

# the threshold 0.5 for prediction.score.max
low.quality <- subset(scATAC.data, subset = prediction.score.max < 0.5)
cluster.nums <- as.data.frame(table(low.quality$seurat_clusters))
cellType.nums <- as.data.frame(table(low.quality$cellType_low))
pdf("3.Integration/PC15.feature3000.stardard/low.quality.pdf")
ggbarplot(cluster.nums, x="Var1", y="Freq", fill = "Var1", color = "Var1",
          sort.by.groups=FALSE, 
          label = T, xlab = "", ylab = "Cell Number") + theme(legend.position="none")
ggbarplot(cellType.nums, x="Var1", y="Freq", fill = "Var1", color = "Var1",
          sort.by.groups=FALSE,
          label = T, xlab = "", ylab = "Cell Number") + theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1)) 
dev.off()

####---- plot the number ratio of each cell type
predict.matrix <- grep("prediction.score", colnames(scATAC.data@meta.data)) 
predict.matrix <- scATAC.data@meta.data[, predict.matrix]

scATAC.data$cellType_low <- factor(scATAC.data$cellType_low, levels = unique(scATAC.data$cellType_low))
max.predict.score <- tapply(scATAC.data$cellType_low, scATAC.data$seurat_clusters, table)
max.score.matrix <- matrix(unlist(max.predict.score), nrow = length(max.predict.score[[1]]))
rownames(max.score.matrix) <- names(max.predict.score[[1]])
colnames(max.score.matrix) <- paste0("Cluster ", 0:(length(max.predict.score)-1))
assign.ratio <- apply(max.score.matrix, 2, function(x){
    return(x/sum(x))
})
pdf("3.Integration/PC15.feature3000.stardard/integration.assign.cellType.ratio.pdf")
Heatmap(assign.ratio, name = "Cell type ratio", col = c("white", "red"), show_row_dend = F, show_column_dend = F, width = unit(10, "cm"), height = unit(8, "cm"))
dev.off()

#remove the cell with predict.score < 0.5
scATAC.data.test <- AddMetaData(scATAC.data, as.character(scATAC.data$cellType_low), "high.class")
index <- which(scATAC.data.test$prediction.score.max<0.5)
scATAC.data.test$high.class[index] <- "Low confidence"
scATAC.test <- scATAC.data.test
predict.matrix <- grep("prediction.score", colnames(scATAC.test@meta.data)) 
predict.matrix <- scATAC.test@meta.data[, predict.matrix]
scATAC.test$high.class <- factor(scATAC.test$high.class, levels = unique(scATAC.test$high.class))
max.predict.score <- tapply(scATAC.test$high.class, scATAC.test$seurat_clusters, table)
max.score.matrix <- matrix(unlist(max.predict.score), nrow = length(max.predict.score[[1]]))
rownames(max.score.matrix) <- names(max.predict.score[[1]])
colnames(max.score.matrix) <- paste0("Cluster ", 0:(length(max.predict.score)-1))
assign.ratio <- apply(max.score.matrix, 2, function(x){
    return(x/sum(x))
})
pdf("3.Integration/PC15.feature3000.stardard/integration.assign.cellType.ratio.high.pdf")
Heatmap(assign.ratio, name = "Cell type ratio", col = c("white", "red"), show_row_dend = F, show_column_dend = F, width = unit(10, "cm"), height = unit(8, "cm"))
dev.off()

scATAC.sub <- subset(scATAC.data, subset = prediction.score.max >= 0.5)
saveRDS(scATAC.sub, file = "scATAC.sub.rds")
