#' @description: analysis macrophage

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
library(reshape2)
set.seed(101)
library(future)
plan("multiprocess", workers = 5)
options(future.globals.maxSize = 50000 * 1024^2) #50G
library(GenomicRanges)
library(ggpubr)
source(file = "/home/longzhilin/Analysis_Code/Plot_colorPaletters.R")
source("/home/longzhilin/Analysis_Code/SingleCell/scATAC.Integrate.multipleSample.R")
source(file = "/home/longzhilin/Analysis_Code/SingleCell/Integrate.scRNA.scATAC.R")
setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC")
source(file = "/home/longzhilin/Analysis_Code/code/ratio.plot.R")

scATAC <- readRDS("scATAC.data.pro.rds")
DefaultAssay(scATAC) <- "ATAC"
cell.type <- "Macrophage"
sub.ATAC <- subset(scATAC, subset = AnnotatedcellType == cell.type)
####--------------------------------------------------------------------------- 1. clustering cell
sub.ATAC <- RunTFIDF(sub.ATAC)
sub.ATAC <- FindTopFeatures(sub.ATAC, min.cutoff = "q1")

pdf("8.Immune/Macrophage/Harmony.integration.pdf")
sub.scATAC.harmony <- Harmony.integration.reduceDimension(scATAC.object = sub.ATAC, set.resolutions = seq(0.2, 1.2, by = 0.1), groups = "dataset", assay = "ATAC", PC = 8, npcs = 10)
dev.off()
sub.scATAC.harmony$seurat_clusters <- sub.scATAC.harmony$ATAC_snn_res.0.2
Idents(sub.scATAC.harmony) <- sub.scATAC.harmony$seurat_clusters
saveRDS(sub.scATAC.harmony, file = "8.Immune/Macrophage/sub.scATAC.harmony.rds")

####--------------------------------------------------------------------------- 2. marker expression
DefaultAssay(sub.scATAC.harmony) <- "ACTIVITY"
pdf("8.Immune/Macrophage/markerExpression.pdf")
VlnPlot(sub.scATAC.harmony, features = c("CD3D", "PLVAP", "CD68", "CD163"), ncol = 2)
FeaturePlot(sub.scATAC.harmony, features = c("CD3D", "PLVAP", "CD68", "CD163"), cols = c("lightgrey", "red"), ncol = 2)
VlnPlot(sub.scATAC.harmony, features = c("APOE", "C1QA", "C1QB", "C1QC"), ncol = 2)
FeaturePlot(sub.scATAC.harmony, features = c("APOE", "C1QA", "C1QB", "C1QC"), cols = c("lightgrey", "red"), ncol = 2)
VlnPlot(sub.scATAC.harmony, features = c("SPP1", "CSTB", "FABP5", "FN1"), ncol = 2)
FeaturePlot(sub.scATAC.harmony, features = c("SPP1", "CSTB", "FABP5", "FN1"), cols = c("lightgrey", "red"), ncol = 2)
VlnPlot(sub.scATAC.harmony, features = c("S100A8", "S100A12", "FCGR3A", "FCN1"), ncol = 2)
FeaturePlot(sub.scATAC.harmony, features = c("S100A8", "S100A12", "FCGR3A", "FCN1"), cols = c("lightgrey", "red"), ncol = 2)
dev.off()
features <- c("CD3D", "CD3E", "PLVAP", "ESM1", "ITGAM", "CD68", "CD163", "CSF1R")
pdf("8.Immune/Macrophage/cluster.markerExpression.dot.pdf", height =  2, width = 5)
DotPlot(sub.scATAC.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "seurat_clusters", dot.scale = 3.5) + theme(axis.text.x = element_text(size = 7, angle = 90)) + theme(axis.text.y = element_text(size = 8))
DotPlot(sub.scATAC.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "seurat_clusters", dot.scale = 3.5) + NoLegend() + theme(axis.text.x = element_text(size = 7)) + theme(axis.text.y = element_text(size = 8))
dev.off()
pdf("8.Immune/Macrophage/cluster.markerExpression.dot.test.pdf", height =  4, width = 5)
DotPlot(sub.scATAC.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "seurat_clusters", dot.scale = 3.5) + theme(axis.text.x = element_text(size = 7, angle = 90)) + theme(axis.text.y = element_text(size = 8))
DotPlot(sub.scATAC.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "seurat_clusters", dot.scale = 3.5) + NoLegend() + theme(axis.text.x = element_text(size = 7, angle = 90)) + theme(axis.text.y = element_text(size = 8))
dev.off()

pdf("8.Immune/Macrophage/cluster.ratio.pdf", height = 3, width = 5)
ratio.plot(seurat.object = sub.scATAC.harmony, id.vars1 = "orig.ident", id.vars2 = "seurat_clusters", angle = 60)
ratio.plot(seurat.object = sub.scATAC.harmony, id.vars1 = "orig.ident", id.vars2 = "seurat_clusters", angle = 60)
dev.off()

pdf("8.Immune/Macrophage/cluster.markerExpression.coveragePlot.pdf")
DefaultAssay(sub.scATAC.harmony) <- "Peaks"
CoveragePlot(
  object = sub.scATAC.harmony,
  region = c("CD3D", "CD3E", "APOE", "C1QC"),
  extend.upstream = 1000,
  extend.downstream = 1000
)
CoveragePlot(
  object = sub.scATAC.harmony,
  region = c("CD68", "CD163", "PLVAP", "ESM1"),
  extend.upstream = 1000,
  extend.downstream = 1000
)
dev.off()

####--------------------------------------------------------------------------- 3.remove the cluster with cell < 100 and confounding cells
sub.scATAC.harmony <- subset(sub.scATAC.harmony, subset = seurat_clusters %in% c(0,1,3))
pdf("8.Immune/Macrophage/Harmony.integration.pro.pdf")
sub.scATAC.harmony <- Harmony.integration.reduceDimension(scATAC.object = sub.scATAC.harmony, set.resolutions = seq(0.1, 1.2, by = 0.1), groups = "dataset", assay = "ATAC", PC = 5, npcs = 10)
dev.off()
sub.scATAC.harmony$seurat_clusters <- sub.scATAC.harmony$ATAC_snn_res.0.1
Idents(sub.scATAC.harmony) <- sub.scATAC.harmony$seurat_clusters
#定义cellType
cellType <- gsub("^0$", "C3", sub.scATAC.harmony$seurat_clusters)
cellType <- gsub("^1$", "C1", cellType)
cellType <- gsub("^2$", "C2", cellType)
sub.scATAC.harmony$cellType2 <- factor(cellType, levels = c("C1", "C2", "C3"))
Idents(sub.scATAC.harmony) <- sub.scATAC.harmony$cellType2

cellType3 <- sub.scATAC.harmony$cellType2
cellType3 <- gsub("^C1$", "TAM-C1QB", cellType3)
cellType3 <- gsub("^C2$", "TAM-RGCC", cellType3)
cellType3 <- gsub("^C3$", "TAM-LGALS3", cellType3)
sub.scATAC.harmony$cellType3 <- factor(cellType3, levels = c("TAM-C1QB", "TAM-RGCC", "TAM-LGALS3"))
Idents(sub.scATAC.harmony) <- sub.scATAC.harmony$cellType3
saveRDS(sub.scATAC.harmony, file = "8.Immune/Macrophage/sub.scATAC.harmony.pro.rds")

#### 选定分类group
sub.scATAC.harmony$cellType <- sub.scATAC.harmony$cellType3
Idents(sub.scATAC.harmony) <- sub.scATAC.harmony$cellType
cellType.colors <- c("#F8766D", "#00BA38", "#619CFF")

DefaultAssay(sub.scATAC.harmony) <- "ACTIVITY"
pdf("8.Immune/Macrophage/cluster.pro.pdf")
DimPlot(object = sub.scATAC.harmony, label = TRUE, group.by = "seurat_clusters", pt.size = 1.5) + NoLegend()
DimPlot(object = sub.scATAC.harmony, label = TRUE, group.by = "orig.ident", pt.size = 1.5) + NoLegend()
dev.off()
pdf("8.Immune/Macrophage/cellType.pdf")
DimPlot(object = sub.scATAC.harmony, label = TRUE, group.by = "cellType", pt.size = 1.5) + NoLegend()
DimPlot(object = sub.scATAC.harmony, label = TRUE, group.by = "orig.ident", pt.size = 1.5) + NoLegend()
dev.off()

features <- c("IER2", "JUN", "F13A1", "APOE", "C1QB", "C1QC", "C1QA",
              "EREG", "NLRP3", "AREG", "SLC2A3", "RGCC", "CLEC5A",
              "LILRB4", "MARCO", "ANXA2", "LGALS3", "GPNMB", "TREM2")
pdf("8.Immune/Macrophage/cellType.markerExpression.dot.pdf", height = 1.5, width = 6)
DotPlot(sub.scATAC.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + theme(axis.text.x = element_text(size = 7, angle = 90)) + theme(axis.text.y = element_text(size = 8))
DotPlot(sub.scATAC.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + NoLegend() + theme(axis.text.x = element_text(size = 7, angle = 90)) + theme(axis.text.y = element_text(size = 8))
dev.off()
pdf("8.Immune/Macrophage/cellType.markerExpression.dot.test.pdf", height = 4, width = 6)
DotPlot(sub.scATAC.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + theme(axis.text.x = element_text(size = 7, angle = 90)) + theme(axis.text.y = element_text(size = 8))
DotPlot(sub.scATAC.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + NoLegend() + theme(axis.text.x = element_text(size = 7, angle = 90)) + theme(axis.text.y = element_text(size = 8))
dev.off()

#### integrate scRNA and scATAC
sub.scATAC.harmony <- AddMetaData(sub.scATAC.harmony, sub.scATAC.harmony$cellType, "ATAC.pupulation")
sub.scRNA <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/5.Immune/Macrophage/sub.scRNA.harmony.pro.rds")
DefaultAssay(sub.RNA) <- "RNA"
sub.scRNA$cellType <- sub.scRNA$cellType3
Idents(sub.scRNA) <- sub.scRNA$cellType
pdf("8.Immune/Macrophage/Integrate.scRNA.scATAC.pdf")
transfer.anchors1 <- Integrate.scRNA.scATAC(scATAC.object = sub.scATAC.harmony, scRNA.object = sub.scRNA,
                                            scATAC.assay = "ACTIVITY", ref.npcs = 30, dims = 30, PC = 5,
                                            nfeatures = 3000, class.label = "cellType", SCT = T)
Coembedding1 <- Coembedding.scATAC.scRNA(transfer.anchors = transfer.anchors1$transfer.anchors, 
                                         scATAC.object = transfer.anchors1$scATAC.object, 
                                         scRNA.object = sub.scRNA, PC = 5, dims = 20, class.label = "cellType")
dev.off()
sub.scATAC.harmony$predicted.id <- Coembedding1$scATAC.object$predicted.id
sub.scATAC.harmony$prediction.score.max <- Coembedding1$scATAC.object$prediction.score.max

####------------------------------------------ 2.3. analysis cell state based on the functional gene sets
## 1.VISION method
source(file = "/home/longzhilin/Analysis_Code/SingleCell/vision_seurat.R")
source(file = "/home/longzhilin/Analysis_Code/SingleCell/vision.plot.R")
require(VISION)
library(reshape2)
options(mc.cores = 36)
signature.genes <- read.table("/data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/M1&M2&angiogenic&phagocytic.Cheng.2021.Cell.txt", header = T, stringsAsFactors = F, sep = "\t")

#### construct the vision calculate mode
# https://yoseflab.github.io/VISION/articles/web_only/Seurat.html
# Expression data should be scaled and normalized, but not log-transformed
# obj@[[assay]]@counts is used as the expression input (after normalizing to a library size of 10,000)
exp.data <- GetAssayData(sub.scATAC.harmony, slot = "data", assay = "ACTIVITY")
#remove the genes expressed less than 3 cells 
idx <- match(signature.genes$Gene, rownames(exp.data))
signature.genes <- signature.genes[which(!is.na(idx)),]
signature.names <- unique(signature.genes$Type)
#构建signature set
signature.genesets <- c()
positive <- 1
for(i in names(table(signature.genes$Type))){
    genes <- signature.genes$Gene[which(signature.genes$Type==i)]
    sigData <- rep(positive, length(genes))
    names(sigData) <- genes
    sig <- createGeneSignature(name = i, sigData = sigData)
    signature.genesets <- c(signature.genesets, sig)
}
vision.obj <- Vision(exp.data, signatures = signature.genesets, min_signature_genes = 4)
vision.obj <- analyze(vision.obj)
vision.obj@metaData$cellType <- sub.scATAC.harmony$cellType
saveRDS(vision.obj, file = "8.Immune/Macrophage/vision_signatureScore.rds")

library(reshape2)
sigScore <- as.data.frame(vision.obj@SigScores)
sigScore$group <- sub.scATAC.harmony$cellType
vision.group.score <- melt(sigScore, id.vars = "group", variable.name = "Type", value.name = "Score")
pdf("8.Immune/Macrophage/vision_signatureScore.pdf")
p1 <- ggboxplot(vision.group.score, x = "group", y = "Score", palette = "npg", color = "Type", xlab = "", title = "", ylab = "Score") + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + stat_compare_means(aes(group = Type), label =  "p.signif")
print(p1)

# heatmap
vision.score <- sigScore[,-ncol(sigScore)]
vision.score.mean <- apply(vision.score, 2, function(x){
    score <- tapply(x, sub.scATAC.harmony$cellType, mean)
    return(score)
})
p <- Heatmap(t(vision.score.mean), width = unit(6, "cm"), height = unit(6, "cm"), name = "Signature score",
        show_row_dend = F, show_column_dend = F, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
print(p)
vision.score.mean.scale <- scale(vision.score.mean) # 比较某个状态在各个细胞类型的情况
p <- Heatmap(t(vision.score.mean.scale), width = unit(6, "cm"), height = unit(6, "cm"), name = "Signature score",
        show_row_dend = F, show_column_dend = F, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
print(p)

# violin
plot.list <- lapply(colnames(vision.score), function(x){
    my_comparisons <- as.list(as.data.frame(combn(levels(sub.scATAC.harmony$cellType),2)))
    a <- sigScore[, c(x, "group")]
    names(a) <- c("Signature score", "group")
    p <- ggboxplot(a, x = "group", y = "Signature score", title = x, color = "black", alpha = 0.8, fill = "group", add = "jitter", add.params = list(fill = "black", size = 0.1)) + 
         xlab("") + rotate_x_text(angle = 45, vjust = 1) + scale_fill_manual(values = cellType.colors) + 
         stat_compare_means(comparisons = my_comparisons, label = "p.signif") + NoLegend()
    return(p)
})
p <- ggarrange(plotlist = plot.list[1:4],ncol = 2, nrow = 2)
print(p)
dev.off()

## check Phagocytosis.marker peak coverage plot
source("/home/longzhilin/Analysis_Code/SingleCell/FindRegion.R")
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
genes <-  c("MRC1", "CD163", "MERTK", "C1QB")
pdf("8.Immune/Macrophage/Phagocytosis.markerExpression.coveragePlot.pdf")
res <- sapply(genes, function(x){
  cat(x, "...\n")
  regions <- FindRegion(object = sub.scATAC.harmony, region = x, assay = "Peaks", extend.upstream = 1000, extend.downstream = 1000)
  idx <- data.frame(findOverlaps(regions, promoter))
  if(nrow(idx)>0){
    p <- CoveragePlot(
      object = sub.scATAC.harmony,
      region = x,
      ranges.title = "MACS2",
      links = F,
      peaks = T,
      extend.upstream = 1000,
      extend.downstream = 1000,
      region.highlight = promoter[idx[,2],])
  }else{
    p <- CoveragePlot(
      object = sub.scATAC.harmony,
      region = x,
      ranges.title = "MACS2",
      links = F,
      extend.upstream = 1000,
      extend.downstream = 1000,
      peaks = T)
  }
  print(p)
  return(regions)
})
dev.off()
pdf("8.Immune/Macrophage/Phagocytosis.markerExpression.coveragePlot.pro.pdf", width = unit(3, "inches"))
res <- sapply(genes, function(x){
  cat(x, "...\n")
  regions <- FindRegion(object = sub.scATAC.harmony, region = x, assay = "Peaks", extend.upstream = 1000, extend.downstream = 1000)
  idx <- data.frame(findOverlaps(regions, promoter))
  if(nrow(idx)>0){
    p <- CoveragePlot(
      object = sub.scATAC.harmony,
      region = x,
      ranges.title = "MACS2",
      links = F,
      peaks = T,
      extend.upstream = 1000,
      extend.downstream = 1000,
      region.highlight = promoter[idx[,2],])
  }else{
    p <- CoveragePlot(
      object = sub.scATAC.harmony,
      region = x,
      ranges.title = "MACS2",
      links = F,
      extend.upstream = 1000,
      extend.downstream = 1000,
      peaks = T)
  }
  print(p)
  return(regions)
})
dev.off()

DefaultAssay(sub.scATAC.harmony) <- "ACTIVITY"
####----plot the expression heatmap for M1 and M2 signature 
library(reshape2)
signature.genes <- signature.genes[which(signature.genes$Type %in% c("M1-state", "M2-state")),]
idx <- match(signature.genes$Gene, rownames(sub.scATAC.harmony))
signature.genes <- signature.genes[which(!is.na(idx)),]
signature.names <- unique(signature.genes$Type)
signature.list <- lapply(signature.names, function(x){
    genes <- signature.genes$Gene[which(signature.genes$Type == x)]
    return(genes)
})
names(signature.list) <- signature.names
#signature gene matrix
exp.matrix <- AverageExpression(sub.scATAC.harmony, signature.genes$Gene, assay = "ACTIVITY")
exp.matrix <- exp.matrix$ACTIVITY
col.split <- signature.genes$Type
exp.matrix.scale <- scale(t(exp.matrix))
pdf("8.Immune/Macrophage/M1&M2.signature.score.pdf")
Heatmap(t(exp.matrix), name = "Expression", cluster_columns = T, column_split = col.split, cluster_rows = T, show_row_dend = F, show_column_dend = F, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 7), width = unit(10, "cm"), height = unit(3, "cm"))
Heatmap(exp.matrix.scale, name = "Expression", cluster_columns = T, column_split = col.split, cluster_rows = T, show_row_dend = F, show_column_dend = F, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 7), width = unit(10, "cm"), height = unit(3, "cm"))
dev.off()

####----other signatrue genes----####
Immune.checkpoint.evasion.genes <- c("ICOSLG", "CD80", "CD86", "VSIR", "VSIG4", "LGALS9", "CD274", "PDCD1LG2", "SIGLEC10")
pdf("8.Immune/Macrophage/Immune.checkpoint.evasion.genes.pdf", height = 1.5, width = 5)
DotPlot(sub.scATAC.harmony, features = Immune.checkpoint.evasion.genes, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + theme(axis.text.x = element_text(size = 7, angle = 90)) + theme(axis.text.y = element_text(size = 8))
DotPlot(sub.scATAC.harmony, features = Immune.checkpoint.evasion.genes, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + NoLegend() + theme(axis.text.x = element_text(size = 7, angle = 90)) + theme(axis.text.y = element_text(size = 8))
dev.off()
pdf("8.Immune/Macrophage/Immune.checkpoint.evasion.genes.test.pdf", height = 4, width = 5)
DotPlot(sub.scATAC.harmony, features = Immune.checkpoint.evasion.genes, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + theme(axis.text.x = element_text(size = 7, angle = 90)) + theme(axis.text.y = element_text(size = 8))
DotPlot(sub.scATAC.harmony, features = Immune.checkpoint.evasion.genes, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + NoLegend() + theme(axis.text.x = element_text(size = 7, angle = 90)) + theme(axis.text.y = element_text(size = 8))
dev.off()
pdf("8.Immune/Macrophage/Immune.checkpoint.coveragePlot.pdf", width = unit(3, "inches"))
DefaultAssay(sub.scATAC.harmony) <- "Peaks"
res <- sapply(Immune.checkpoint.evasion.genes[-c(4)], function(x){
  cat(x, "...\n")
  regions <- FindRegion(object = sub.scATAC.harmony, region = x, assay = "Peaks", extend.upstream = 1000, extend.downstream = 1000)
  idx <- data.frame(findOverlaps(regions, promoter))
  if(nrow(idx)>0){
    p <- CoveragePlot(
      object = sub.scATAC.harmony,
      region = x,
      links = F,
      peaks = T,
      extend.upstream = 1000,
      extend.downstream = 1000,
      region.highlight = promoter[idx[,2],])
  }else{
    p <- CoveragePlot(
      object = sub.scATAC.harmony,
      region = x,
      links = F,
      extend.upstream = 1000,
      extend.downstream = 1000,
      peaks = T)
  }
  print(p)
  return(regions)
})
dev.off()


###--- MHC molecular
signature.genes <- read.table("/data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/antigenPresenting.LeiZhang.2021.Cell.txt", header = T, stringsAsFactors = F, sep = "\t")

idx <- match(signature.genes$Gene, rownames(sub.scATAC.harmony))
signature.genes <- signature.genes[which(!is.na(idx)),]
exp.matrix <- AverageExpression(sub.scATAC.harmony, signature.genes$Gene, assay = "ACTIVITY")
exp.matrix <- exp.matrix$ACTIVITY
exp.matrix.scale <- scale(t(exp.matrix))
col.split <- signature.genes$Type
pdf("8.Immune/Macrophage/MHC.expression.pdf")
Heatmap(exp.matrix, name = "Expression", cluster_rows = T, row_split = col.split, cluster_columns = T, show_row_dend = F, show_column_dend = F, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8), height = unit(8, "cm"), width = unit(3, "cm"))
Heatmap(t(exp.matrix.scale), name = "Expression", cluster_rows = T, row_split = col.split, cluster_columns = T, show_row_dend = F, show_column_dend = F, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8), height = unit(8, "cm"), width = unit(3, "cm"))
DotPlot(sub.scATAC.harmony, features = signature.genes$Gene, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 4) + RotatedAxis() + theme(axis.text.y = element_text(size = 8)) + labs(x= "", y = "") + theme(axis.text.x = element_text(size = 8))
dev.off()

####--------------------------------------------------------------------------- 4.motif analysis
motif.info <- data.frame(originName = names(sub.scATAC.harmony@assays$Peaks@motifs@motif.names), TF = unlist(sub.scATAC.harmony@assays$Peaks@motifs@motif.names))
rownames(motif.info) <- NULL
motif.info$originName <- gsub("_", "-", motif.info$originName)
library(openxlsx)
library(chromVARmotifs)
library(BSgenome.Hsapiens.UCSC.hg38)
data("human_pwms_v2")
source(file = "/home/longzhilin/Analysis_Code/Combined.P.FC.R")
DefaultAssay(sub.scATAC.harmony) <- 'chromvar'
library(tidyverse)
GetChromvarActivities <- function(cellType, scATAC.object, motif.info) {
  print(paste0("Finding chromVAR activities for: ",cellType))
  differential.activity <- FindMarkers(scATAC.object, 
                     ident.1 = cellType,
                     test.use = 'LR',
                     logfc.threshold = 0,
                     latent.vars = "nCount_Peaks") 
  motifs <- gsub("-", "_", rownames(differential.activity))
  motifNames <- sapply(motifs, function(x) motif.info@motif.names[[x]])
  return(cbind(differential.activity, gene = motifNames))
}
idents <- as.character(levels(Idents(sub.scATAC.harmony)))
sub.motifs.chromVAR <- lapply(idents, function(x) GetChromvarActivities(x, scATAC.object = sub.scATAC.harmony, motif.info = sub.scATAC.harmony@assays$Peaks@motifs))
names(sub.motifs.chromVAR) <- idents
sub.motifs.chromVAR <- lapply(sub.motifs.chromVAR, function(x){
  pi <- Combined.P.FC(x[,c("avg_log2FC", "p_val_adj")], log10P = F)
  x$pi <- pi$pi
  up.x <- x %>% filter(avg_log2FC>0) %>% arrange(desc(avg_log2FC))
  down.x <- x %>% filter(avg_log2FC<=0) %>% arrange(avg_log2FC)
  x <- rbind(up.x, down.x)
  return(x)
})
names(sub.motifs.chromVAR) <- idents
write.xlsx(sub.motifs.chromVAR, file = "8.Immune/Macrophage/sub.motifs.chromVAR.human_pwms_v2.xlsx", sheetName = idents, rowNames = T)
saveRDS(sub.motifs.chromVAR, file = "8.Immune/Macrophage/sub.motifs.chromVAR.human_pwms_v2.rds")

sub.motifs.chromVAR.sig <- lapply(sub.motifs.chromVAR, function(x){
  x <- x %>% filter(avg_log2FC>0.25 & p_val_adj<0.05) %>% arrange(desc(pi))
  return(x)
})
names(sub.motifs.chromVAR.sig) <- names(sub.motifs.chromVAR)
write.xlsx(sub.motifs.chromVAR.sig, file = "8.Immune/Macrophage/sub.motifs.chromVAR.sig.human_pwms_v2.xlsx", sheetName = idents, rowNames = T)

sub.motifs.chromVAR.sig <- lapply(names(sub.motifs.chromVAR.sig), function(x){
    sub.motifs.chromVAR.sig[[x]]$Type <- rep(x, nrow(sub.motifs.chromVAR.sig[[x]]))
    return(sub.motifs.chromVAR.sig[[x]])
})
names(sub.motifs.chromVAR.sig) <- names(sub.motifs.chromVAR)
sig.chromVAR.TFs <- sub.motifs.chromVAR.sig %>% bind_rows()

####---- Extract TF specific to each cell type
specific.TFs <- lapply(names(sub.motifs.chromVAR.sig), function(x){
  TFs <- sig.chromVAR.TFs[which(sig.chromVAR.TFs$Type == x),]
  others <- sig.chromVAR.TFs[which(sig.chromVAR.TFs$Type %in% setdiff(levels(sub.scATAC.harmony), x)),]
  genes <- setdiff(TFs$gene, others$gene)
  idx <- which(genes %in% TFs$gene)
  return(TFs[idx,])
})
top10 <- sapply(specific.TFs, function(x){
    return(x$gene[1:10])
})
top10 <- unique(as.character(top10))

TFs.set <- unique(sig.chromVAR.TFs$gene)
idx <- match(TFs.set, motif.info$TF)
sig.motifNames <- motif.info[idx,]

motifs.avgExp <- AverageExpression(sub.scATAC.harmony, features = sig.motifNames$originName, assays = "chromvar")
motifs.avgExp <- motifs.avgExp$chromvar
rownames(motifs.avgExp) <- sig.motifNames$TF
zScore <- function(x){(x - mean(x)) /sd(x)}
motifs.avgExp.scale <- apply(motifs.avgExp, 1, zScore) %>% t() # row: TF; column: cell type
mark.idx <- match(top10, rownames(motifs.avgExp.scale))
pdf("8.Immune/Macrophage/TAM.sig.motifs.heatmap.pdf")
Heatmap(motifs.avgExp.scale, name = "Deviation score", 
        width = unit(3, "cm"), height = unit(12, "cm"),
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
        show_row_dend = F, show_column_dend = F, show_column_names = T, show_row_names = F, 
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 8), by_row = T))+ 
rowAnnotation(link = anno_mark(at = mark.idx, labels = top10, link_width = unit(2, "mm"), labels_gp = gpar(fontsize = 5), padding = unit(1, "mm")))
dev.off()

# combined the DEGs
sub.motifs.chromVAR.sig <- lapply(sub.motifs.chromVAR, function(x){
  x <- x %>% filter(abs(avg_log2FC)>0.25 & p_val_adj<0.05) %>% arrange(desc(pi))
  return(x)
})
names(sub.motifs.chromVAR.sig) <- names(sub.motifs.chromVAR)
sub.scRNA.DEGs <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/5.Immune/Macrophage/cellType.DE.pro.rds")
sub.scRNA.DEGs <- filter(sub.scRNA.DEGs, abs(avg_log2FC) > 0.25 & p_val_adj < 0.05)
common.DEGs.TFs <- lapply(names(sub.motifs.chromVAR.sig), function(x){
  TFs <- sub.motifs.chromVAR.sig[[x]]
  idx <- which(sub.scRNA.DEGs$cluster == x)
  sub.scRNA.DEGs.cluster <- sub.scRNA.DEGs[idx,]
  all.info <- merge(TFs, sub.scRNA.DEGs.cluster, by = "gene", all = F) %>% arrange(desc(pi.x))
  colnames(all.info) <- gsub("\\.x", ".scATAC", colnames(all.info))
  colnames(all.info) <- gsub("\\.y", ".scRNA", colnames(all.info))
  return(all.info)
})
names(common.DEGs.TFs) <- names(sub.motifs.chromVAR.sig)
common.DEGs.TFs <- common.DEGs.TFs %>% bind_rows()

pdf("8.Immune/Macrophage/scATAC.scRNA.TFs.pdf")
row_split <- common.DEGs.TFs$cluster
rownames(common.DEGs.TFs) <- paste0(common.DEGs.TFs$cluster, "_", common.DEGs.TFs$gene)
col_fun <- colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))
p1 <- Heatmap(common.DEGs.TFs[, "avg_log2FC.scATAC",drop = F], col = col_fun, name = "avg_log2FC (scATAC)", row_split = row_split, row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 8), cluster_rows = T, show_row_dend = F, cluster_columns = F, width = unit(0.5, "cm"), height = unit(8, "cm"))
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
p2 <- Heatmap(common.DEGs.TFs[, "avg_log2FC.scRNA",drop = F], name = "avg_log2FC (scRNA)", , row_split = row_split, row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 8), cluster_rows = T, show_row_dend = F, cluster_columns = F, width = unit(0.5, "cm"), height = unit(8, "cm"))
p1+p2
dev.off()

# Choosing a consistent change of TF
genes <- c("SOX9", "MEF2C", "NFKB1", "RUNX3", "RELA", "ENO1", "CEBPB", "CEBPD")
interest.signature <- as.list(genes)
names(interest.signature) <- genes

##### survival analysis in TCGA-KIRC
source(file = "/home/longzhilin/Analysis_Code/code/analysis.diff.survival.TCGA.R")
DESeq2.normalized_counts <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/DESeq2.normalized_counts.rds")
DESeq2.normalized_counts <- log2(DESeq2.normalized_counts+1)
DESeq2.result <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/DESeq2.result.rds")
clin.data <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/clin.data.rds")
pdf("8.Immune/Macrophage/common.TFs.survival.pdf")
res <- sapply(interest.signature, function(x){
  a <- TCGA.TAM.C1QB <- analysis.diff.survival.TCGA(interest.gene = x, diff.gene.pro = DESeq2.result, exp.data.process = DESeq2.normalized_counts, clin.data = clin.data, EnhancedVolcano.plot = F, main = "TAM-C1QB", Box.plot = T, meta.signature = F, single.signature = T)
})
dev.off()

##### survival analysis in ICB data
library(ggpubr)
source(file = "/home/longzhilin/Analysis_Code/SurvivalAnalysis/Cox.function.R")
source(file = "/home/longzhilin/Analysis_Code/code/RCC.ICB.analysis.R")
normalized_expression <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/ICB.therapy/normalized_expression.rds")
patient.info.RNA <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/ICB.therapy/patient.info.RNA.rds")
#sex：F - 0; M - 1
patient.info.RNA$Sex <- gsub("^F|FEMALE|Female$", 0, patient.info.RNA$Sex)
patient.info.RNA$Sex <- gsub("^M|Male|MALE$", 1, patient.info.RNA$Sex)
patient.info.RNA$Tumor_Sample_Primary_or_Metastasis <- gsub("PRIMARY", 0, patient.info.RNA$Tumor_Sample_Primary_or_Metastasis)
patient.info.RNA$Tumor_Sample_Primary_or_Metastasis <- gsub("METASTASIS", 1, patient.info.RNA$Tumor_Sample_Primary_or_Metastasis)
pdf("8.Immune/Macrophage/common.TFs.ICB.survival.pdf")
ICB.res <- RCC.icb.analysis(signature.list = interest.signature, expresionMatrix = normalized_expression, clincal.info = patient.info.RNA)
dev.off()
