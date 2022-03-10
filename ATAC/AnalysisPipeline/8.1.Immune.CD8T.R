#' @description: analysis CD8+ T cell

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
cell.type <- "CD8+ T cell"
sub.ATAC <- subset(scATAC, subset = AnnotatedcellType == cell.type)

####----------------------------------------------------- 1. clustering cells
sub.ATAC <- RunTFIDF(sub.ATAC)
sub.ATAC <- FindTopFeatures(sub.ATAC, min.cutoff = "q1")

pdf("8.Immune/CD8T/Harmony.integration.PC20.pdf")
sub.scATAC.harmony <- Harmony.integration.reduceDimension(scATAC.object = sub.ATAC, set.resolutions = seq(0.2, 1.2, by = 0.1), group.by.vars = "dataset", assay = "ATAC", PC = 20, npcs = 30)
dev.off()
sub.scATAC.harmony$seurat_clusters <- sub.scATAC.harmony$ATAC_snn_res.0.4
Idents(sub.scATAC.harmony) <- sub.scATAC.harmony$seurat_clusters
saveRDS(sub.scATAC.harmony, file = "8.Immune/CD8T/sub.scATAC.harmony.rds")

#### integrate scRNA and scATAC
source(file = "/home/longzhilin/Analysis_Code/SingleCell/Integrate.scRNA.scATAC.R")
source(file = "/home/longzhilin/Analysis_Code/SingleCell/variableFeatureSelection.R")
sub.scRNA.harmony <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/5.Immune/CD8T/sub.scRNA.harmony.pro.rds")
DefaultAssay(sub.scRNA.harmony) <- "RNA"
VariableFeatures(sub.scRNA.harmony) <- sub.scRNA.harmony@assays$SCT@var.features
pdf("8.Immune/CD8T/Integrate.scRNA.scATAC.PC20.pdf")
transfer.anchors <- Integrate.scRNA.scATAC(scATAC.object = sub.scATAC.harmony, scRNA.object = sub.scRNA.harmony,
                                            scATAC.assay = "ACTIVITY", scRNA.assay = "RNA", ref.npcs = 30, feature.assay = "SCT", dims = 20,
                                            nfeatures = 3000, class.label = "cellType3", scATAC.reduction = "harmony", PC = 20)
Coembedding <- Coembedding.scATAC.scRNA(transfer.anchors = transfer.anchors$transfer.anchors, 
                                         scATAC.object = transfer.anchors$scATAC.object, 
                                         scATAC.reduction = "harmony", PC = 20, dims = 30, observe.vars = c("dataset", "type"),
                                         scRNA.object = sub.scRNA.harmony, class.label = "cellType3")
dev.off()

source(file = "/home/longzhilin/Analysis_Code/code/ratio.plot.R")
sub.scATAC.harmony$scRNA.id <- Coembedding$scATAC.object$predicted.id
sub.scATAC.harmony$scATAC.id <- sub.scATAC.harmony$seurat_clusters
pdf("8.Immune/CD8T/scRNA.scATAC.comparsion.pdf")
ratio.plot(seurat.object = sub.scATAC.harmony, id.vars1 = "scRNA.id", id.vars2 = "scATAC.id", angle = 60)
dev.off()

####----------------------------------------------------- 2. differential peaks
Idents(sub.scATAC.harmony) <- sub.scATAC.harmony$seurat_clusters
idents <- as.character(levels(sub.scATAC.harmony))
cluster.DARs <- FindAllMarkers(sub.scATAC.harmony, 
                                test.use = 'LR',
                                logfc.threshold=0.25, 
                                min.pct = 0.05, # often necessary to lower the min.pct threshold
                                latent.vars = "peak_region_fragments")
cf <- ClosestFeature(sub.scATAC.harmony, regions = rownames(cluster.DARs)) # Find the closest feature to a given set of genomic regions
cluster.DARs <- cbind(cluster.DARs, gene=cf$gene_name, gene_biotype = cf$gene_biotype, type = cf$type, distance=cf$distance)
colnames(cluster.DARs)[6:7] <- c("cluster", "genomicRegion")
saveFormat <- lapply(idents, function(x){
  index <- which(cluster.DARs$cluster == x)
  DARs <- cluster.DARs[index,]
  DARs.up <- DARs %>% filter(avg_log2FC>0) %>% arrange(p_val_adj)
  DARs.down <- DARs %>% filter(avg_log2FC<0) %>% arrange(desc(p_val_adj))
  DARs <- rbind(DARs.up, DARs.down)
  return(DARs)
})
write.xlsx(saveFormat, file = "8.Immune/CD8T/cluster.all.DARs.xlsx", sheetName = idents, rowNames = F)
saveRDS(cluster.DARs, file = "8.Immune/CD8T/cluster.all.DARs.rds")

####----------------------------------------------------- 2. annotated cell type
DefaultAssay(sub.scATAC.harmony) <- "ACTIVITY"
features <- c("CD3E", "CD3D", "CD8A", "APOE", "C1QC", "PLVAP", "CD4", "FGFBP2", "KLRD1")
pdf("8.Immune/CD8T/markerExpression.pdf")
VlnPlot(sub.scATAC.harmony, features = c("CD3D", "CD3E", "C1QC", "APOE"), ncol = 2)
FeaturePlot(sub.scATAC.harmony, features = c("CD3D", "CD3E", "C1QC", "APOE"), cols = c("lightgrey", "red"), ncol = 2)
VlnPlot(sub.scATAC.harmony, features = c("CD8A", "CD8B", "PLVAP", "ESM1"), ncol = 2)
FeaturePlot(sub.scATAC.harmony, features = c("CD8A", "CD8B", "PLVAP", "ESM1"), cols = c("lightgrey", "red"), ncol = 2)
VlnPlot(sub.scATAC.harmony, features = c("PRF1", "TCF7", "TOX", "PDCD1"), ncol = 2)
FeaturePlot(sub.scATAC.harmony, features = c("PRF1", "TCF7", "TOX", "PDCD1"), cols = c("lightgrey", "red"), ncol = 2)
dev.off()
pdf("8.Immune/CD8T/cluster.markerExpression.dot.pdf", height = 3, width = 5)
DotPlot(sub.scATAC.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "seurat_clusters", dot.scale = 3.5) + RotatedAxis() + theme(axis.text.x = element_text(size = 7)) + theme(axis.text.y = element_text(size = 8))
DotPlot(sub.scATAC.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "seurat_clusters", dot.scale = 3.5) + RotatedAxis() + NoLegend() + theme(axis.text.x = element_text(size = 7)) + theme(axis.text.y = element_text(size = 8))
dev.off()

pdf("8.Immune/CD8T/cluster.markerExpression.coveragePlot.pdf")
DefaultAssay(sub.scATAC.harmony) <- "Peaks"
CoveragePlot(
  object = sub.scATAC.harmony,
  region = c("CD3D", "CD3E", "APOE", "C1QC"),
  extend.upstream = 3000,
  extend.downstream = 3000
)
CoveragePlot(
  object = sub.scATAC.harmony,
  region = c("CD8A", "CD8B", "PLVAP", "ESM1"),
  extend.upstream = 3000,
  extend.downstream = 3000
)
region.highlight <- c("chr4-108157400-108174201", "chr5-35851000-35865000", "chr17-40559800-40566601", "chr1-169706600-169712401")
CoveragePlot(
  object = sub.scATAC.harmony,
  region = region.highlight,
  extend.upstream = 3000,
  extend.downstream = 3000
)
region.highlight <- c("chr8-59115000-59121000", "chr2-241856401-241860601", "chr10-95755400-95758401", "chr12-68157600-68159800")
CoveragePlot(
  object = sub.scATAC.harmony,
  region = region.highlight,
  extend.upstream = 3000,
  extend.downstream = 3000
)
dev.off()

##Plot--- cell number
pdf("8.Immune/CD8T/cluster.ratio.pdf", height = 4, width = 6)
ratio.plot(seurat.object = sub.scATAC.harmony, id.vars1 = "orig.ident", id.vars2 = "seurat_clusters", angle = 60)
dev.off()

####----------------------------------------------------- 3. remove the condfounding cell type
# 1. cell number < 100
# 2. condfounding cell type
sub.scATAC.harmony <- subset(sub.scATAC.harmony, subset = seurat_clusters %in% c(0:4))
sub.scATAC.harmony$seurat_clusters <- as.character(sub.scATAC.harmony$seurat_clusters)
#consistent with scRNA-seq
cellType2 <- sub.scATAC.harmony$seurat_clusters
cellType2 <- gsub("^0$", "C3", cellType2)
cellType2 <- gsub("^1$", "C2", cellType2)
cellType2 <- gsub("^2$", "C1", cellType2)
cellType2 <- gsub("^3$", "C4", cellType2)
cellType2 <- gsub("^4$", "C5", cellType2)
sub.scATAC.harmony$cellType2 <- factor(sub.scATAC.harmony$cellType2, levels = c("C1", "C2", "C3", "C4", "C5"))
Idents(sub.scATAC.harmony) <- sub.scATAC.harmony$cellType2

cellType3 <- sub.scATAC.harmony$cellType2
cellType3 <- gsub("^C1$", "Unknown", cellType3)
cellType3 <- gsub("^C2$", "Tissue-resident", cellType3)
cellType3 <- gsub("^C3$", "Exhaustion", cellType3)
cellType3 <- gsub("^C4$", "Exhausted IEG.C2", cellType3)
cellType3 <- gsub("^C5$", "Exhausted IEG.C1", cellType3)
sub.scATAC.harmony$cellType3 <- factor(cellType3, levels = c("Unknown", "Tissue-resident", "Exhausted IEG.C1", "Exhausted IEG.C2", "Exhaustion"))
Idents(sub.scATAC.harmony) <- sub.scATAC.harmony$cellType3
saveRDS(sub.scATAC.harmony, file = "8.Immune/CD8T/sub.scATAC.harmony.pro.rds")

#### set idents
sub.scATAC.harmony$cellType <- sub.scATAC.harmony$cellType3
Idents(sub.scATAC.harmony) <- sub.scATAC.harmony$cellType
cellType.colors <- c("#00BF7D", "#A3A500", "#00BFC4", "#00B0F6", "#C77CFF")

pdf("8.Immune/CD8T/cluster.pro.pdf")
DimPlot(object = sub.scATAC.harmony, label = TRUE, group.by = "seurat_clusters", pt.size = 1.5) + NoLegend()
DimPlot(object = sub.scATAC.harmony, label = TRUE, group.by = "orig.ident",  pt.size = 1.5) + NoLegend()
DimPlot(object = sub.scATAC.harmony, label = TRUE, group.by = "seurat_clusters", pt.size = 1.5, reduction = "tsne") + NoLegend()
DimPlot(object = sub.scATAC.harmony, label = TRUE, group.by = "orig.ident",  pt.size = 1.5, reduction = "tsne") + NoLegend()
dev.off()

pdf("8.Immune/CD8T/cellType.pdf")
DimPlot(object = sub.scATAC.harmony, cols = cellType.colors, label = TRUE, group.by = "cellType", pt.size = 1.5) + NoLegend()
DimPlot(object = sub.scATAC.harmony, label = TRUE, group.by = "orig.ident",  pt.size = 1.5) 
DimPlot(object = sub.scATAC.harmony, cols = cellType.colors, label = TRUE, group.by = "cellType", pt.size = 1.5, reduction = "tsne") + NoLegend()
DimPlot(object = sub.scATAC.harmony, label = TRUE, group.by = "orig.ident",  pt.size = 1.5, reduction = "tsne")
dev.off()

##Plot--- cell number
pdf("8.Immune/CD8T/cellType.ratio.pdf", height = 4, width = 6)
ratio.plot(seurat.object = sub.scATAC.harmony, id.vars1 = "orig.ident", id.vars2 = "cellType", angle = 60)
dev.off()
####----------------------------------------------------- 4. marker activity
DefaultAssay(sub.scATAC.harmony) <- "ACTIVITY"
features <- c("LEF1", "IL7R", "CCR7", "SELL", "TCF7", "TGFB1", "CD44", "CD69", "ZNF683", "ITGAE", "ITGA1",  
              "TNF", "IFNG", "KLRG1", "GZMA", "GZMH",
              "NR4A1", "JUNB", "FOS", "ATF3", "DNAJB1", "HSPA1A", "EOMES", 
              "GZMK", "GZMB", "PRF1", "TNFRSF9", "TOX", "ENTPD1", "PDCD1", "CTLA4", "TIGIT", "LAG3", "HAVCR2")
pdf("8.Immune/CD8T/markerExpression.pro.pdf")
VlnPlot(sub.scATAC.harmony, features = c("LEF1", "IL7R", "TCF7", "CD44"), ncol = 2)
FeaturePlot(sub.scATAC.harmony, features = c("LEF1", "IL7R", "TCF7", "CD44"), cols = c("lightgrey", "red"), ncol = 2)
VlnPlot(sub.scATAC.harmony, features = c("CD69", "ZNF683", "ITGA1", "ITGAE"), ncol = 2)
FeaturePlot(sub.scATAC.harmony, features = c("CD69", "ZNF683", "ITGA1", "ITGAE"), cols = c("lightgrey", "red"), ncol = 2)
VlnPlot(sub.scATAC.harmony, features = c("JUNB", "FOS", "DNAJB1", "HSPA1A"), ncol = 2)
FeaturePlot(sub.scATAC.harmony, features = c("JUNB", "FOS", "DNAJB1", "HSPA1A"), cols = c("lightgrey", "red"), ncol = 2)
VlnPlot(sub.scATAC.harmony, features = c("PDCD1", "TOX", "ENTPD1", "HAVCR2"), ncol = 2)
FeaturePlot(sub.scATAC.harmony, features = c("PDCD1", "TOX", "ENTPD1", "HAVCR2"), cols = c("lightgrey", "red"), ncol = 2)
dev.off()

# vlnplot
pdf("8.Immune/CD8T/markerExpression.vlnplot.pdf", width = 5)
p1 <- VlnPlot(sub.scATAC.harmony,features = features[1:17], group.by="cellType", same.y.lims=T,flip = T, stack = T) & xlab("") & ylab("Log-normalized activity") & theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.5))
p2 <- VlnPlot(sub.scATAC.harmony,features = features[18:34], group.by="cellType", same.y.lims=T,flip = T, stack = T) & xlab("") & ylab("Log-normalized activity") & theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.5))
ggarrange(p1,p2,ncol=2)
dev.off()

avg.expression <- AverageExpression(sub.scATAC.harmony, features, assays = "ACTIVITY", slot = "data")
avg.expression <- scale(t(avg.expression$ACTIVITY))
pdf("8.Immune/CD8T/cellType.markerExpression.heatmap.pdf")
Heatmap(t(avg.expression), cluster_rows = F, show_column_dend = F, name = "Gene activity",
        width = unit(2.5, "cm"), height = unit(8, "cm"),
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
dev.off()

source("/home/longzhilin/Analysis_Code/SingleCell/FindRegion.R")
DefaultAssay(sub.scATAC.harmony) <- "Peaks"
peak.info <- readRDS("4.Peak/peak.annotation.simple.ChIPseeker.rds")
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
pdf("8.Immune/CD8T/markerExpression.pro.peaks.pdf")
res <- sapply(features, function(x){
  cat(x, "...\n")
  regions <- FindRegion(object = sub.scATAC.harmony, region = x, assay = "Peaks", extend.upstream = 1000, extend.downstream = 1000)
  idx <- data.frame(findOverlaps(regions, promoter))
  if(nrow(idx)>0){
    p <- CoveragePlot(
      object = sub.scATAC.harmony,
      region = x,
      links = F,
      peaks = F,
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
      peaks = F)
  }
  print(p)
  return(regions)
})
dev.off()

Naive.Marker <- c("LEF1", "IL7R", "CCR7", "SELL")
pdf("8.Immune/CD8T/markerExpression.Naive.peaks.pdf")
res <- lapply(Naive.Marker, function(x){
  cat(x, "...\n")
  regions <- FindRegion(object = sub.scATAC.harmony, region = x, assay = "Peaks", extend.upstream = 3000, extend.downstream = 3000)
  idx <- data.frame(findOverlaps(regions, promoter))
  if(nrow(idx)>0){
    p <- CoveragePlot(
      object = sub.scATAC.harmony,
      region = x,
      links = F,
      peaks = F,
      extend.upstream = 0,
      extend.downstream = 0,
      region.highlight = promoter[idx[,2],])
  }else{
    p <- CoveragePlot(
      object = sub.scATAC.harmony,
      region = x,
      links = F,
      extend.upstream = 0,
      extend.downstream = 0,
      peaks = F)
  }
  return(p)
})
ggarrange(res[[1]], res[[2]], res[[3]], res[[4]], ncol = 2, nrow = 2)
dev.off()
pdf("8.Immune/CD8T/markerExpression.Naive.maker.pdf", height = 2)
DotPlot(sub.scATAC.harmony, features = c(Naive.Marker, "IL2RA", "CD44", "CD69", "TOX", "ENTPD1", "PDCD1", "CTLA4", "TIGIT", "LAG3"), cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1), axis.text.y = element_text(size = 8))
dev.off()

DefaultAssay(sub.scATAC.harmony) <- "ACTIVITY"
pdf("8.Immune/CD8T/cellType.marker.dot.pdf", height = 2)
DotPlot(sub.scATAC.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1)) + theme(axis.text.y = element_text(size = 8))
DotPlot(sub.scATAC.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + NoLegend() + theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1)) + theme(axis.text.y = element_text(size = 8))
dev.off()
pdf("8.Immune/CD8T/cellType.marker.dot.test.pdf", height = 4)
DotPlot(sub.scATAC.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1)) + theme(axis.text.y = element_text(size = 8))
DotPlot(sub.scATAC.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + NoLegend() + theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1)) + theme(axis.text.y = element_text(size = 8))
dev.off()

TFs <- c("FOXP1", "TBX21", "PRDM1", "TCF7", "LEF1", "RUNX3", "GNLY", "GZMA", "GZMB", "GZMK", "GZMM", "PRF1", "PDCD1", "TOX")
pdf("8.Immune/CD8T/cellType.marker.TF.dot.pdf", height = 2)
DotPlot(sub.scATAC.harmony, features = TFs, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1)) + theme(axis.text.y = element_text(size = 8))
DotPlot(sub.scATAC.harmony, features = TFs, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + NoLegend() + theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1)) + theme(axis.text.y = element_text(size = 8))
dev.off()

####----------------------------------------------------- 3. analysis cell state based on the functional gene sets
## 1.load signature
DefaultAssay(sub.scATAC.harmony) <- "ACTIVITY"
T.signature1 <- read.table("/data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/Tsignature.Braun.CancerCell.txt", header = T, stringsAsFactors = F, sep = "\t")
T.signature2 <- read.table("/data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/CD8.Exhausted(anti-pd1).Bi.2021.CancerCell&Sade-Feldman et al. 2018.Cell.txt", header = T, stringsAsFactors = F, sep = "\t")
T.signature3 <- read.table("/data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/FunctionalStateOfTcell.Mathewson.2021.Cell.txt", header = T, stringsAsFactors = F, sep = "\t")
T.signature4 <- read.table("/data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/ImmuneSignatureGeneSet.T.Chung.2017.NatureComm.txt", header = T, stringsAsFactors = F, sep = "\t")
T.signature5 <- read.table("/data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/CD8.Tstate.vanderLeun.2020.NatureReviewsCancer.txt", header = T, stringsAsFactors = F, sep = "\t")
T.signature6 <- read.table("/data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/CD8.Resident&Exhausted.MomenehForoutan.2020.BioRXiv.txt", header = T, stringsAsFactors = F, sep = "\t")
signature.genes <- rbind(T.signature1, T.signature2)
signature.genes <- rbind(signature.genes, T.signature3)
signature.genes <- rbind(signature.genes, T.signature4)
signature.genes <- rbind(signature.genes, T.signature5)
signature.genes <- rbind(signature.genes, T.signature6)
signature.genes <- signature.genes[which(signature.genes$Type %in% c("Naive", "Effector memory", "Central memory", "Res", "Exh",
                                                                     "Cell_stress", "Cytotoxicity Signature", "Exhaustion",
                                                                     "Terminal_differentiation", "Progenitor Exhausted CD8", "Terminally Exhausted CD8")),]
signature.genes$Type <- gsub("Cytotoxicity Signature", "Cytotoxic", signature.genes$Type)
signature.genes$Type <- gsub("_", " ", signature.genes$Type)
signature.genes$Type <- gsub(" CD8", "", signature.genes$Type)
signature.genes$Type <- gsub("Exhausted", "Exhaustion", signature.genes$Type)
signature.genes$Type <- gsub("Terminally", "Terminal", signature.genes$Type)

##-- 2.VISION method
source(file = "/home/longzhilin/Analysis_Code/SingleCell/vision_seurat.R")
source(file = "/home/longzhilin/Analysis_Code/SingleCell/vision.plot.R")
require(VISION)
library(reshape2)
options(mc.cores = 36)
## construct the vision calculate mode
# https://yoseflab.github.io/VISION/articles/web_only/Seurat.html
# Expression data should be scaled and normalized, but not log-transformed
# obj@[[assay]]@counts is used as the expression input (after normalizing to a library size of 10,000)
exp.data <- GetAssayData(sub.scATAC.harmony, slot = "data", assay = "ACTIVITY")
idx <- match(signature.genes$Gene, rownames(exp.data))
signature.genes <- signature.genes[which(!is.na(idx)),]
signature.names <- unique(signature.genes$Type)
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
saveRDS(vision.obj, file = "8.Immune/CD8T/vision_signatureScore.rds")

library(reshape2)
sigScore <- as.data.frame(vision.obj@SigScores)
sigScore$group <- sub.scATAC.harmony$cellType
vision.group.score <- melt(sigScore, id.vars = "group", variable.name = "Type", value.name = "Score")
pdf("8.Immune/CD8T/vision_signatureScore.pdf")
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
    my_comparisons <- as.list(as.data.frame(combn(levels(sub.scATAC.harmony$cellType)[-1],2)))
    a <- sigScore[, c(x, "group")]
    names(a) <- c("Signature score", "group")
    p <- ggviolin(a, x = "group", y = "Signature score", title = x, color = "black", alpha = 0.8, fill = "group", add = "boxplot", add.params = list(fill = "white", size = 0.05)) + 
         xlab("") + rotate_x_text(angle = 45, vjust = 1) + scale_fill_manual(values = cellType.colors) + 
         stat_compare_means(comparisons = my_comparisons, label = "p.signif") + NoLegend()
    return(p)
})
p <- ggarrange(plotlist = plot.list[1:4],ncol = 2, nrow = 2)
print(p)
p <- ggarrange(plotlist = plot.list[5:8],ncol = 2, nrow = 2)
print(p)
p <- ggarrange(plotlist = plot.list[9:11],ncol = 2, nrow = 2)
print(p)

# comparison between groups
vision.score.scale <- as.data.frame(scale(vision.score[, c("Progenitor Exhaustion", "Terminal Exhaustion")]))
vision.score.scale$group <- sub.scATAC.harmony$cellType3
group.signature.data <- pivot_longer(vision.score.scale[,c("group", "Progenitor Exhaustion", "Terminal Exhaustion")], cols = 2:3, names_to = "Type")
p <- ggviolin(group.signature.data, x = "group", y = "value",
          color = "Type", palette = c("#00A087FF", "#F39B7FFF"), fill = "white",
          add = "jitter", add.params = list(size = 0.1)) + ylab("Signature score") + xlab("")
p <- p + rotate_x_text(angle = 45, vjust = 1) + stat_compare_means(aes(group = Type), label = "p.signif")
print(p)
dev.off()

####----------------------------------------------------- 4. peak and motif analysis
DefaultAssay(sub.scATAC.harmony) <- "Peaks"
# add motif information
sub.scATAC.harmony <- AddMotifs(
  object = sub.scATAC.harmony,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = human_pwms_v2
)

####---- differential Motifs
library(openxlsx)
library(chromVARmotifs)
library(BSgenome.Hsapiens.UCSC.hg38)
data("human_pwms_v2")
Idents(sub.scATAC.harmony) <- sub.scATAC.harmony$cellType
scATAC.harmony <- RunChromVAR(
  object = sub.scATAC.harmony,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
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
  up.x <- x %>% filter(avg_log2FC>0) %>% arrange(desc(avg_log2FC))
  down.x <- x %>% filter(avg_log2FC<=0) %>% arrange(avg_log2FC)
  x <- rbind(up.x, down.x)
  return(x)
})
names(sub.motifs.chromVAR) <- idents
write.xlsx(sub.motifs.chromVAR, file = "8.Immune/CD8T/sub.motifs.chromVAR.human_pwms_v2.xlsx", sheetName = idents, rowNames = T, overwrite=T)
saveRDS(sub.motifs.chromVAR, file = "8.Immune/CD8T/sub.motifs.chromVAR.human_pwms_v2.rds")

sub.motifs.chromVAR.sig <- lapply(sub.motifs.chromVAR, function(x){
  x <- x %>% filter(avg_log2FC>0.5 & p_val_adj<0.05) %>% arrange(desc(avg_log2FC))
  return(x)
})
write.xlsx(sub.motifs.chromVAR.sig, file = "8.Immune/CD8T/sub.motifs.chromVAR.sig.human_pwms_v2.xlsx", sheetName = idents, rowNames = T, overwrite=T)
####-- plot chromVAR deviation for this TF:
top.TFs <- lapply(names(sub.motifs.chromVAR.sig), function(x){
    sub.motifs.chromVAR.sig[[x]]$Type <- rep(x, nrow(sub.motifs.chromVAR.sig[[x]]))
    return(sub.motifs.chromVAR.sig[[x]]$gene[1:5])
})
top.TFs <- unique(unlist(top.TFs))
top.TFs[11:15] <- c("STAT3", "NFATC3", "NR4A1", "NFATC2", "NFATC1") 
top.TFs[21:25] <- c("RELA", "REL", "NFKB1", "BATF", "EOMES") 

enrich <- lapply(sub.motifs.chromVAR.sig, function(x){
    idx <- which(x$gene %in% top.TFs)
    return(x[idx, c("gene", "avg_log2FC", "p_val_adj")])
})                    
names(enrich) <- names(sub.motifs.chromVAR.sig)
enrich <- lapply(names(enrich), function(x){
    enrich[[x]]$Type <- rep(x, nrow(enrich[[x]]))
    return(enrich[[x]])
})  
top.TFs.res <- Reduce(function(x,y) rbind(x,y), enrich)
rownames(top.TFs.res) <- NULL

#排序
top.TFs.res$p.adjust <- -log10(top.TFs.res$p_val_adj)
top.TFs.res$Type <- factor(top.TFs.res$Type, levels = c("Naive", "Tissue-resident", "Exhausted IEG.C1", "Exhausted IEG.C2", "Exhausted"))
top.TFs.res$gene <- factor(top.TFs.res$gene , levels = unique(top.TFs.res$gene))
top.TFs.res$avg_log2FC <- top.TFs.res$avg_log2FC
library(ggthemes)
pdf("8.Immune/CD8T/top.TF.plot.pdf", width = 3.5)
p1 <- ggplot(top.TFs.res, 
            aes(x = Type, 
                y = gene, 
                size = avg_log2FC,
                color = p.adjust,
                fill = p.adjust))
p2 <- p1 + guides(color=FALSE) + geom_point(shape = 21) + theme_bw() + scale_size("avg_log2FC", breaks = c(0, 1, 2, 3), labels = c(0, 2, 4, 6), range = c(1.5, 5.5)) + scale_color_gradientn(colours = c("#fed71b", "#fcc41c", "#f27620", "#e92925"), breaks = c(10, 30, 50, 70), labels = c(1e-10, 1e-30, 1e-50, 1e-70)) + scale_fill_gradientn(colours = c("#fed71b", "#fcc41c", "#f27620", "#e92925"), breaks = c(10, 30, 50, 70), labels = c(1e-10, 1e-30, 1e-50, 1e-70))
p2 <- p2 + theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1, size = 8), axis.text.y = element_text(size = 8)) + xlab("") + ylab("")
print(p2)

#heatmap
pathway.p <- top.TFs.res[, c("gene", "Type", "avg_log2FC")]
pathway.data <- pivot_wider(pathway.p, names_from = "Type", values_from = "avg_log2FC")
pathway.data <- as.data.frame(pathway.data)
rownames(pathway.data) <- pathway.data$gene
pathway.data <- as.data.frame(pathway.data[,-1])
cols <- colorRamp2(c(0, 2, 4), c("white", "#ffff66", "red"))
Heatmap(pathway.data, name = "avg_log2FC", show_column_dend = F, show_row_dend = F, col = cols,  
        na_col = "#f2f3f4", border = "grey", border_gp = gpar(col = "grey"), rect_gp = gpar(col = "grey"),
        width = unit(3, "cm"), height = unit(10, "cm"), cluster_rows = F, cluster_columns = F,
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
dev.off()

####Plot --- cellType specific motif heatmap
top5 <- sapply(names(sub.motifs.chromVAR.sig), function(x){
    sub.motifs.chromVAR.sig[[x]]$Type <- rep(x, nrow(sub.motifs.chromVAR.sig[[x]]))
    return(sub.motifs.chromVAR.sig[[x]]$gene[1:5])
})
names(top5) <- names(sub.motifs.chromVAR.sig)
top5 <- unique(as.character(top5))
sig.motifs <- sub.motifs.chromVAR.sig %>% 
                                bind_rows %>%
                                select(gene) %>%
                                distinct()
idx <- match(sig.motifs$gene, motif.info$TF)
sig.motifNames <- motif.info$originName[idx]

motifs.avgExp <- AverageExpression(sub.scATAC.harmony, features = sig.motifNames, assays = "chromvar")
motifs.avgExp <- motifs.avgExp$chromvar
rownames(motifs.avgExp) <- sig.motifs$gene
zScore <- function(x){(x - mean(x)) /sd(x)}
motifs.avgExp.scale <- apply(motifs.avgExp, 1, zScore) %>% t() # row: TF; column: cell type
interest.genes <- c(top5, "NR4A1", "NR4A2", "NR4A3", "NFAT5", "NFATC3", "NFATC4", "NFATC2", "NFATC1", "FOXO1", "PRDM1", "EOMES")
mark.idx <- match(interest.genes, rownames(motifs.avgExp.scale))
pdf("8.Immune/CD8T/chromVAR.cellType.sig.motifs.heatmap.pdf")
Heatmap(motifs.avgExp.scale, name = "Deviation score", 
        width = unit(5, "cm"), height = unit(12, "cm"),
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
        show_row_dend = F, show_column_dend = F, show_column_names = T, show_row_names = F, 
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 8), by_row = T)) + 
rowAnnotation(link = anno_mark(at = mark.idx, labels = interest.genes, link_width = unit(2, "mm"), labels_gp = gpar(fontsize = 5), padding = unit(1, "mm")))
dev.off()

# combined the DEGs
sub.scRNA.DEGs <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/5.Immune/CD8T/cellType.DEGs.rds")
common.DEGs <- sapply(names(sub.motifs.chromVAR.sig), function(x){
  TFs <- sub.motifs.chromVAR.sig[[x]]$gene
  name <- gsub("\\..*", "", x)
  idx <- grep(name, names(sub.scRNA.DEGs))
  if(length(idx)>0){
    inter.genes <- sapply(idx, function(y){
      a <- intersect(sub.scRNA.DEGs[[y]]$gene, TFs)
      return(a)
    })
  }else{
    return(NULL)
  }
})
