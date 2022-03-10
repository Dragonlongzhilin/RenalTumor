#' @description: analysis CD8+ T population

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
setwd(dir = "/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA")
source(file = "/home/longzhilin/Analysis_Code/Visualization/colorPalettes.R")
source(file = "/home/longzhilin/Analysis_Code/SingleCell/scRNA.Integrate.multipleSample.R")
source(file = "/home/longzhilin/Analysis_Code/SingleCell/variableFeatureSelection.R")
source(file = "/home/longzhilin/Analysis_Code/code/ratio.plot.R")

data.merge <- readRDS("data.merge.pro.rds")
cell.Type <- "CD8+ T cell"
sub.scRNA <- subset(data.merge, subset = cellType_low == cell.Type)
DefaultAssay(sub.scRNA) <- "RNA"
sub.scRNA <- DietSeurat(sub.scRNA, assays = "RNA")
# observe the batch effect
# sub.scRNA <- SCTransform(sub.scRNA, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
# sub.scRNA <- RunPCA(sub.scRNA, npcs = 30, verbose = FALSE) %>%
# RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
# FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
# FindClusters(resolution = seq(0.2, 1.2, by = 0.1), verbose = FALSE)

# split the dataset into a list of two seurat objects (stim and CTRL)
sub.list <- SplitObject(sub.scRNA, split.by = "orig.ident")
sub.list <- lapply(X = sub.list, FUN = function(x) {
    x <- SCTransform(x, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
})
sub.scRNA <- merge(sub.list[[1]], y = sub.list[2:length(sub.list)], project = "Renal_CD8")
DefaultAssay(sub.scRNA) <- "SCT"
seurat.features.SCT <- SelectIntegrationFeatures(object.list = sub.list, nfeatures = 3000)
VariableFeatures(sub.scRNA) <- seurat.features.SCT
pdf("5.Immune/CD8T/SCT.Harmony.Integration.PC20.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "SCT", set.resolutions = seq(0.2, 1.2, by = 0.1), PC = 20, npcs = 30)
dev.off()
sub.scRNA.harmony$seurat_clusters <- sub.scRNA.harmony$SCT_snn_res.0.5
Idents(sub.scRNA.harmony) <- sub.scRNA.harmony$seurat_clusters
DefaultAssay(sub.scRNA.harmony) <- "RNA"
sub.scRNA.harmony <- NormalizeData(sub.scRNA.harmony, verbose = FALSE)
sub.scRNA.harmony <- ScaleData(sub.scRNA.harmony, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mt"), features = rownames(sub.scRNA.harmony))
saveRDS(sub.scRNA.harmony, file = "5.Immune/CD8T/sub.scRNA.harmony.rds")

####-------------------------------------------------1.the characteristics of each cluster
#### 1.marker 
features <- c("CD3E", "CD3D", "CD8A", "APOE", "C1QC", "PLVAP", "CD4", "FGFBP2", "KLRD1")
pdf("5.Immune/CD8T/markerExpression.pdf")              
VlnPlot(sub.scRNA.harmony, features = c("CD3D", "CD3E", "C1QC", "APOE"), ncol = 2)
FeaturePlot(sub.scRNA.harmony, features = c("CD3D", "CD3E", "C1QC", "APOE"), cols = c("lightgrey", "red"), ncol = 2)
VlnPlot(sub.scRNA.harmony, features = c("CD8A", "CD8B", "PLVAP", "ESM1"), ncol = 2)
FeaturePlot(sub.scRNA.harmony, features = c("CD8A", "CD8B", "PLVAP", "ESM1"), cols = c("lightgrey", "red"), ncol = 2)
VlnPlot(sub.scRNA.harmony, features = c("PRF1", "IFNG", "TOX", "PDCD1"), ncol = 2)
FeaturePlot(sub.scRNA.harmony, features = c("PRF1", "IFNG", "TOX", "PDCD1"), cols = c("lightgrey", "red"), ncol = 2)
dev.off()

pdf("5.Immune/CD8T/cluster.markerExpression.dot.pdf", height = 3, width = 5)
DotPlot(sub.scRNA.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "seurat_clusters", dot.scale = 3.5) + xlab("") + ylab("Cluser") + RotatedAxis() + theme(axis.text.x = element_text(size = 8)) + theme(axis.text.y = element_text(size = 8))
DotPlot(sub.scRNA.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "seurat_clusters", dot.scale = 3.5) + xlab("") + ylab("Cluser") + RotatedAxis() + NoLegend() + theme(axis.text.x = element_text(size = 8)) + theme(axis.text.y = element_text(size = 8))
dev.off()

pdf("5.Immune/CD8T/cluster.pdf")
DimPlot(object = sub.scRNA.harmony, reduction = 'umap', pt.size = 1.5, label = TRUE, group.by = "orig.ident")
DimPlot(object = sub.scRNA.harmony, reduction = 'tsne', pt.size = 1.5, label = TRUE, group.by = "orig.ident")
DimPlot(object = sub.scRNA.harmony, reduction = 'tsne', pt.size = 1.5, label = TRUE, group.by = "seurat_clusters")+NoLegend()
DimPlot(object = sub.scRNA.harmony, reduction = 'umap', pt.size = 1.5, label = TRUE, group.by = "seurat_clusters")+NoLegend()
dev.off()

##Plot--- cell number
pdf("5.Immune/CD8T/cluster.ratio.pdf", height = 4, width = 6)
ratio.plot(seurat.object = sub.scRNA.harmony, id.vars1 = "orig.ident", id.vars2 = "seurat_clusters", angle = 60)
dev.off()

####-------------------------------------------------2.Filter out mixed cell populations
sub.scRNA.harmony <- subset(sub.scRNA.harmony, subset = seurat_clusters %in% 0:3) # 4078 cells
sub.scRNA.harmony$seurat_clusters <- as.character(sub.scRNA.harmony$seurat_clusters)
#定义cellType
sub.scRNA.harmony$cellType2 <- paste0("C", as.numeric(sub.scRNA.harmony$seurat_clusters)+1)
sub.scRNA.harmony$cellType2 <- factor(sub.scRNA.harmony$cellType2, levels = c("C1", "C2", "C3", "C4"))

cellType3 <- sub.scRNA.harmony$cellType2
cellType3 <- gsub("^C1$", "Exhausted IEG", cellType3)
cellType3 <- gsub("^C2$", "Tissue-resident.C1", cellType3)
cellType3 <- gsub("^C3$", "Exhausted", cellType3)
cellType3 <- gsub("^C4$", "Tissue-resident.C2", cellType3)
sub.scRNA.harmony$cellType3 <- factor(cellType3, levels = c("Tissue-resident.C1", "Tissue-resident.C2", "Exhausted IEG", "Exhausted"))
Idents(sub.scRNA.harmony) <- sub.scRNA.harmony$cellType3
saveRDS(sub.scRNA.harmony, file = "5.Immune/CD8T/sub.scRNA.harmony.pro.rds")

cellType3 <- sub.scRNA.harmony$cellType2
cellType3 <- gsub("^C2$", "Effector memory", cellType3)
cellType3 <- gsub("^C1$", "Exhausted", cellType3)
cellType3 <- gsub("^C3$", "Terminally exhausted", cellType3)
cellType3 <- gsub("^C4$", "Exhausted", cellType3)
sub.scRNA.harmony$cellType3 <- factor(cellType3, levels = c("Effector memory", "Exhausted", "Terminally exhausted"))

#### set Idents
sub.scRNA.harmony$cellType <- sub.scRNA.harmony$cellType3
Idents(sub.scRNA.harmony) <- sub.scRNA.harmony$cellType
cellType.colors <- c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")

DefaultAssay(sub.scRNA.harmony) <- "RNA"
pdf("5.Immune/CD8T/cluster.pro.pdf")
DimPlot(object = sub.scRNA.harmony, reduction = 'umap', pt.size = 1.5, label = TRUE, group.by = "orig.ident")
DimPlot(object = sub.scRNA.harmony, reduction = 'tsne', pt.size = 1.5, label = TRUE, group.by = "orig.ident")
DimPlot(object = sub.scRNA.harmony, reduction = 'tsne', pt.size = 1.5, label.size = 6, label = TRUE, group.by = "seurat_clusters")+NoLegend()
DimPlot(object = sub.scRNA.harmony, reduction = 'umap', pt.size = 1.5, label.size = 6, label = TRUE, group.by = "seurat_clusters")+NoLegend()
dev.off()

DefaultAssay(sub.scRNA.harmony) <- "RNA"
pdf("5.Immune/CD8T/cellType.pdf")
DimPlot(object = sub.scRNA.harmony, reduction = 'umap', pt.size = 1.5, label.size = 6, label = TRUE, group.by = "cellType")+NoLegend()
DimPlot(object = sub.scRNA.harmony, reduction = 'umap', pt.size = 1.5, label.size = 6, label = TRUE, group.by = "cellType")
dev.off()

pdf("5.Immune/CD8T/cellType.ratio.pdf", height = 4, width = 6)
ratio.plot(seurat.object = sub.scRNA.harmony, id.vars1 = "orig.ident", id.vars2 = "cellType", angle = 60)
ratio.plot(seurat.object = sub.scRNA.harmony, id.vars1 = "orig.ident", id.vars2 = "cellType", angle = 60)
dev.off()

#### singleR
library(SingleR)
ImmGen.se <- ImmGenData()
DatabaseImmuneCell.se <- DatabaseImmuneCellExpressionData()
NovershternHematopoietic.se <- NovershternHematopoieticData()
MonacoImmune.se <- MonacoImmuneData()
singler.ImmGen.predict.cluster <- SingleR(test = sub.scRNA.harmony@assays$RNA@data, 
                                        ref = ImmGen.se, 
                                        labels = ImmGen.se$label.fine,
                                        clusters = as.character(sub.scRNA.harmony@meta.data$seurat_clusters))
singler.ImmGen.predict.cluster$labels <- rownames(singler.ImmGen.predict.cluster)                         
singler.DatabaseImmuneCell.predict.cluster <- SingleR(test = sub.scRNA.harmony@assays$RNA@data, 
                                      ref = DatabaseImmuneCell.se,
                                      labels = DatabaseImmuneCell.se$label.fine,
                                      clusters = as.character(sub.scRNA.harmony@meta.data$seurat_clusters))
singler.DatabaseImmuneCell.predict.cluster$labels <- rownames(singler.DatabaseImmuneCell.predict.cluster)                                    
singler.NovershternHematopoietic.predict.cluster <- SingleR(test = sub.scRNA.harmony@assays$RNA@data, 
                                   ref = NovershternHematopoietic.se, 
                                   labels = NovershternHematopoietic.se$label.fine,
                                   clusters = as.character(sub.scRNA.harmony@meta.data$seurat_clusters))
singler.NovershternHematopoietic.predict.cluster$labels <- rownames(singler.NovershternHematopoietic.predict.cluster) 
singler.MonacoImmune.predict.cluster <- SingleR(test = sub.scRNA.harmony@assays$RNA@data, 
                                   ref = MonacoImmune.se, 
                                   labels = MonacoImmune.se$label.fine,
                                   clusters = as.character(sub.scRNA.harmony@meta.data$seurat_clusters))
singler.MonacoImmune.predict.cluster$labels <- rownames(singler.MonacoImmune.predict.cluster) 
pdf("5.Immune/CD8T/singleR.predicted.score.pdf")
plotScoreHeatmap(singler.ImmGen.predict.cluster)
plotScoreHeatmap(singler.DatabaseImmuneCell.predict.cluster)
plotScoreHeatmap(singler.NovershternHematopoietic.predict.cluster)
plotScoreHeatmap(singler.MonacoImmune.predict.cluster)
dev.off()

####-------------------------------------------------3.marker expression anlalysis
features <- c("LEF1", "IL7R", "CCR7", "SELL", "TCF7", "TGFB1", "CD44", "CD69", "ZNF683", "ITGAE", "ITGA1",  
              "TNF", "IFNG", "KLRG1", "GZMA", "GZMH",
              "NR4A1", "JUNB", "FOS", "ATF3", "DNAJB1", "HSPA1A", "EOMES", 
              "GZMK", "GZMB", "PRF1", "TNFRSF9", "TOX", "ENTPD1", "PDCD1", "CTLA4", "TIGIT", "LAG3", "HAVCR2")
pdf("5.Immune/CD8T/cellType.markerExpression.dot.pdf", height = 2)
DotPlot(sub.scRNA.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + xlab("") + ylab("") + theme(axis.text.x = element_text(size = 8, angle = 90)) + theme(axis.text.y = element_text(size = 8))
DotPlot(sub.scRNA.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + xlab("") + ylab("") + NoLegend() + theme(axis.text.x = element_text(size = 8, angle = 90)) + theme(axis.text.y = element_text(size = 8))
dev.off()
pdf("5.Immune/CD8T/cellType.markerExpression.dot.test.pdf", height = 4)
DotPlot(sub.scRNA.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + xlab("") + ylab("") + theme(axis.text.x = element_text(size = 8, angle = 90)) + theme(axis.text.y = element_text(size = 8))
DotPlot(sub.scRNA.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + xlab("") + ylab("") + NoLegend() + theme(axis.text.x = element_text(size = 8, angle = 90)) + theme(axis.text.y = element_text(size = 8))
dev.off()
avg.expression <- AverageExpression(sub.scRNA.harmony, features, assays = "RNA", slot = "data")
avg.expression.scale <- scale(t(avg.expression$RNA))
avg.expression.max <- decostand(avg.expression$RNA, "range", 1)
pdf("5.Immune/CD8T/cellType.markerExpression.heatmap.pdf")
Heatmap(t(avg.expression.scale), cluster_rows = F, show_column_dend = F, name = "RNA expression",
        width = unit(2, "cm"), height = unit(8, "cm"),
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
Heatmap(avg.expression.max, cluster_rows = F, show_column_dend = F, name = "RNA expression",
        width = unit(2, "cm"), height = unit(8, "cm"),
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))        
dev.off()

pdf("5.Immune/CD8T/markerExpression.pro.pdf")
VlnPlot(sub.scRNA.harmony, features = c("LEF1", "IL7R", "TCF7", "CD44"), ncol = 2)
FeaturePlot(sub.scRNA.harmony, features = c("LEF1", "IL7R", "TCF7", "CD44"), cols = c("lightgrey", "red"), ncol = 2)
VlnPlot(sub.scRNA.harmony, features = c("CD69", "ZNF683", "ITGA1", "ITGAE"), ncol = 2)
FeaturePlot(sub.scRNA.harmony, features = c("CD69", "ZNF683", "ITGA1", "ITGAE"), cols = c("lightgrey", "red"), ncol = 2)
VlnPlot(sub.scRNA.harmony, features = c("JUNB", "FOS", "DNAJB1", "HSPA1A"), ncol = 2)
FeaturePlot(sub.scRNA.harmony, features = c("JUNB", "FOS", "DNAJB1", "HSPA1A"), cols = c("lightgrey", "red"), ncol = 2)
VlnPlot(sub.scRNA.harmony, features = c("PDCD1", "TOX", "ENTPD1", "HAVCR2"), ncol = 2)
FeaturePlot(sub.scRNA.harmony, features = c("PDCD1", "TOX", "ENTPD1", "HAVCR2"), cols = c("lightgrey", "red"), ncol = 2)
dev.off()

# vlnplot
pdf("5.Immune/CD8T/markerExpression.vlnplot.pdf", width = 5)
p1 <- VlnPlot(sub.scRNA.harmony,features = features[1:17], group.by="cellType", same.y.lims=T,flip = T, stack = T) & xlab("") & ylab("Log-normalized expression") & theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.5))
p2 <- VlnPlot(sub.scRNA.harmony,features = features[18:34], group.by="cellType", same.y.lims=T,flip = T, stack = T) & xlab("") & ylab("Log-normalized expression") & theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.5))
ggarrange(p1,p2,ncol=2)
dev.off()
features <- c("CCR7","TCF7","PDCD1","HAVCR2","TOX","CXCR6","XCL1","GZMB","GZMK","IFNG")

####-------------------------------------------------4.Difference analysis
DefaultAssay(sub.scRNA.harmony) <- "RNA"
cluster.DE <- FindAllMarkers(sub.scRNA.harmony, 
                            group.by = "cellType", logfc.threshold = 0, min.pct = 0.1, 
                            test.use = "MAST", latent.vars = "orig.ident")
idents <- levels(sub.scRNA.harmony)
saveFormat <- lapply(idents, function(x){
  index <- which(cluster.DE$cluster == x)
  DEGs <- cluster.DE[index,]
  DEGs.up <- DEGs %>% filter(avg_log2FC>0) %>% arrange(desc(avg_log2FC))
  DEGs.down <- DEGs %>% filter(avg_log2FC<0) %>% arrange(avg_log2FC)
  DEGs <- rbind(DEGs.up, DEGs.down)
  return(DEGs)
})
write.xlsx(saveFormat, file = "5.Immune/CD8T/cellType.DE.xlsx", sheetName = idents, rowNames = F)
saveRDS(cluster.DE, file = "5.Immune/CD8T/cellType.DE.rds")
top.genes <- cluster.DE %>% filter(p_val_adj<0.05 & avg_log2FC>0.25) %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
pdf("5.Immune/CD8T/cellType.topgenes.pdf")
DoHeatmap(sub.scRNA.harmony, features = unique(top.genes$gene), size = 2) + NoLegend()
dev.off()

##-- Functional enrichment analysis
DEGs <- lapply(saveFormat, function(x){
    x <- x %>% filter(p_val_adj<0.05 & avg_log2FC>0.5) %>% arrange(desc(avg_log2FC))
    return(x)
})
names(DEGs) <- idents
write.xlsx(DEGs, file = "5.Immune/CD8T/cellType.DEGs.xlsx", sheetName = idents, rowNames = F)
saveRDS(DEGs, file = "5.Immune/CD8T/cellType.DEGs.rds")

DEGs <- lapply(DEGs, function(x){
    return(x$gene)
})
source("/home/longzhilin/Analysis_Code/PathwayEnrichment/clusterProfiler.enricher.R")
pdf("5.Immune/CD8T/program.pathway.enrichment.pdf")
res <- lapply(names(DEGs), function(x){
    y <- DEGs[[x]]
    res <- cluterProfiler.enricher(gene = y, geneType = "SYMBOL", db.type = "MsigDB", saveDir = paste0(getwd(),"/5.Immune/CD8T/cluterProfiler_MsigDB"),
                            title = x, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
    # gene ratio
    res <- res$em.res.genesymbol@result %>% filter(p.adjust<0.05) #fdr adjust
    pathway.geneNum <- unlist(strsplit(res$BgRatio, "/"))[seq(1, nrow(res),by=2)]
    gene.ratio <- as.numeric(res$Count)/as.numeric(pathway.geneNum)
    res$gene.ratio <- gene.ratio
    return(res)
})
dev.off()
names(res) <- names(DEGs)
write.xlsx(res, file = "5.Immune/CD8T/cluterProfiler_MsigDB/clusterProfiler.enricher.result.xlsx", sheetName = idents, rowNames = F)
#### show cluster profiler result
wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}
library(ggthemes)
enrich.pathways <- c(res[[1]]$ID[c(1, 12, 20, 32, 47, 48, 51)],
                    res[[2]]$ID[c(1, 2, 7, 11, 12, 35, 42, 101)],
                    res[[3]]$ID[c(1, 3, 8, 16, 17, 25, 29, 32, 41, 66)],
                    res[[4]]$ID[c(12, 16, 23, 38, 48, 53, 60, 97, 186)])
enrich.pathways <- unique(enrich.pathways)       
enrich <- lapply(res, function(x){
    idx <- which(x$ID %in% enrich.pathways)
    return(x[idx, c("ID", "p.adjust", "Count")])
})                    
names(enrich) <- names(DEGs)
enrich <- lapply(names(enrich), function(x){
    enrich[[x]]$Type <- rep(x, nrow(enrich[[x]]))
    return(enrich[[x]])
})  
enrich.res <- Reduce(function(x,y) rbind(x,y), enrich)
rownames(enrich.res) <- NULL

enrich.res$p.adjust <- -log10(enrich.res$p.adjust)
enrich.res$wrap <- wrapText(enrich.res$ID, 45)
pdf("5.Immune/CD8T/cluterProfiler_MsigDB/cluterProfiler_MsigDB.enrichment.pdf")
p1 <- ggplot(enrich.res, 
            aes(x = Type, 
                y = wrap, 
                size = Count, 
                color = p.adjust,
                fill = p.adjust))
p2 <- p1 + guides(color=FALSE) + geom_point(shape = 21, alpha = 0.7) + theme_few() + scale_color_gradient(low = "#FF9375", high = "red") + scale_fill_gradient(low = "#FF9375", high = "red", breaks = c(1.3, 2, 3, 4, 5, 6), labels = c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))
p2 <- p2 + theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1, size = 8), axis.text.y = element_text(size = 8)) + xlab("") + ylab("")
print(p2)

#heatmap
pathway.p <- enrich.res[, c("ID", "Type", "p.adjust")]
pathway.data <- pivot_wider(pathway.p, names_from = "Type", values_from = "p.adjust")
pathway.data <- as.data.frame(pathway.data)
rownames(pathway.data) <- pathway.data$ID
pathway.data <- as.data.frame(pathway.data[,-1])
cols <- colorRamp2(c(0, 7.5, 15), c("white", "#ffff66", "red"))
Heatmap(pathway.data, name = "-log10(FDR)", show_column_dend = F, show_row_dend = F, col = cols,  
        na_col = "#f2f3f4", border = "grey", border_gp = gpar(col = "grey"), rect_gp = gpar(col = "grey"),
        width = unit(3, "cm"), height = unit(10, "cm"), cluster_rows = F, cluster_columns = F,
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
dev.off()

##### TCGA survival
source(file = "/home/longzhilin/Analysis_Code/code/analysis.diff.survival.TCGA.R")
DESeq2.normalized_counts <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/DESeq2.normalized_counts.rds")
DESeq2.normalized_counts <- log2(DESeq2.normalized_counts+1)
DESeq2.result <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/DESeq2.result.rds")
clin.data <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/clin.data.rds")
pdf("5.Immune/CD8T/DEGs.survival.pdf")
TCGA.Tissue.resident.C1 <- analysis.diff.survival.TCGA(interest.gene = DEGs[["Tissue-resident.C1"]], diff.gene.pro = DESeq2.result, exp.data.process = DESeq2.normalized_counts, clin.data = clin.data, EnhancedVolcano.plot = F, main = "Tissue-resident.C1", Box.plot = F, meta.signature = T, single.signature = F)
TCGA.Tissue.resident.C2 <- analysis.diff.survival.TCGA(interest.gene = DEGs[["Tissue-resident.C2"]], diff.gene.pro = DESeq2.result, exp.data.process = DESeq2.normalized_counts, clin.data = clin.data, EnhancedVolcano.plot = F, main = "Tissue-resident.C2", Box.plot = F, meta.signature = T, single.signature = F)
TCGA.Exhausted.IEG <- analysis.diff.survival.TCGA(interest.gene = DEGs[["Exhausted IEG"]], diff.gene.pro = DESeq2.result, exp.data.process = DESeq2.normalized_counts, clin.data = clin.data, EnhancedVolcano.plot = F, main = "Exhausted IEG", Box.plot = F, meta.signature = T, single.signature = F)
TCGA.Exhausted <- analysis.diff.survival.TCGA(interest.gene = DEGs[["Exhausted"]], diff.gene.pro = DESeq2.result, exp.data.process = DESeq2.normalized_counts, clin.data = clin.data, EnhancedVolcano.plot = F, main = "Exhausted", Box.plot = F, meta.signature = T, single.signature = F)
dev.off()

##### ICB survial
library(ggpubr)
source(file = "/home/longzhilin/Analysis_Code/SurvivalAnalysis/Cox.function.R")
source(file = "/home/longzhilin/Analysis_Code/code/RCC.ICB.analysis.R")
normalized_expression <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/ICB.therapy/normalized_expression.rds")
patient.info.RNA <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/ICB.therapy/patient.info.RNA.rds")

patient.info.RNA$Sex <- gsub("^F|FEMALE|Female$", 0, patient.info.RNA$Sex)
patient.info.RNA$Sex <- gsub("^M|Male|MALE$", 1, patient.info.RNA$Sex)
patient.info.RNA$Tumor_Sample_Primary_or_Metastasis <- gsub("PRIMARY", 0, patient.info.RNA$Tumor_Sample_Primary_or_Metastasis)
patient.info.RNA$Tumor_Sample_Primary_or_Metastasis <- gsub("METASTASIS", 1, patient.info.RNA$Tumor_Sample_Primary_or_Metastasis)

pdf("5.Immune/CD8T/DEGs.ICB.survival.pdf")
ICB.res <- RCC.icb.analysis(signature.list = DEGs, expresionMatrix = normalized_expression, clincal.info = patient.info.RNA)
dev.off()

####-------------------------------------------------5.Define cell state based on signature
## 1. load signature
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
## 2.1 VISION method
source(file = "/home/longzhilin/Analysis_Code/SingleCell/vision_seurat.R")
source(file = "/home/longzhilin/Analysis_Code/SingleCell/vision.plot.R")
vision_state <- vision_seurat(seurat.object = sub.scRNA.harmony, customize.signature = signature.genes, mc.cores = 36, min_signature_genes = 4)
saveRDS(vision_state, file = "5.Immune/CD8T/vision_signatureScore.rds")

pdf("5.Immune/CD8T/vision_signatureScore.pdf")
res <- vision.plot(vision.obj = vision_state, groupName = "cellType")
res <- vision.plot(vision.obj = vision_state, groupName = "cellType", signature.names = c("Progenitor Exhaustion", "Terminal Exhaustion"), title = "Progenitor & Terminally")
# heatmap
vision.score <- as.data.frame(vision_state@SigScores)
vision.score.mean <- apply(vision.score, 2, function(x){
    score <- tapply(x, sub.scRNA.harmony$cellType, mean)
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
vision.score$group <- sub.scRNA.harmony$cellType
plot.list <- lapply(colnames(vision.score), function(x){
    my_comparisons <- as.list(as.data.frame(combn(levels(sub.scRNA.harmony$cellType),2)))
    a <- vision.score[, c(x, "group")]
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

vision.score.scale <- as.data.frame(scale(vision.score[, c("Progenitor Exhaustion", "Terminal Exhaustion")]))
vision.score.scale$group <- sub.scRNA.harmony$cellType
group.signature.data <- pivot_longer(vision.score.scale, cols = 1:2, names_to = "Type")
p <- ggviolin(group.signature.data, x = "group", y = "value",
          color = "Type", palette = c("#00A087FF", "#F39B7FFF"), fill = "white",
          add = "jitter", add.params = list(size = 0.1)) + ylab("Signature score") + xlab("")
p <- p + rotate_x_text(angle = 45, vjust = 1)+ stat_compare_means(aes(group = Type), label = "p.signif")
print(p)
dev.off()


#### terminal VS progenitor
source(file = "/home/longzhilin/Analysis_Code/Visualization/Plot.EnhancedVolcano.R")
DefaultAssay(sub.scRNA.harmony) <- "RNA"
Progenitor.Terminal.Exhausted <- FindMarkers(sub.scRNA.harmony,
                          ident.1 = "C3",
                          ident.2 = "C2",
                          logfc.threshold = 0, min.pct = 0.1, 
                          test.use = "MAST", latent.vars = "orig.ident")
up.x <- Progenitor.Terminal.Exhausted %>% filter(avg_log2FC>0) %>% arrange(desc(avg_log2FC))
down.x <- Progenitor.Terminal.Exhausted %>% filter(avg_log2FC<=0) %>% arrange(avg_log2FC)
Progenitor.Terminal.Exhausted <- rbind(up.x, down.x)
Progenitor.Terminal.Exhausted$gene <- rownames(Progenitor.Terminal.Exhausted)
write.xlsx(Progenitor.Terminal.Exhausted, file = "5.Immune/CD8T/Progenitor.Terminal.Exhausted.xlsx", rowNames = F)
saveRDS(Progenitor.Terminal.Exhausted, file = "5.Immune/CD8T/Progenitor.Terminal.Exhausted.rds")

pdf("5.Immune/CD8T/Progenitor.Terminal.Exhausted.pdf")
p <- Plot.EnhancedVolcano(Progenitor.Terminal.Exhausted, x = "avg_log2FC", y = "p_val_adj", id.column = "gene", 
                            FC.cutoff = 0.25, selected.showType = c("FC"), select.num = 10, 
                            drawConnectors = T, title = "Terminal VS Progenitor")
print(p)                            
# top30 or down 30 --- log2FC
a <- arrange(Progenitor.Terminal.Exhausted, desc(avg_log2FC))
Progenitor.Terminal.Exhausted.top30 <- a[c(1:30, (nrow(a)-29):nrow(a)), ]
rownames(Progenitor.Terminal.Exhausted.top30) <- Progenitor.Terminal.Exhausted.top30$gene
Progenitor.Terminal.Exhausted.top30 <- Progenitor.Terminal.Exhausted.top30[,"avg_log2FC",drop = F]
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
Heatmap(Progenitor.Terminal.Exhausted.top30, name = "avg_log2FC (scRNA)", row_names_side = "left", col = col_fun, row_names_gp = gpar(fontsize = 9), cluster_rows = F, cluster_columns = F, width = unit(1, "cm"))
dev.off()