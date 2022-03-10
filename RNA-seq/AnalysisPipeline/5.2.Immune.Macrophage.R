#' @description: analysis Macrophage

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
source(file = "/home/longzhilin/Analysis_Code/Plot_colorPaletters.R")
source(file = "/home/longzhilin/Analysis_Code/SingleCell/scRNA.Integrate.multipleSample.R")
source(file = "/home/longzhilin/Analysis_Code/SingleCell/variableFeatureSelection.R")
source(file = "/home/longzhilin/Analysis_Code/code/ratio.plot.R")

data.merge <- readRDS("data.merge.pro.rds")
DefaultAssay(data.merge) <- "RNA"
cell.Type <- c("Macrophage")
sub.scRNA <- subset(data.merge, subset = cellType == cell.Type)
sub.scRNA <- DietSeurat(sub.scRNA, assays = "RNA")
index <- match(paste0("SCT_snn_res.", seq(0.2, 1.2, by=0.1)), colnames(sub.scRNA@meta.data))
sub.scRNA@meta.data <- sub.scRNA@meta.data[,-index]

####----------------------------------------------------- 1. clustring cell
# split the dataset into a list of two seurat objects (stim and CTRL)
sub.list <- SplitObject(sub.scRNA, split.by = "orig.ident")
sub.list <- lapply(X = sub.list, FUN = function(x) {
    x <- SCTransform(x, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
})
sub.scRNA <- merge(sub.list[[1]], y = sub.list[2:length(sub.list)], project = "Renal_Macrophage")
DefaultAssay(sub.scRNA) <- "SCT"
seurat.features.SCT <- SelectIntegrationFeatures(object.list = sub.list, nfeatures = 3000)
VariableFeatures(sub.scRNA) <- seurat.features.SCT
pdf("5.Immune/Macrophage/SCT.Harmony.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "SCT", set.resolutions = seq(0.2, 1.2, by = 0.1), PC = 30, npcs = 50)
dev.off()
sub.scRNA.harmony$seurat_clusters <- sub.scRNA.harmony$SCT_snn_res.0.3
Idents(sub.scRNA.harmony) <- sub.scRNA.harmony$seurat_clusters
DefaultAssay(sub.scRNA.harmony) <- "RNA"
sub.scRNA.harmony <- NormalizeData(sub.scRNA.harmony, verbose = FALSE)
sub.scRNA.harmony <- ScaleData(sub.scRNA.harmony, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mt"), features = rownames(sub.scRNA.harmony))
saveRDS(sub.scRNA.harmony, file = "5.Immune/Macrophage/sub.scRNA.harmony.rds")

####----------------------------------------------------- 2. defined the cell population
DefaultAssay(sub.scRNA.harmony) <- "RNA"
#### 1.marker 
pdf("5.Immune/Macrophage/markerExpression.pdf")              
VlnPlot(sub.scRNA.harmony, features = c("CD3D", "PLVAP", "CD68", "CD163"), ncol = 2)
FeaturePlot(sub.scRNA.harmony, features = c("CD3D", "PLVAP", "CD68", "CD163"), cols = c("lightgrey", "red"), ncol = 2)
VlnPlot(sub.scRNA.harmony, features = c("APOE", "C1QA", "C1QB", "C1QC"), ncol = 2)
FeaturePlot(sub.scRNA.harmony, features = c("APOE", "C1QA", "C1QB", "C1QC"), cols = c("lightgrey", "red"), ncol = 2)
VlnPlot(sub.scRNA.harmony, features = c("SPP1", "CSTB", "FABP5", "FN1"), ncol = 2)
FeaturePlot(sub.scRNA.harmony, features = c("SPP1", "CSTB", "FABP5", "FN1"), cols = c("lightgrey", "red"), ncol = 2)
VlnPlot(sub.scRNA.harmony, features = c("S100A8", "S100A12", "FCGR3A", "FCN1"), ncol = 2)
FeaturePlot(sub.scRNA.harmony, features = c("S100A8", "S100A12", "FCGR3A", "FCN1"), cols = c("lightgrey", "red"), ncol = 2)
dev.off()

features <- c("CD3D", "CD3E", "PLVAP", "ESM1", "ITGAM", "CD68", "CD163", "CSF1R")
pdf("5.Immune/Macrophage/cluster.markerExpression.dot.pdf", height =  2, width = 5)
DotPlot(sub.scRNA.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "seurat_clusters", dot.scale = 3.5) + RotatedAxis() + theme(axis.text.x = element_text(size = 7)) + theme(axis.text.y = element_text(size = 8))
DotPlot(sub.scRNA.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "seurat_clusters", dot.scale = 3.5) + RotatedAxis() + NoLegend() + theme(axis.text.x = element_text(size = 7)) + theme(axis.text.y = element_text(size = 8))
dev.off()
pdf("5.Immune/Macrophage/cluster.markerExpression.dot.test.pdf", height =  4, width = 5)
DotPlot(sub.scRNA.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "seurat_clusters", dot.scale = 3.5) + RotatedAxis() + theme(axis.text.x = element_text(size = 7)) + theme(axis.text.y = element_text(size = 8))
DotPlot(sub.scRNA.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "seurat_clusters", dot.scale = 3.5) + RotatedAxis() + NoLegend() + theme(axis.text.x = element_text(size = 7)) + theme(axis.text.y = element_text(size = 8))
dev.off()

pdf("5.Immune/Macrophage/cluster.ratio.pdf", height = 3, width = 5)
ratio.plot(seurat.object = sub.scRNA.harmony, id.vars1 = "orig.ident", id.vars2 = "seurat_clusters", angle = 60)
ratio.plot(seurat.object = sub.scRNA.harmony, id.vars1 = "orig.ident", id.vars2 = "seurat_clusters", angle = 60)
dev.off()

#### 2. differential expression
DefaultAssay(sub.scRNA.harmony) <- "RNA"
cluster.DE <- FindAllMarkers(sub.scRNA.harmony, only.pos = T,
                            group.by = "seurat_clusters", logfc.threshold = 0.25, min.pct = 0.1, 
                            test.use = "MAST", latent.vars = "orig.ident")
idents <- as.character(unique(cluster.DE$cluster))
saveFormat <- lapply(idents, function(x){
  index <- which(cluster.DE$cluster == x)
  DEGs <- cluster.DE[index,]
  DEGs.up <- DEGs %>% filter(avg_log2FC>0) %>% arrange(desc(avg_log2FC))
  DEGs.down <- DEGs %>% filter(avg_log2FC<0) %>% arrange(avg_log2FC)
  DEGs <- rbind(DEGs.up, DEGs.down)
  return(DEGs)
})
write.xlsx(saveFormat, file = "5.Immune/Macrophage/cluster.DE.xlsx", sheetName = idents, rowNames = F)
saveRDS(cluster.DE, file = "5.Immune/Macrophage/cluster.DE.rds")

####----------------------------------------------------- 3. Remove miscellaneous
sub.scRNA.harmony <- subset(sub.scRNA.harmony, subset = seurat_clusters %in% c(0,1,2))
DefaultAssay(sub.scRNA.harmony) <- "SCT"
pdf("5.Immune/Macrophage/SCT.Harmony.Integration.pro.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA.harmony, assay = "SCT", set.resolutions = seq(0.2, 1.2, by = 0.1), PC = 30, npcs = 50)
dev.off()
sub.scRNA.harmony$seurat_clusters <- sub.scRNA.harmony$SCT_snn_res.0.2
Idents(sub.scRNA.harmony) <- sub.scRNA.harmony$seurat_clusters
sub.scRNA.harmony$cellType2 <- paste0("C", as.numeric(sub.scRNA.harmony$seurat_clusters))
Idents(sub.scRNA.harmony) <- factor(sub.scRNA.harmony$cellType2, levels = c("C1", "C2", "C3"))

cellType3 <- sub.scRNA.harmony$cellType2
cellType3 <- gsub("^C1$", "TAM-C1QB", cellType3)
cellType3 <- gsub("^C2$", "TAM-RGCC", cellType3)
cellType3 <- gsub("^C3$", "TAM-LGALS3", cellType3)
sub.scRNA.harmony$cellType3 <- factor(cellType3, levels = c("TAM-C1QB", "TAM-RGCC", "TAM-LGALS3"))
Idents(sub.scRNA.harmony) <- sub.scRNA.harmony$cellType3
saveRDS(sub.scRNA.harmony, file = "5.Immune/Macrophage/sub.scRNA.harmony.pro.rds")

sub.scRNA.harmony$cellType <- sub.scRNA.harmony$cellType3
Idents(sub.scRNA.harmony) <- sub.scRNA.harmony$cellType
cellType.colors <- c("#F8766D", "#00BA38", "#619CFF")

pdf("5.Immune/Macrophage/cluster.pro.pdf")
DimPlot(object = sub.scRNA.harmony, reduction = 'tsne', pt.size = 1.5, label = TRUE, group.by = "seurat_clusters")+NoLegend()
DimPlot(object = sub.scRNA.harmony, reduction = 'umap', pt.size = 1.5, label = TRUE, group.by = "seurat_clusters")+NoLegend()
dev.off()
pdf("5.Immune/Macrophage/cellType.pdf")
DimPlot(object = sub.scRNA.harmony, reduction = 'tsne', pt.size = 1.5, label = TRUE, group.by = "cellType")+NoLegend()
DimPlot(object = sub.scRNA.harmony, reduction = 'umap', pt.size = 1.5, label = TRUE, group.by = "cellType")+NoLegend()
dev.off()

features <- c("IER2", "JUN", "F13A1", "APOE", "C1QB", "C1QC", "C1QA",
              "EREG", "NLRP3", "AREG", "SLC2A3", "RGCC", "CLEC5A",
              "LILRB4", "MARCO", "ANXA2", "LGALS3", "GPNMB", "TREM2")
pdf("5.Immune/Macrophage/cluster.markerExpression.dot.pro.pdf", height = 1.5, width = 6)
DotPlot(sub.scRNA.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + theme(axis.text.x = element_text(size = 7, angle = 90)) + theme(axis.text.y = element_text(size = 8))
DotPlot(sub.scRNA.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + NoLegend() + theme(axis.text.x = element_text(size = 7, angle = 90)) + theme(axis.text.y = element_text(size = 8))
dev.off()
pdf("5.Immune/Macrophage/cluster.markerExpression.dot.pro.test.pdf", height = 4, width = 6)
DotPlot(sub.scRNA.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + theme(axis.text.x = element_text(size = 7, angle = 90)) + theme(axis.text.y = element_text(size = 8))
DotPlot(sub.scRNA.harmony, features = features, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + NoLegend() + theme(axis.text.x = element_text(size = 7, angle = 90)) + theme(axis.text.y = element_text(size = 8))
dev.off()

pdf("5.Immune/Macrophage/markerExpression.pro.pdf")              
VlnPlot(sub.scRNA.harmony, features = c("IL1A", "IL1B", "IFITM1", "IFITM2"), ncol = 2)
FeaturePlot(sub.scRNA.harmony, features = c("IL1A", "IL1B", "IFITM1", "IFITM2"), cols = c("lightgrey", "red"), ncol = 2)
VlnPlot(sub.scRNA.harmony, features = c("FLT1", "SPARC", "RGS5", "MRC1"), ncol = 2)
FeaturePlot(sub.scRNA.harmony, features = c("FLT1", "SPARC", "RGS5", "MRC1"), cols = c("lightgrey", "red"), ncol = 2)
dev.off()

####----------------------------------------------------- 2.differential expression
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
names(saveFormat) <- idents
write.xlsx(saveFormat, file = "5.Immune/Macrophage/cellType.DE.pro.xlsx", sheetName = idents, rowNames = F)
saveRDS(cluster.DE, file = "5.Immune/Macrophage/cellType.DE.pro.rds")

####-- pathway enrichment
DEGs <- lapply(saveFormat, function(x){
    x <- x %>% filter(p_val_adj<0.05 & avg_log2FC>0.5) %>% arrange(desc(avg_log2FC))
    return(x)
})
names(DEGs) <- idents
write.xlsx(DEGs, file = "5.Immune/Macrophage/cellType.DEGs.xlsx", sheetName = idents, rowNames = F)
saveRDS(DEGs, file = "5.Immune/Macrophage/cellType.DEGs.rds")

DEGs <- lapply(DEGs, function(x){
    return(x$gene)
})
source("/home/longzhilin/Analysis_Code/PathwayEnrichment/clusterProfiler.enricher.R")
pdf("5.Immune/Macrophage/program.pathway.enrichment.pdf")
res <- lapply(names(DEGs), function(x){
    y <- DEGs[[x]]
    res <- cluterProfiler.enricher(gene = y, geneType = "SYMBOL", db.type = "MsigDB", saveDir = paste0(getwd(),"/5.Immune/Macrophage/cluterProfiler_MsigDB"),
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
library(openxlsx)
write.xlsx(res, file = "5.Immune/Macrophage/cluterProfiler_MsigDB/clusterProfiler.enricher.result.xlsx", sheetName = idents, rowNames = F)
#### show the cluster profiler result
wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}
library(ggthemes)
enrich.pathways <- c(res[[1]]$ID[c(1, 3, 17, 26, 27, 30, 32, 46, 68, 73, 75, 99, 119, 196, 216, 397, 430)],
                    res[[2]]$ID[c(1:5, 7, 8, 10, 18, 19, 61, 72, 122, 127)],
                    res[[3]]$ID[c(1, 2, 3, 5, 6, 8, 9, 38, 74, 173)])
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
pdf("5.Immune/Macrophage/cluterProfiler_MsigDB/cluterProfiler_MsigDB.enrichment.pdf")
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
cols <- colorRamp2(c(0, 10, 20), c("white", "#ffff66", "red"))
Heatmap(pathway.data, name = "-log10(FDR)", show_column_dend = F, show_row_dend = F, col = cols,  
        na_col = "#f2f3f4", border = "grey", border_gp = gpar(col = "grey"), rect_gp = gpar(col = "grey"),
        width = unit(3, "cm"), height = unit(12, "cm"), cluster_rows = F, cluster_columns = F,
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
dev.off()

##### TCGA survival
source(file = "/home/longzhilin/Analysis_Code/code/analysis.diff.survival.TCGA.R")
DESeq2.normalized_counts <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/DESeq2.normalized_counts.rds")
DESeq2.normalized_counts <- log2(DESeq2.normalized_counts+1)
DESeq2.result <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/DESeq2.result.rds")
clin.data <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/clin.data.rds")
pdf("5.Immune/Macrophage/DEGs.survival.pdf")
TCGA.TAM.C1QB <- analysis.diff.survival.TCGA(interest.gene = DEGs[["TAM-C1QB"]], diff.gene.pro = DESeq2.result, exp.data.process = DESeq2.normalized_counts, clin.data = clin.data, EnhancedVolcano.plot = F, main = "TAM-C1QB", Box.plot = F, meta.signature = T, single.signature = F)
TCGA.TAM.RGCC <- analysis.diff.survival.TCGA(interest.gene = DEGs[["TAM-RGCC"]], diff.gene.pro = DESeq2.result, exp.data.process = DESeq2.normalized_counts, clin.data = clin.data, EnhancedVolcano.plot = F, main = "TAM-RGCC", Box.plot = F, meta.signature = T, single.signature = F)
TCGA.TAM.LGALS3 <- analysis.diff.survival.TCGA(interest.gene = DEGs[["TAM-LGALS3"]], diff.gene.pro = DESeq2.result, exp.data.process = DESeq2.normalized_counts, clin.data = clin.data, EnhancedVolcano.plot = F, main = "TAM-LGALS3", Box.plot = F, meta.signature = T, single.signature = F)
dev.off()

##### ICB survival
library(ggpubr)
source(file = "/home/longzhilin/Analysis_Code/SurvivalAnalysis/Cox.function.R")
source(file = "/home/longzhilin/Analysis_Code/code/RCC.ICB.analysis.R")
normalized_expression <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/ICB.therapy/normalized_expression.rds")
patient.info.RNA <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/ICB.therapy/patient.info.RNA.rds")
patient.info.RNA$Sex <- gsub("^F|FEMALE|Female$", 0, patient.info.RNA$Sex)
patient.info.RNA$Sex <- gsub("^M|Male|MALE$", 1, patient.info.RNA$Sex)
patient.info.RNA$Tumor_Sample_Primary_or_Metastasis <- gsub("PRIMARY", 0, patient.info.RNA$Tumor_Sample_Primary_or_Metastasis)
patient.info.RNA$Tumor_Sample_Primary_or_Metastasis <- gsub("METASTASIS", 1, patient.info.RNA$Tumor_Sample_Primary_or_Metastasis)

pdf("5.Immune/Macrophage/DEGs.ICB.survival.pdf")
ICB.res <- RCC.icb.analysis(signature.list = DEGs, expresionMatrix = normalized_expression, clincal.info = patient.info.RNA)
dev.off()


####----------------------------------------------------- 3.Define cell state based on signature
## M1 和 M2
source(file = "/home/longzhilin/Analysis_Code/SingleCell/vision_seurat.R")
source(file = "/home/longzhilin/Analysis_Code/SingleCell/vision.plot.R")
signature.genes <- read.table("/data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/M1&M2&angiogenic&phagocytic.Cheng.2021.Cell.txt", header = T, stringsAsFactors = F, sep = "\t")
vision_state <- vision_seurat(seurat.object = sub.scRNA.harmony, customize.signature = signature.genes, mc.cores = 36, min_signature_genes = 4)
saveRDS(vision_state, file = "5.Immune/Macrophage/vision_signatureScore.rds")

pdf("5.Immune/Macrophage/vision_signatureScore.pdf")
res <- vision.plot(vision.obj = vision_state, groupName = "cellType")
res <- vision.plot(vision.obj = vision_state, groupName = "cellType", signature.names = c("M1-state", "M2-state"), title = "M1 & M2")
# heatmap
vision.score <- as.data.frame(vision_state@SigScores)
vision.score.mean <- apply(vision.score, 2, function(x){
    score <- tapply(x, sub.scRNA.harmony$cellType, mean)
    return(score)
})
p <- Heatmap(t(vision.score.mean), width = unit(6, "cm"), height = unit(6, "cm"), name = "Signature score",
        show_row_dend = F, show_column_dend = F, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
print(p)
vision.score.mean.scale <- scale(vision.score.mean) 
p <- Heatmap(t(vision.score.mean.scale), width = unit(6, "cm"), height = unit(6, "cm"), name = "Signature score",
        show_row_dend = F, show_column_dend = F, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
print(p)

# violin
vision.score$group <- sub.scRNA.harmony$cellType
plot.list <- lapply(colnames(vision.score), function(x){
    my_comparisons <- as.list(as.data.frame(combn(levels(sub.scRNA.harmony$cellType),2)))
    a <- vision.score[, c(x, "group")]
    names(a) <- c("Signature score", "group")
    p <- ggboxplot(a, x = "group", y = "Signature score", title = x, color = "black", fill = "group", add = "jitter", add.params = list(color = "black", size = 0.1)) + 
         xlab("") + rotate_x_text(angle = 45, vjust = 1) + scale_fill_manual(values = cellType.colors) + 
         stat_compare_means(comparisons = my_comparisons, label = "p.signif") + NoLegend()
    return(p)
})
p <- ggarrange(plotlist = plot.list[1:4],ncol = 2, nrow = 2)
print(p)
dev.off()

#### M1 and M2 expression matrix
library(reshape2)
signature.genes <- signature.genes[which(signature.genes$Type %in% c("M1-state", "M2-state")),]
idx <- match(signature.genes$Gene, rownames(sub.scRNA.harmony))
signature.genes <- signature.genes[which(!is.na(idx)),]
signature.names <- unique(signature.genes$Type)
signature.list <- lapply(signature.names, function(x){
    genes <- signature.genes$Gene[which(signature.genes$Type == x)]
    return(genes)
})
names(signature.list) <- signature.names
#signature gene matrix
exp.matrix <- AverageExpression(sub.scRNA.harmony, signature.genes$Gene, assay = "RNA", slot = "scale.data")
exp.matrix <- exp.matrix$RNA
col.split <- signature.genes$Type
exp.matrix.scale <- scale(t(exp.matrix))
pdf("5.Immune/Macrophage/M1&M2.signature.score.pdf")
Heatmap(t(exp.matrix), name = "Expression", cluster_columns = T, column_split = col.split, cluster_rows = T, show_row_dend = F, show_column_dend = F, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 7), width = unit(10, "cm"), height = unit(3, "cm"))
Heatmap(exp.matrix.scale, name = "Expression", cluster_columns = T, column_split = col.split, cluster_rows = T, show_row_dend = F, show_column_dend = F, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 7), width = unit(10, "cm"), height = unit(3, "cm"))
dev.off()

####----Observe other signatrue patterns----####
Immune.checkpoint.evasion.genes <- c("ICOSLG", "CD80", "CD86", "VSIR", "VSIG4", "LGALS9", "CD274", "PDCD1LG2", "SIGLEC10")
pdf("5.Immune/Macrophage/Immune.checkpoint.evasion.genes.pdf", height = 1.5, width = 5)
DotPlot(sub.scRNA.harmony, features = Immune.checkpoint.evasion.genes, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + theme(axis.text.x = element_text(size = 7, angle = 90)) + theme(axis.text.y = element_text(size = 8))
DotPlot(sub.scRNA.harmony, features = Immune.checkpoint.evasion.genes, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + NoLegend() + theme(axis.text.x = element_text(size = 7, angle = 90)) + theme(axis.text.y = element_text(size = 8))
dev.off()
pdf("5.Immune/Macrophage/Immune.checkpoint.evasion.genes.test.pdf", height = 4, width = 5)
DotPlot(sub.scRNA.harmony, features = Immune.checkpoint.evasion.genes, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + theme(axis.text.x = element_text(size = 7, angle = 90)) + theme(axis.text.y = element_text(size = 8))
DotPlot(sub.scRNA.harmony, features = Immune.checkpoint.evasion.genes, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + NoLegend() + theme(axis.text.x = element_text(size = 7, angle = 90)) + theme(axis.text.y = element_text(size = 8))
dev.off()

###--- MHC molecular
signature.genes <- read.table("/data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/antigenPresenting.LeiZhang.2021.Cell&Guillermo.2014.frontiersImmunology.txt", header = T, stringsAsFactors = F, sep = "\t")
idx <- match(signature.genes$Gene, rownames(sub.scRNA.harmony))
signature.genes <- signature.genes[which(!is.na(idx)),]
signature.genes <- signature.genes[-c(4, 9, 15, 21, 23, 25, 26, 32:35),]
exp.matrix <- AverageExpression(sub.scRNA.harmony, signature.genes$Gene, assay = "RNA", slot = "data")
exp.matrix <- exp.matrix$RNA
exp.matrix.scale <- scale(t(exp.matrix))
col.split <- factor(signature.genes$Type, levels = c("MCH-I", "MCH-II", "Pro-inflammatory cytokine", "Anti-inflammatory cytokine", "Chemokine", "Other"))
pdf("5.Immune/Macrophage/MHC.expression_cytokines.pdf")
Heatmap(exp.matrix, name = "Expression", cluster_rows = T, row_split = col.split, cluster_columns = T, show_row_dend = F, show_column_dend = F, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8), height = unit(12, "cm"), width = unit(1.5, "cm"))
Heatmap(t(exp.matrix.scale), name = "Expression", cluster_rows = T, row_split = col.split, cluster_columns = T, show_row_dend = F, show_column_dend = F, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8), height = unit(12, "cm"), width = unit(1.5, "cm"))
DotPlot(sub.scRNA.harmony, features = signature.genes$Gene, cols = c("#1e90ff", "#F15F30"), group.by = "cellType", dot.scale = 3.5) + RotatedAxis() + theme(axis.text.y = element_text(size = 8)) + labs(x= "", y = "") + theme(axis.text.x = element_text(size = 8))
dev.off()


## version 
signature.genes <- read.table("/data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/M1&M2&angiogenic&phagocytic.Cheng.2021.Cell.txt", header = T, stringsAsFactors = F, sep = "\t")
signature.list <- lapply(unique(signature.genes$Type), function(x){
    idx <- which(signature.genes$Type == x)
    return(signature.genes$Gene[idx])
})
names(signature.list) <- unique(signature.genes$Type)
sub.scRNA.harmony <- AddModuleScore(sub.scRNA.harmony, features = signature.list)
signature.scores <- sub.scRNA.harmony@meta.data[, grep("Cluster", colnames(sub.scRNA.harmony@meta.data))]
colnames(signature.scores) <- names(signature.list)

# violin
signature.scores$group <- factor(sub.scRNA.harmony$cellType, levels = c("TAM-C1QB", "TAM-RGCC", "TAM-LGALS3"))
pdf("5.Immune/Macrophage/seurat.signatureScore.pdf")
plot.list <- lapply(colnames(signature.scores)[1:4], function(x){
    my_comparisons <- as.list(as.data.frame(combn(levels(signature.scores$group),2)))
    a <- signature.scores[, c(x, "group")]
    names(a) <- c("Signature score", "group")
    p <- ggboxplot(a, x = "group", y = "Signature score", title = x, color = "black", fill = "group", add = "jitter", add.params = list(color = "black", size = 0.1)) + 
         xlab("") + rotate_x_text(angle = 45, vjust = 1) + scale_fill_manual(values = cellType.colors) + 
         stat_compare_means(comparisons = my_comparisons, label = "p.signif") + NoLegend()
    return(p)
})
p <- ggarrange(plotlist = plot.list[1:4],ncol = 2, nrow = 2)
print(p)
dev.off()
