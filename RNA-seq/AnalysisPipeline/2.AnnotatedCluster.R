#' @description: annotated the cell type

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
library(openxlsx)
library(future)
plan("multiprocess", workers = 10) 
options(future.globals.maxSize = 50000 * 1024^2) # set 50G RAM
setwd(dir = "/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA")
source(file = "/home/longzhilin/Analysis_Code/code/ratio.plot.R")

#### Harmony corrected result
data.merge <- readRDS("data.merge.harmony.PC40.SCT.feature3000.rds")
DefaultAssay(data.merge) <- "RNA"

pdf("2.Cluster/cluster.pdf")
DimPlot(object = data.merge, reduction = 'tsne',label = TRUE, group.by = "seurat_clusters")+NoLegend()
DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "seurat_clusters")+NoLegend()
dev.off()

#### Plot the ratio of each cell type and sample situation
pdf("2.Cluster/cluster.number.ratio.pdf", height = 4, width = 7)
ratio.plot(seurat.object = data.merge, id.vars1 = "orig.ident", id.vars2 = "seurat_clusters")
dev.off()

expMatrix <- GetAssayData(data.merge, slot = "scale.data")
highVariableGenes <- VariableFeatures(data.merge)
expMatrix.high <- expMatrix[highVariableGenes,]
meanExpCluster <- apply(expMatrix.high, 1, function(x){
    mean.value <- tapply(x, data.merge$seurat_clusters, mean)
    return(mean.value)
})

corrMatrix <- (1- cor(t(meanExpCluster), method="pearson"))/2
library(ape)
## dd <- dist(M)
hc <- hclust(as.dist(corrMatrix),method="complete")
pdf("2.Cluster/clustersimilarity.HVG.pdf")
plot.phylo(as.phylo(hc), type = "fan", cex = 0.8,no.margin = FALSE)
dev.off()

#### Differential expression
Idents(data.merge) <- data.merge$seurat_clusters
cluster.all.markers <- FindAllMarkers(data.merge, only.pos = TRUE, group.by = "seurat_clusters", test.use = "MAST", latent.vars = "orig.ident")
cluster.sig.markers <- cluster.all.markers[which(cluster.all.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "2.Cluster/AnnotateCellType/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "2.Cluster/AnnotateCellType/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.table(cluster.sig.markers.top, file = "2.Cluster/AnnotateCellType/cluster.sig.markers.top30.txt", col.names = T, row.names = F, sep = "\t", quote = F)
top.genes <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
pdf("2.Cluster/AnnotateCellType/cluster.topgenes.pdf")
DoHeatmap(data.merge, features = unique(top.genes$gene), size = 2) + NoLegend()
dev.off()

####Method 1: The expression of the classic marker
cell.type.markers <- read.table(file = "2.Cluster/AnnotateCellType/cellMarker.txt", header = T, stringsAsFactors = F, sep = "\t")
exp.matrix <- GetAssayData(data.merge, slot = "data")
index <- match(cell.type.markers$Gene, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]
cell.type.markers <- cell.type.markers[which(!is.na(index)),]

# Calculate the average expression of the cell type of each cluster
cluster.score <- apply(gene.matrix, 1, function(x){
    a <- tapply(x, data.merge$seurat_clusters, mean)
    return(a)
})
cluster.score.normailzed <- decostand(cluster.score, "range", 2) 
cellType.cluster.score <- apply(cluster.score, 1, function(x){
    a <- tapply(x, cell.type.markers$CellType, mean)
    return(a)
})
cellType.cluster.score.normailzed <- decostand(cellType.cluster.score, "range", 1)
annotation.colors <- Palettes$stallion2[1:length(unique(cell.type.markers$CellType))]
names(annotation.colors) <- unique(cell.type.markers$CellType)
row.annotations <- rowAnnotation(Type = factor(cell.type.markers$CellType, 
                                 levels = unique(cell.type.markers$CellType)),
                                 col = list(Type = annotation.colors))
pdf("2.Cluster/AnnotateCellType/cluster.signature.expression.pdf")
col_fun1 <- colorRamp2(c(0, 3), c("white", "#ff5a36"))
col_fun2 <- colorRamp2(c(0, 0.5, 1), c("#1e90ff", "white", "#ff5a36"))

row_split <- factor(cell.type.markers$CellType, levels = unique(cell.type.markers$CellType))
Heatmap(t(cluster.score.normailzed), col = col_fun2, row_split = row_split, width = unit(10, "cm"), height = unit(12, "cm"), cluster_columns = F, cluster_rows = F, show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6),
        heatmap_legend_param = list(
                title = "Expression", at = c(0, 1), 
                labels = c("min", "max")))
Heatmap(cellType.cluster.score, name = "Expression", col = col_fun1, width = unit(8, "cm"), height = unit(8, "cm"), cluster_columns = T , cluster_rows = T, show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 6))
Heatmap(cellType.cluster.score.normailzed, col = col_fun2, width = unit(8, "cm"), height = unit(8, "cm"), cluster_columns = T , cluster_rows = F, show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 6), 
        heatmap_legend_param = list(
                title = "Expression", at = c(0, 1), 
                labels = c("min", "max")))

a <- data.merge
a$seurat_clusters <- factor(a$seurat_clusters, levels = rownames(cluster.score.normailzed))
DotPlot(a, cols = c("#1e90ff", "#ff5a36"), features = cell.type.markers$Gene, group.by = "seurat_clusters", dot.scale = 4) + RotatedAxis() + theme(axis.text.x = element_text(size = 6)) + theme(axis.text.y = element_text(size = 6))
DotPlot(a, cols = c("#1e90ff", "#ff5a36"), features = cell.type.markers$Gene, group.by = "seurat_clusters", dot.scale = 4) + RotatedAxis() + NoLegend()+ theme(axis.text.x = element_text(size = 8)) + theme(axis.text.y = element_text(size = 8))
dev.off()

#####annotated the cell type
#RNA low resolution
cluster.label <- data.merge@meta.data$seurat_clusters
cluster.label <- gsub("^0$", "NK/NKT cell", cluster.label)
cluster.label <- gsub("^1$", "Macrophage", cluster.label)
cluster.label <- gsub("^2$", "CD8+ T cell", cluster.label)
cluster.label <- gsub("^3$", "CD4+ T cell", cluster.label)
cluster.label <- gsub("^4$", "Tumor", cluster.label)
cluster.label <- gsub("^5$", "Endothelium (VCAM1-)", cluster.label)
cluster.label <- gsub("^6$", "Mesangial cell", cluster.label)
cluster.label <- gsub("^7$", "Monocyte", cluster.label) # --- monocyte CD14
cluster.label <- gsub("^8$", "Mast cell", cluster.label)
cluster.label <- gsub("^9$", "B cell", cluster.label)
cluster.label <- gsub("^10$", "B cell", cluster.label)
cluster.label <- gsub("^11$", "Monocyte", cluster.label) # --- monocyte CD16
cluster.label <- gsub("^12$", "Endothelium (VCAM1+)", cluster.label)
cluster.label <- gsub("^13$", "Treg", cluster.label)
cluster.label <- gsub("^14$", "Dendritic cell", cluster.label)
cluster.label <- gsub("^15$", "Unknown1", cluster.label)
cluster.label <- gsub("^16$", "Proliferative CD8+ T cell", cluster.label)
cluster.label <- gsub("^17$", "Neutrophil", cluster.label)
cluster.label <- gsub("^18$", "Unknown2", cluster.label)
data.merge <- AddMetaData(data.merge, cluster.label, col.name = "cellType_low")
# Remove cells that cannot be annotated or the cell number less than 100 cells
data.merge.pro <- subset(data.merge, subset = seurat_clusters %in% c(0:14, 16,17))
data.merge.pro$seurat_clusters <- factor(data.merge.pro$seurat_clusters, levels = c(0:14, 16,17))
saveRDS(data.merge.pro, "data.merge.pro.rds")

pdf("2.Cluster/AnnotateCellType/cellType.pro.pdf")
DimPlot(object = data.merge, reduction = 'tsne',label = TRUE, group.by = "cellType_low")+NoLegend()
DimPlot(object = data.merge, reduction = 'tsne',label = TRUE, group.by = "cellType_low")
DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "cellType_low")+NoLegend()
DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "cellType_low")
dev.off()

##Plot--- celltype marker plot
plot.cellType <- c("CD4+ T cell", "Treg", "CD8+ T cell", "NK/NKT cell", "Proliferative CD8+ T cell", "B cell", "Macrophage", "Monocyte", "Dendritic cell", "Neutrophil", "Mast cell", "Endothelium (VCAM1+)", "Endothelium (VCAM1-)", "Mesangial cell", "Tumor")
pdf("2.Cluster/AnnotateCellType/cellType.ratio.pdf", height = 4, width = 7)
ratio.plot(seurat.object = data.merge, id.vars1 = "orig.ident", id.vars2 = "cellType_low", angle = 60)
dev.off()

pdf("2.Cluster/AnnotateCellType/cellType.marker.dotplot.pdf", height = 4)
data.merge$cellType_low <- factor(data.merge$cellType_low, levels = plot.cellType)
DotPlot(data.merge, cols = c("#1e90ff", "#F15F30"), features = cell.type.markers$Gene, group.by = "cellType_low", dot.scale = 3.5)  + NoLegend() + theme(axis.text.y = element_text(size = 4)) + labs(x= "", y = "") + theme(axis.text.x = element_text(size = 6, angle = 90))
DotPlot(data.merge, cols = c("#1e90ff", "#F15F30"), features = cell.type.markers$Gene, group.by = "cellType_low", dot.scale = 3.5)  + theme(axis.text.y = element_text(size = 4)) + labs(x= "", y = "") + theme(axis.text.x = element_text(size = 5, angle = 90))
dev.off()

#### Cell type specific gene
Idents(data.merge) <- data.merge$cellType_low
idents <- as.character(levels(data.merge))
cellType.all.markers <- FindAllMarkers(data.merge, 
                                       group.by = "cellType_low", 
                                       logfc.threshold = 0, 
                                       min.pct = 0.1, 
                                       test.use = "MAST", 
                                       latent.vars = "orig.ident")
saveFormat <- lapply(idents, function(x){
  index <- which(cellType.all.markers$cluster == x)
  DEGs <- cellType.all.markers[index,]
  DEGs.up <- DEGs %>% filter(avg_log2FC>0) %>% arrange(desc(avg_log2FC))
  DEGs.down <- DEGs %>% filter(avg_log2FC<0) %>% arrange(avg_log2FC)
  DEGs <- rbind(DEGs.up, DEGs.down)
  return(DEGs)
})
write.xlsx(saveFormat, file = "2.Cluster/AnnotateCellType/celltype.all.DEGs.xlsx", sheetName = idents, rowNames = F)
saveRDS(cellType.all.markers, file = "2.Cluster/AnnotateCellType/cellType.all.DEGs.rds")

#require logfc.threshold >= 0.25 & p_val_adj < 0.05
cellType.sig.DEGs <- cellType.all.markers %>% filter(avg_log2FC >=0.25 & p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))
saveFormat <- lapply(idents, function(x){
  index <- which(cellType.sig.DEGs$cluster == x)
  DEGs <- cellType.sig.DEGs[index,]
  DEGs <- DEGs %>% arrange(desc(avg_log2FC))
  return(DEGs)
})
write.xlsx(saveFormat, file = "2.Cluster/AnnotateCellType/cellType.sig.pos.DEGs.xlsx", sheetName = idents, rowNames = F)
saveRDS(cellType.sig.DEGs, file = "2.Cluster/AnnotateCellType/cellType.sig.pos.DEGs.rds")
top.genes <- cellType.sig.DEGs %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
pdf("2.Cluster/AnnotateCellType/cellType.topgenes.pdf")
DoHeatmap(data.merge, features = unique(top.genes$gene), size = 2) + NoLegend()
dev.off()

Tumor.info <- data.merge@meta.data[which(data.merge@meta.data$cellType_low=="Tumor"),]
cell.number <- as.data.frame(table(Tumor.info$orig.ident))
pdf("2.Cluster/AnnotateCellType/Tumor.patient.pdf")
ggbarplot(cell.number, x="Var1", y="Freq", fill = "Var1", color = "Var1", palette = Palettes$group_pal[1:length(unique(Tumor.info$orig.ident))],
          sort.by.groups=FALSE,
          label = T, xlab = "", ylab = "Cell Number") + theme(legend.position="none") 
dev.off()

cellType.ratio <- as.data.frame(table(data.merge$cellType_low))
cellType.ratio$Type <- rep("Lymphoid cells", nrow(cellType.ratio))
idx <- which(cellType.ratio$Var1 %in% c("Macrophage", "Monocyte", "Dendritic cell", "Neutrophil", "Mast cell"))
cellType.ratio$Type[idx] <- "Myeloid cells"
idx <- which(cellType.ratio$Var1 %in% c("Endothelium (VCAM1+)", "Endothelium (VCAM1-)", "Mesangial cell"))
cellType.ratio$Type[idx] <- "Other cells"
idx <- which(cellType.ratio$Var1 == "Tumor")
cellType.ratio$Type[idx] <- "Tumor cells"
group.ratio <- tapply(cellType.ratio$Freq, cellType.ratio$Type, sum)
group.ratio <- data.frame(Type = names(group.ratio), Ratio = group.ratio)
labs <- paste0(group.ratio$Type, " (", round(group.ratio$Ratio/sum(group.ratio$Ratio), 4)*100, "%)")
pdf("2.Cluster/AnnotateCellType/lineage.ratio.pdf")
p <- ggpie(group.ratio, "Ratio", label = labs,
      fill = "Type", color = "white", lab.pos = "in",
      palette = "npgs")
print(p)
dev.off()

#### functional enrichment
## GSEA software --- recalculate the result
a <- data.frame(geneName = Tumor.DEGs$gene, Rank = Tumor.DEGs$avg_log2FC)
write.csv(a, file = "2.Cluster/AnnotateCellType/GSEA/Tumor.DEGs.csv")

#### The proportion of various cell types in each patient
library(plotly)
pdf("2.Cluster/AnnotateCellType/cellType.patient.ratio.pdf")
patients <- unique(data.merge$orig.ident)
res <- sapply(patients, function(x){
  a <- subset(data.merge, subset = orig.ident==x)
  cellType.ratio <- as.data.frame(table(as.character(a$cellType_low)))
  group.ratio <- cellType.ratio$Freq/sum(cellType.ratio$Freq)
  group.ratio <- data.frame(Type = cellType.ratio$Var1, Ratio = group.ratio)
  labs <- paste0(round(group.ratio$Ratio, 4)*100, "%")
  p <- ggpie(group.ratio, "Ratio", label = labs,
        fill = "Type", color = "white", lab.pos = "out",
        palette = "npgs")
  print(p)
})
dev.off()
