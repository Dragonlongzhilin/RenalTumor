#' @description peak calling

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86) #---GRCh38 (hg38)
library(ggplot2)
library(patchwork)
set.seed(101)
library(GenomicRanges)
library(ggpubr)
library(tibble)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(openxlsx)

source(file = "/home/longzhilin/Analysis_Code/Visualization/colorPalettes.R")
source(file = "/home/longzhilin/Analysis_Code/Combined.P.FC.R")
setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC")

scATAC.data <- readRDS("scATAC.data.rds")
Idents(scATAC.data) <- scATAC.data$AnnotatedcellType
DefaultAssay(scATAC.data) <- "ATAC"

####To call peaks on each annotated cell type, we can use the group.by argument
peaks <- CallPeaks(
  object = scATAC.data,
  group.by = "AnnotatedcellType",
  macs2.path = "/home/longzhilin/miniconda3/envs/SingleCell/bin/macs2",
  outdir = "/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210721/scATAC/4.Peak"
)
saveRDS(peaks, "4.Peak/cellType.peak.rds") 

library(future)
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 100000 * 1024^2) 
# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(scATAC.data), # from cellranger fragment result
  features = peaks,
  cells = colnames(scATAC.data)
)
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"
saveRDS(annotation, file = "4.Peak/annotation.rds")
# create a new assay using the MACS2 peak set and add it to the Seurat object
scATAC.data[["Peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(scATAC.data),
  annotation = annotation,
  genome = "hg38"
)
DefaultAssay(scATAC.data) <- "Peaks"
gene.activities <- GeneActivity(scATAC.data)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
scATAC.data[['Macs2ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
scATAC.data <- NormalizeData(
  object = scATAC.data,
  assay = 'Macs2ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(scATAC.data$nCount_Peaks)
)
saveRDS(scATAC.data, "scATAC.data.pro.rds") #### macs2 calling

pdf("4.Peak/CA9.macs2.pdf")
CoveragePlot(
  object = scATAC.data,
  region = "CA9",
  assay = "Peaks",
  features = "CA9",
  ranges.title = "MACS2",
  expression.assay = "Macs2ACTIVITY",
  annotation = TRUE,
  peaks = F,
  links = F
)
tile_plot <- TilePlot(
  object = scATAC.data,
  region = "CA9"
)
print(tile_plot)
dev.off()

pdf("4.Peak/CA9.pdf")
CoveragePlot(
  object = scATAC.data,
  assay = "ATAC",
  expression.assay = "ACTIVITY",
  region = "CA9",
  features = "CA9",
  annotation = TRUE,
  peaks = TRUE,
  links = TRUE
)
dev.off()

###coverage plot of marker genes
source("/home/longzhilin/Analysis_Code/SingleCell/FindRegion.R")
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
plot.cellType <- rev(c("CD4+ T cell", "Treg", "CD8+ T cell", "NK/NKT cell", "B cell", "Macrophage", "Monocyte", "Mast cell", "Endothelium (VCAM1+)", "Endothelium (VCAM1-)", "Mesangial cell", "Tumor"))
scATAC.data$AnnotatedcellType <- factor(scATAC.data$AnnotatedcellType, levels = plot.cellType)
DefaultAssay(scATAC.data) <- "Peaks"
Idents(scATAC.data) <- scATAC.data$AnnotatedcellType

cell.type.markers <- read.table(file = "/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/2.Cluster/AnnotateCellType/cellMarker.txt", header = T, stringsAsFactors = F, sep = "\t")
genes <- cell.type.markers$Gene
genes <- genes[-19]

genes <- c("CD8A", "CD4", "FGFBP2", "MS4A1", "CSF1R", "CD163", "S100A12", "TPSAB1", "KDR", "VCAM1", "PDGFRB", "CA9")
pdf("2.Cluster/AnnotateCellType/cellType.coverage.plot.origin2.pdf", height = unit(3, "inches"))
res <- sapply(genes, function(x){
  cat(x, "...\n")
  regions <- FindRegion(object = scATAC.data, region = x, assay = "Peaks", extend.upstream = 1000, extend.downstream = 1000)
  idx <- data.frame(findOverlaps(regions, promoter))
  if(nrow(idx)>0){
    p <- CoveragePlot(
      object = scATAC.data,
      region = x,
      ranges.title = "MACS2",
      links = F,
      peaks = T,
      extend.upstream = 1000,
      extend.downstream = 1000,
      region.highlight = promoter[idx[,2],])
  }else{
    p <- CoveragePlot(
      object = scATAC.data,
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
#genes <- c("CD8A", "CD4", "FOXP3", "FGFBP2", "MS4A1", "CD79A", "CSF1R", "CD163", "S100A12", "FCGR3A", "TPSAB1", "KDR", "VCAM1", "PDGFRB", "CA9")
pdf("2.Cluster/AnnotateCellType/cellType.coverage.plot2.pdf", width = unit(3, "inches"), height = unit(3, "inches"))
res <- sapply(genes, function(x){
  cat(x, "...\n")
  regions <- FindRegion(object = scATAC.data, region = x, assay = "Peaks", extend.upstream = 1000, extend.downstream = 1000)
  idx <- data.frame(findOverlaps(regions, promoter))
  if(nrow(idx)>0){
    p <- CoveragePlot(
      object = scATAC.data,
      region = x,
      ranges.title = "MACS2",
      links = F,
      peaks = F,
      extend.upstream = 1000,
      extend.downstream = 1000,
      region.highlight = promoter[idx[,2],])
  }else{
    p <- CoveragePlot(
      object = scATAC.data,
      region = x,
      ranges.title = "MACS2",
      links = F,
      extend.upstream = 1000,
      extend.downstream = 1000,
      peaks = F)
  }
  print(p)
  return(regions)
})
dev.off()

# Additional modifications based on the promoter core sequence in ensemble
#### 修改region长度 绘制
pdf("2.Cluster/AnnotateCellType/cellType.coverage.user.pdf", width = unit(3, "inches"), height = unit(3, "inches"))
# CD8A
# CD8 <- FindRegion(object = scATAC.data, region = "CD8A", assay = "Peaks", extend.upstream = 0, extend.downstream = 0)
CD8 <- StringToGRanges("chr2-86784610-86808396")
regions <- FindRegion(object = scATAC.data, region = "CD8A", assay = "Peaks", extend.upstream = 1000, extend.downstream = 200)
idx <- data.frame(findOverlaps(regions, promoter))
regions <- promoter[idx[,2],]
regions <- StringToGRanges(c("chr2-86787200-86793401", "chr2-86805800-86810600")) # core pormoter region from ensembel 
p <- CoveragePlot(
  object = scATAC.data,
  region = CD8,
  links = F,
  peaks = F,
  extend.upstream = 0,
  extend.downstream = 1000,
  region.highlight = regions)
print(p)

# GNLY
# GNLY <- FindRegion(object = scATAC.data, region = "GNLY", assay = "Peaks", extend.upstream = 0, extend.downstream = 0)
GNLY <- StringToGRanges("chr2-85685175-85698854")
# regions <- FindRegion(object = scATAC.data, region = "FGFBP2", assay = "Peaks", extend.upstream = 1000, extend.downstream = 200)
# idx <- data.frame(findOverlaps(regions, promoter))
# regions <- promoter[idx[,2],]
regions <- StringToGRanges("chr2-85693002-85695201") # extend the promoter flank region
p <- CoveragePlot(
  object = scATAC.data,
  region = GNLY,
  links = F,
  peaks = F,
  extend.upstream = 0,
  extend.downstream = 0,
  region.highlight = regions)
print(p)

# CD4
# CD4 <- FindRegion(object = scATAC.data, region = "CD4", assay = "Peaks", extend.upstream = 0, extend.downstream = 0)
CD4 <- StringToGRanges("chr12-6786858-6800000")
# regions <- FindRegion(object = scATAC.data, region = "CD4", assay = "Peaks", extend.upstream = 1000, extend.downstream = 200)
# idx <- data.frame(findOverlaps(regions, promoter))
# regions <- promoter[idx[,2],]
regions <- StringToGRanges(c("chr12-6786600-6787201", "chr12-6788200-6793801"))
p <- CoveragePlot(
  object = scATAC.data,
  region = CD4,
  links = F,
  peaks = F,
  extend.upstream = 0,
  extend.downstream = 0,
  region.highlight = regions)
print(p)

# MS4A1 <- FindRegion(object = scATAC.data, region = "MS4A1", assay = "Peaks", extend.upstream = 0, extend.downstream = 0)
MS4A1 <- StringToGRanges("chr11-60455752-60460000")
regions <- FindRegion(object = scATAC.data, region = "MS4A1", assay = "Peaks", extend.upstream = 1000, extend.downstream = 200)
idx <- data.frame(findOverlaps(regions, promoter))
regions <- promoter[idx[,2],]
#regions <- StringToGRanges(c("chr11-60455600-60458400"))
p <- CoveragePlot(
  object = scATAC.data,
  region = MS4A1,
  links = F,
  peaks = F,
  extend.upstream = 1000,
  extend.downstream = 0,
  region.highlight = regions)
print(p)

# # CSF1R
# # CSF1R <- FindRegion(object = scATAC.data, region = "CSF1R", assay = "Peaks", extend.upstream = 0, extend.downstream = 0)
# CSF1R <- StringToGRanges("chr5-150060000-150090000")
# # regions <- FindRegion(object = scATAC.data, region = "CSF1R", assay = "Peaks", extend.upstream = 1000, extend.downstream = 200)
# # idx <- data.frame(findOverlaps(regions, promoter))
# # regions <- promoter[idx[,2],]
# regions <- StringToGRanges(c("chr5-150081600-150086801"))
# p <- CoveragePlot(
#   object = scATAC.data,
#   region = CSF1R,
#   links = F,
#   peaks = F,
#   extend.upstream = 0,
#   extend.downstream = 0,
#   region.highlight = regions)
# print(p)

# CD163
# CD163 <- FindRegion(object = scATAC.data, region = "CD163", assay = "Peaks", extend.upstream = 0, extend.downstream = 0)
CD163 <- StringToGRanges("chr12-7490000-7503893")
# regions <- FindRegion(object = scATAC.data, region = "CD163", assay = "Peaks", extend.upstream = 1000, extend.downstream = 200)
# idx <- data.frame(findOverlaps(regions, promoter))
# regions <- promoter[idx[,2],]
regions <- StringToGRanges(c("chr12-7502893-7504893"))  # from overlap to promoter regions
p <- CoveragePlot(
  object = scATAC.data,
  region = CD163,
  links = F,
  peaks = F,
  extend.upstream = 0,
  extend.downstream = 1000,
  region.highlight = regions)
print(p)

# S100A12 <- FindRegion(object = scATAC.data, region = "S100A12", assay = "Peaks", extend.upstream = 0, extend.downstream = 0)
S100A12 <- StringToGRanges("chr1-153373706-153375649")
regions <- FindRegion(object = scATAC.data, region = "S100A12", assay = "Peaks", extend.upstream = 1000, extend.downstream = 200)
idx <- data.frame(findOverlaps(regions, promoter))
regions <- promoter[idx[,2],]
p <- CoveragePlot(
  object = scATAC.data,
  region = S100A12,
  links = F,
  peaks = F,
  extend.upstream = 0,
  extend.downstream = 500,
  region.highlight = regions)
print(p)

# TPSAB1 <- FindRegion(object = scATAC.data, region = "TPSAB1", assay = "Peaks", extend.upstream = 0, extend.downstream = 0)
TPSAB1 <- StringToGRanges("chr16-1240696-1242554")
regions <- FindRegion(object = scATAC.data, region = "TPSAB1", assay = "Peaks", extend.upstream = 1000, extend.downstream = 200)
idx <- data.frame(findOverlaps(regions, promoter))
regions <- promoter[idx[,2],]
p <- CoveragePlot(
  object = scATAC.data,
  region = TPSAB1,
  links = F,
  peaks = F,
  extend.upstream = 1000,
  extend.downstream = 0,
  region.highlight = regions)
print(p)

# PECAM1 <- FindRegion(object = scATAC.data, region = "PECAM1", assay = "Peaks", extend.upstream = 0, extend.downstream = 0)
regions <- StringToGRanges(c("chr17-64386600-64393401", "chr17-64411800-64414001"))
PECAM1 <- StringToGRanges("chr17-64385600-64413776")
p <- CoveragePlot(
  object = scATAC.data,
  region = PECAM1,
  links = F,
  peaks = F,
  extend.upstream = 0,
  extend.downstream = 1000,
  region.highlight = regions)
print(p)

# PDGFRB <- FindRegion(object = scATAC.data, region = "PDGFRB", assay = "Peaks", extend.upstream = 0, extend.downstream = 0)
PDGFRB <- StringToGRanges("chr5-150145000-150155872")
# regions <- FindRegion(object = scATAC.data, region = "PDGFRB", assay = "Peaks", extend.upstream = 1000, extend.downstream = 200)
# idx <- data.frame(findOverlaps(regions, promoter))
# regions <- promoter[idx[,2],]
regions <- StringToGRanges("chr5-150154872-150156872") # from overlap to promoter regions
p <- CoveragePlot(
  object = scATAC.data,
  region = PDGFRB,
  links = F,
  peaks = F,
  extend.upstream = 0,
  extend.downstream = 500,
  region.highlight = regions)
print(p)

# CA9 <- FindRegion(object = scATAC.data, region = "CA9", assay = "Peaks", extend.upstream = 0, extend.downstream = 0)
CA9 <- StringToGRanges("chr9-35673856-35678000")
# regions <- FindRegion(object = scATAC.data, region = "CA9", assay = "Peaks", extend.upstream = 1000, extend.downstream = 200)
# idx <- data.frame(findOverlaps(regions, promoter))
# regions <- promoter[idx[,2],]
regions <- StringToGRanges("chr9-35672918-35674918") # from overlap to promoter regions
p <- CoveragePlot(
  object = scATAC.data,
  region = CA9,
  links = F,
  peaks = F,
  extend.upstream = 1000,
  extend.downstream = 0,
  region.highlight = regions)
print(p)
dev.off()


#### NK markers
NK.markers <- c("KLRB1", "KLRD1", "NKG7", "GNLY")
pdf("2.Cluster/AnnotateCellType/NK.marker.coverage.plot.origin.pdf", height = unit(3, "inches"))
res <- sapply(NK.markers, function(x){
  cat(x, "...\n")
  regions <- FindRegion(object = scATAC.data, region = x, assay = "Peaks", extend.upstream = 1000, extend.downstream = 1000)
  idx <- data.frame(findOverlaps(regions, promoter))
  if(nrow(idx)>0){
    p <- CoveragePlot(
      object = scATAC.data,
      region = x,
      ranges.title = "MACS2",
      links = F,
      peaks = T,
      extend.upstream = 1000,
      extend.downstream = 1000,
      region.highlight = promoter[idx[,2],])
  }else{
    p <- CoveragePlot(
      object = scATAC.data,
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

#### endothelia cells
pdf("2.Cluster/AnnotateCellType/endothelia.marker.coverage.plot.origin.pdf", height = unit(3, "inches"))
regions <- StringToGRanges(c("chr17-64386600-64393401", "chr17-64411800-64414001"))
PECAM1 <- StringToGRanges("chr17-64385600-64413776")
p <- CoveragePlot(
  object = scATAC.data,
  region = PECAM1,
  ranges.title = "MACS2",
  links = F,
  peaks = F,
  extend.upstream = 0,
  extend.downstream = 1000,
  region.highlight = regions)
print(p)

############################# identify differentially accessible chromatin regions between celltypes
DefaultAssay(scATAC.data) <- "Peaks"
Idents(scATAC.data) <- scATAC.data$AnnotatedcellType
idents <- as.character(levels(scATAC.data))
cellType.DARs <- FindAllMarkers(scATAC.data, 
                                test.use = 'LR',
                                logfc.threshold=0, 
                                min.pct = 0.05, # often necessary to lower the min.pct threshold
                                latent.vars = "peak_region_fragments")
cf <- ClosestFeature(scATAC.data, regions = rownames(cellType.DARs)) # Find the closest feature to a given set of genomic regions
cellType.DARs <- cbind(cellType.DARs, gene=cf$gene_name, gene_biotype = cf$gene_biotype, type = cf$type, distance=cf$distance)
colnames(cellType.DARs)[6:7] <- c("cellType", "genomicRegion")
pi <- Combined.P.FC(cellType.DARs[,c("avg_log2FC", "p_val_adj")], log10P = F)
cellType.DARs$pi <- pi$pi
saveFormat <- lapply(idents, function(x){
  index <- which(cellType.DARs$cellType == x)
  DARs <- cellType.DARs[index,]
  DARs.up <- DARs %>% filter(avg_log2FC>0) %>% arrange(desc(pi))
  DARs.down <- DARs %>% filter(avg_log2FC<0) %>% arrange(pi)
  DARs <- rbind(DARs.up, DARs.down)
  return(DARs)
})
write.xlsx(saveFormat, file = "4.Peak/celltype.all.DARs.xlsx", sheetName = idents, rowNames = F)
saveRDS(cellType.DARs, file = "4.Peak/cellType.all.DARs.rds")

#require logfc.threshold = 0.25 & p_val_adj < 0.05
cellType.sig.pos.DARs <- cellType.DARs %>% filter(avg_log2FC >=0.25 & p_val_adj < 0.05) %>% arrange(desc(pi)) # 31925 peaks
saveFormat <- lapply(idents, function(x){
  index <- which(cellType.sig.pos.DARs$cellType == x)
  DARs <- cellType.sig.pos.DARs[index,]
  DARs <- DARs %>% arrange(desc(pi))
  return(DARs)
})
names(saveFormat) <- idents
write.xlsx(saveFormat, file = "4.Peak/celltype.sig.pos.DARs.xlsx", sheetName = idents, rowNames = F)
saveRDS(cellType.sig.pos.DARs, file = "4.Peak/cellType.sig.pos.DARs.rds")

#plot--- differentially accessible chromatin regions heatmap
sig.region <- cellType.sig.pos.DARs %>% select(genomicRegion) %>% distinct() 
#average fragment of each peak in each cell type
sig.region.mean <- AverageExpression(scATAC.data, features = sig.region$genomicRegion, assays = "Peaks")
sig.region.mean.scale <- scale(t(sig.region.mean$Peaks))
pdf("4.Peak/cellType.sig.pos.DAR.pdf")
Heatmap(sig.region.mean.scale, name = "z-score", show_column_dend = F, show_row_dend = F, show_column_names = F, row_names_gp = gpar(fontsize = 10), width = unit(10, "cm"), height = unit(8, "cm"))
dev.off()

##plot---DAR distribution
cellType.sig.pos.DARs.ratio <- as.data.frame(table(cellType.sig.pos.DARs$cellType))
cellType.sig.pos.DARs.ratio$Type <- rep("Lymphoid", nrow(cellType.sig.pos.DARs.ratio))
cellType.sig.pos.DARs.ratio$Type[c(4, 6, 12)] <- "Myeloid"
cellType.sig.pos.DARs.ratio$Type[c(3, 7, 8)] <- "Other"
cellType.sig.pos.DARs.ratio$Type[c(9)] <- "Tumor"
cellType.sig.pos.DARs.ratio$Type <- factor(cellType.sig.pos.DARs.ratio$Type, levels = c("Lymphoid", "Myeloid", "Tumor", "Other"))
pdf("4.Peak/cellType.sig.pos.DAR.ratio.pdf")
ggbarplot(cellType.sig.pos.DARs.ratio, x="Var1", y="Freq", fill = "Type", color = "Type",
          sort.by.groups=FALSE, sort.val = "desc", palette = colors <- c("#00A087", "#4DBBD5", "#E64B35", "#3C5488"),#不按组排序
          label = T, xlab = "", ylab = "Number of DAR") + rotate_x_text(60)
dev.off()

############################# cell tpye differentially accessible chromatin and genes
## load DEGs
scRNA.DEGs <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/2.Cluster/AnnotateCellType/cellType.sig.pos.DEGs.rds")
scRNA.DEGs <- scRNA.DEGs %>% filter(avg_log2FC >= 0.25 & p_val_adj < 0.05)
compared.idents <- as.character(levels(scATAC.data))

# calculate the intersection gene between scRNA-seq and scATAC-seq in same cell type
overlap.gene.list <- sapply(compared.idents, function(x){
  idx <- which(scRNA.DEGs$cluster == x)
  DEGs <- scRNA.DEGs$gene[idx]

  idx <- which(cellType.sig.pos.DARs$cellType == x)
  DARs <- unique(cellType.sig.pos.DARs$gene[idx])

  overlap <- intersect(DEGs, DARs)
  return(overlap)
})
saveRDS(overlap.gene.list, file = "4.Peak/DEG.DAR.overlap.genes.rds")

# calculate the differential ration
overlap.ratio <- sapply(names(overlap.gene.list), function(x){
  genes <- overlap.gene.list[[x]]

  #scRNA-seq
  idx <- which(scRNA.DEGs$cluster == x)
  DEGs <- length(scRNA.DEGs$gene[idx])
  DEGs.with.DARs <- length(genes)
  Prop.DEGs.with.DARs <- DEGs.with.DARs/DEGs

  #scATAC-seq
  idx <- which(cellType.sig.pos.DARs$cellType == x)
  DARs <- length(cellType.sig.pos.DARs$genomicRegion[idx])
  DARs.near.DEGs <- length(which(cellType.sig.pos.DARs$gene[idx] %in% genes))
  Prop.DARs.near.DEGs <- DARs.near.DEGs/DARs
  return(c(DEGs, DEGs.with.DARs, Prop.DEGs.with.DARs, DARs, DARs.near.DEGs, Prop.DARs.near.DEGs))
})
overlap.ratio <- t(overlap.ratio)
colnames(overlap.ratio) <- c("DEGs", "DEGs with DARs", "Prop DEGs with DARs", "DARs", "DARs near DEGs", "Prop DARs near DEGs")

#calculate the min max mean sd
res <- apply(overlap.ratio, 2, function(x){
  return(c(min(x), max(x), mean(x), sd(x)))
})
rownames(res) <- c("min", "max", "mean", "sd")
write.xlsx(list(overlap.ratio, res), file = "4.Peak/overlap.ratio.xlsx", sheetName = c("overlap of DEGs and DARs", "Statistics"), rowNames = T)

############################# cell type sepecific peaks
scRNA.data <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/data.merge.pro.rds")
DefaultAssay(scRNA.data) <- "RNA"
# 展示top5
pdf("4.Peak/cell.type.sepecific.peaks.pdf")
res <- sapply(names(saveFormat),function(x){
  cat("Ploting", x, " top 5 peaks\n")
  y <- saveFormat[[x]]
  top5 <- y$genomicRegion[1:5]
  top5.gene <- y$gene[1:5]

  for(i in 1:length(top5)){
    cat("Plot: ", top5[i], "\n")
    p <- CoveragePlot(object = scATAC.data,
                   region = top5[i],
                   assay = "Peaks",
                   annotation = TRUE,
                   extend.upstream = 5000,
                   extend.downstream = 5000,
                   peaks = TRUE)
    p <- p + ggtitle(paste0(x," ", top5[i]))
    print(p)

    idx <- match(top5.gene[i], rownames(scRNA.data))
    if(!is.na(idx)){
      p <- VlnPlot(scRNA.data, pt.size = 0, features = top5.gene[i], group.by = "cellType_low") + xlab("") + NoLegend() + ggtitle(paste0(x," ", top5.gene[i]))
      print(p)
    }
  }
})
dev.off()
