#' @description: analyisis extra data

library(Seurat)
library(harmony)
library(clustree)
library(ggpubr)
library(dplyr)
library(patchwork)
library(vegan)
set.seed(101)
library(future)
library(ComplexHeatmap)
library(circlize)
library(corrplot)
library(cowplot)
library(openxlsx)
plan("multiprocess", workers = 10) 
options(future.globals.maxSize = 100000 * 1024^2) # set 50G RAM
setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/7.TF.analysis/ExtraData")
tumor.TFs <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/5.Motif/Analysis/tumor.specific.TFs.rds")
four.TFs <- c("HOXC5", "VENTX", "ISL1", "OTP")

target.genes.list <- lapply(four.TFs, function(x){
  targets <- read.xlsx(paste0("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/7.TF.analysis/targets/", x, "_targetGenes.xlsx"))
  targets <- unique(targets$to)
  return(targets)
})
names(target.genes.list) <- four.TFs
target.genes.list$AllTFs <- tumor.TFs$Name
TF.names <- names(target.genes.list)

############################################1.Young_2018_Science###############################################
Young.Science.proj.QC <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/Young_2018_Science_scRNA_kidney/Young.Science.proj.QC.rds")
# 72501 cells

#### extract RCC tumor
RCC.object <- subset(Young.Science.proj.QC, subset = PatientDiseaseState %in% c("RCC", "VHL_RCC"))
# remove RCC3 samples
RCC.source <- unique(RCC.object$Source[grep("RCC1|RCC2|VHL", RCC.object$Source)])
RCC.object <- subset(RCC.object, subset = Source %in% RCC.source) # 30815 cells

# mapping the celltype
cellType <- unique(setdiff(RCC.object$Cell_type1, c(NA, "Junk", "MNP PapRCC", "Wilms_tumour", "Wilms_tumour_and_fibroblast", "Private")))
RCC.object <- subset(RCC.object, subset = Cell_type1 %in% cellType)
cellType <- RCC.object$Cell_type1
cellType <- gsub("NK cell 1", "NK cell", cellType)
cellType <- gsub("NK cell 2", "NK cell", cellType)
cellType <- gsub("Proliferating NK cell", "NK cell", cellType)
cellType[which(cellType %in% c("MNP1-like", "MNP2-like"))] <- "Mononuclear phagocyte-like"
cellType[which(cellType %in% c("MNP RCC2", "MNP RCC1", "MNP1", "MNP2", "MNP3"))] <- "Mononuclear phagocyte"
cellType <- gsub("T regulatory", "Treg", cellType)
cellType <- gsub("Mesangial cells", "Mesangial cell", cellType)
cellType <- gsub("Renal_cell_carcinoma", "Tumor", cellType)
cellType <- gsub("Plasmacytoid DC", "Plasmacytoid dendritic cell", cellType)
RCC.object <- AddMetaData(RCC.object, cellType, "cellType")
# re-assign sub cell type
extra.cellType <- c("Normal_cell")
for(i in extra.cellType){
    idx <- which(RCC.object$cellType == i)
    
    index1 <- which(RCC.object$Cell_type2[idx] != "-")

    if(length(index1)>0){
        RCC.object$cellType[idx][index1] <- RCC.object$Cell_type2[idx][index1]
        RCC.object$cellType[idx][index1] <- gsub("_", " ", RCC.object$cellType[idx][index1])
    }
}
RCC.object$cellType <- gsub("_", " ", RCC.object$cellType)

# remove normal cell , Junk and smaller population with less than 50 cells 
cell.Type <- unique(setdiff(RCC.object$cellType, c("Normal cell", "Ureter others", "Erythroblast")))
RCC.object <- subset(RCC.object, subset = cellType %in% cell.Type)

# extract tumor region
RCC.object.tumor <- subset(RCC.object, subset = Category %in% c("Kidney_tumour", "Kidney_tumour_immune"))
RCC.object <- NormalizeData(RCC.object, verbose = FALSE) # 20499 cells
RCC.object.tumor <- NormalizeData(RCC.object.tumor, verbose = FALSE) # 14681 cells

pdf("Young.Science.2018.ccRCCSample.tumorTF.expression.pdf")
DotPlot(RCC.object, cols = c("#1e90ff", "#ff5a36"), features = tumor.TFs$Name, group.by = "cellType", dot.scale = 4) + theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 8))
dev.off()
pdf("Young.Science.2018.ccRCCSample.FourTF.expression.pdf", width = 5)
DotPlot(RCC.object, cols = c("#1e90ff", "#ff5a36"), features = four.TFs, group.by = "cellType", dot.scale = 4) + RotatedAxis() + theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8))
dev.off()

RCC.object <- AddModuleScore(RCC.object, features = target.genes.list, name = TF.names)
colnames(RCC.object@meta.data)[(ncol(RCC.object@meta.data)-4):ncol(RCC.object@meta.data)] <- paste0(TF.names, "_targets")
RCC.object$cellType <- factor(RCC.object$cellType, levels = c(setdiff(names(table(RCC.object$cellType)), "Tumor"), "Tumor"))

signature.scores <- RCC.object@meta.data
pdf("Young.Science.2018.ccRCCSample.targetScore.pdf", height = 4)
p <- lapply(TF.names, function(x){
    p <- ggviolin(signature.scores, y=paste0(x, '_targets'), x='cellType', pt.size=0, fill = 'cellType', add = "mean_sd", error.plot = "crossbar") +
    xlab('') + ylab(paste0(x, ' target score')) + ggtitle('') +
    theme(plot.margin = unit(c(0, 0, 0, 0.1), "in"), axis.title.y=element_text(face='bold'), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    NoLegend() + stat_compare_means(ref.group = "Tumor", label = "p.signif")
    print(p)
    return(p)
})
dev.off()

pdf("Young.Science.2018.ccRCCRegion.tumorTF.expression.pdf")
DotPlot(RCC.object.tumor, cols = c("#1e90ff", "#ff5a36"), features = tumor.TFs$Name, group.by = "cellType", dot.scale = 4) + theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 8))
dev.off()
pdf("Young.Science.2018.ccRCCRegion.FourTF.expression.pdf", width = 5)
DotPlot(RCC.object.tumor, cols = c("#1e90ff", "#ff5a36"), features = four.TFs, group.by = "cellType", dot.scale = 4) + RotatedAxis() + theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8))
dev.off()

RCC.object.tumor <- AddModuleScore(RCC.object.tumor, features = target.genes.list, name = TF.names)
colnames(RCC.object.tumor@meta.data)[(ncol(RCC.object.tumor@meta.data)-4):ncol(RCC.object.tumor@meta.data)] <- paste0(TF.names, "_targets")
RCC.object.tumor$cellType <- factor(RCC.object.tumor$cellType, levels = c(setdiff(names(table(RCC.object.tumor$cellType)), "Tumor"), "Tumor"))

signature.scores <- RCC.object.tumor@meta.data
pdf("Young.Science.2018.ccRCCRegion.targetScore.pdf", height = 4)
p <- lapply(TF.names, function(x){
    p <- ggviolin(signature.scores, y=paste0(x, '_targets'), x='cellType', pt.size=0, fill = 'cellType', add = "mean_sd", error.plot = "crossbar") +
    xlab('') + ylab(paste0(x, ' target score')) + ggtitle('') +
    theme(plot.margin = unit(c(0, 0, 0, 0.1), "in"), axis.title.y=element_text(face='bold'), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    NoLegend() + stat_compare_means(ref.group = "Tumor", label = "p.signif")
    print(p)
    return(p)
})
dev.off()

############################################2.Bi.CancerCell.2021###############################################
Bi.CancerCell.2021 <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/Bi_2021_CancerCell_scRNA_ccRCC/ICB.ccRCC.rds")
Bi.CancerCell.2021 <- subset(Bi.CancerCell.2021, subset = disease__ontology_label == "clear cell renal carcinoma")

# merge cell type
cellType <- Bi.CancerCell.2021$FinalCellType
cellType <- gsub("41BB-Hi CD8\\+ T cell|41BB-Lo CD8\\+ T cell|Cycling CD8\\+ T cell|MitoHigh CD8\\+ T cell|MX1-Hi CD8\\+ T cell", "CD8+ T cell", cellType)
cellType <- gsub("CD16- Monocyte|CD16\\+ Monocyte", "Monocyte", cellType)
cellType <- gsub("CD1C\\+ DC|CLEC9A\\+ DC", "Dendritic cell", cellType)
cellType <- gsub("CXCL10-Hi TAM|Cycling TAM|FOLR2-Hi TAM|GPNMB-Hi TAM|VSIR-Hi TAM", "TAM", cellType)
cellType <- gsub("Memory T-Helper|MitoHigh T-Helper|Effector T-Helper", "T-Helper", cellType)
cellType <- gsub("FGFBP2\\- NK|FGFBP2\\+ NK|MitoHigh NK", "NK", cellType)
cellType <- gsub("TP1|TP2|Cycling Tumor", "Tumor", cellType)
Bi.CancerCell.2021 <- AddMetaData(Bi.CancerCell.2021, cellType, "cellType")
cell.Type <- setdiff(unique(Bi.CancerCell.2021$cellType), c("Misc/Undetermined")) # 31599 cells
Bi.CancerCell.2021 <- subset(Bi.CancerCell.2021, subset = cellType %in% cell.Type)

pdf("Bi.CancerCell.2021.Four.TF.expression.pdf", width = 5)
DotPlot(Bi.CancerCell.2021, cols = c("#1e90ff", "#ff5a36"), features = four.TFs, group.by = "cellType", dot.scale = 4) + RotatedAxis() + theme(axis.text.x = element_text(size = 8)) + theme(axis.text.y = element_text(size = 8))
dev.off()
pdf("Bi.CancerCell.2021.tumor.TF.expression.pdf")
DotPlot(Bi.CancerCell.2021, cols = c("#1e90ff", "#ff5a36"), features = tumor.TFs$Name, group.by = "cellType", dot.scale = 4) + theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.text.y = element_text(size = 8))
dev.off()

Bi.CancerCell.2021 <- AddModuleScore(Bi.CancerCell.2021, target.genes.list, name = TF.names)
colnames(Bi.CancerCell.2021@meta.data)[(ncol(Bi.CancerCell.2021@meta.data)-4):ncol(Bi.CancerCell.2021@meta.data)] <- paste0(TF.names, '_targets')
Bi.CancerCell.2021$cellType <- factor(Bi.CancerCell.2021$cellType, levels = c(setdiff(unique(Bi.CancerCell.2021$cellType), c("Tumor")), c("Tumor")))

signature.scores <- Bi.CancerCell.2021@meta.data
pdf("Bi.CancerCell.2021.targetScore.pdf", height = 4)
p <- lapply(TF.names, function(x){
    p <- ggviolin(signature.scores, y=paste0(x, '_targets'), x='cellType', pt.size=0, fill = 'cellType', add = "mean_sd", error.plot = "crossbar") +
    xlab('') + ylab(paste0(x, ' target score')) + ggtitle('') +
    theme(plot.margin = unit(c(0, 0, 0, 0.1), "in"), axis.title.y=element_text(face='bold'), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    NoLegend() + stat_compare_means(ref.group = "Tumor", label = "p.signif")
    print(p)
    return(p)
})
dev.off()

############################################3.Braun.CancerCell.2021###############################################
Braun.CancerCell.2021 <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/Braun_2021_CancerCell_scRNA_ccRCC/diff.state.ccRCC.rds")
cellType <- Braun.CancerCell.2021@meta.data$ClusterName_AllCells
cellType <- gsub("Tumor cell.*", "Tumor", cellType)
cellType <- gsub("CD8\\+ T cell.*", "CD8+ T cell", cellType)
cellType <- gsub("CD8\\+ Tcell\\.3", "CD8+ T cell", cellType)
cellType <- gsub("CD4\\+ T cell.*", "CD4+ T cell", cellType)
cellType <- gsub("Myeloid cell.*", "Myeloid cell", cellType)
cellType <- gsub("NK cell.*", "NK cell", cellType)
cellType <- gsub("^T cell mixed|Tumor-immune doublet|Immune doublet\\.1|Immune doublet\\.2|Immune doublet\\.3$", "doublet", cellType)
Braun.CancerCell.2021 <- AddMetaData(Braun.CancerCell.2021, cellType, "cellType")
Braun.CancerCell.2021 <- subset(Braun.CancerCell.2021, subset = cellType %in% setdiff(unique(Braun.CancerCell.2021$cellType), "doublet"))

# select tumor sample
Braun.CancerCell.2021.tumor <- subset(Braun.CancerCell.2021, subset = Tumor_Normal == "T")

Braun.CancerCell.2021 <- NormalizeData(Braun.CancerCell.2021, verbose = FALSE) # 157136 cells
Braun.CancerCell.2021.tumor <- NormalizeData(Braun.CancerCell.2021.tumor, verbose = FALSE) # 116171 cells

pdf("Braun.CancerCell.2021.FourTF.expression.pdf", width = 5)
DotPlot(Braun.CancerCell.2021, cols = c("#1e90ff", "#ff5a36"), features = four.TFs, group.by = "cellType", dot.scale = 4) + RotatedAxis() + theme(axis.text.x = element_text(size = 8)) + theme(axis.text.y = element_text(size = 8))
dev.off()
pdf("Braun.CancerCell.2021.FourTF.expression.Vlnplot.pdf", width = 5)
VlnPlot(Braun.CancerCell.2021, features = four.TFs, group.by="cellType", same.y.lims=T,flip = T, stack = T) & xlab("") & ylab("Log-normalized expression") & theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.5))
dev.off()
pdf("Braun.CancerCell.2021.tumorTF.expression.pdf")
DotPlot(Braun.CancerCell.2021, cols = c("#1e90ff", "#ff5a36"), features = tumor.TFs$Name, group.by = "cellType", dot.scale = 4) + theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.text.y = element_text(size = 8))
dev.off()

Braun.CancerCell.2021 <- AddModuleScore(Braun.CancerCell.2021, target.genes.list, name = TF.names)
colnames(Braun.CancerCell.2021@meta.data)[(ncol(Braun.CancerCell.2021@meta.data)-4):ncol(Braun.CancerCell.2021@meta.data)] <- paste0(TF.names, '_targets')
Braun.CancerCell.2021$cellType <- factor(Braun.CancerCell.2021$cellType, levels = c(setdiff(names(table(Braun.CancerCell.2021$cellType)), "Tumor cell"), "Tumor cell"))

signature.scores <- Braun.CancerCell.2021@meta.data
pdf("Braun.CancerCell.2021.ccRCCSample.targetScore.pdf", height = 4)
p <- lapply(TF.names, function(x){
    p <- ggviolin(signature.scores, y=paste0(x, '_targets'), x='cellType', pt.size=0, fill = 'cellType', add = "mean_sd", error.plot = "crossbar") +
    xlab('') + ylab(paste0(x, ' target score')) + ggtitle('') +
    theme(plot.margin = unit(c(0, 0, 0, 0.1), "in"), axis.title.y=element_text(face='bold'), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    NoLegend() + stat_compare_means(ref.group = "Tumor", label = "p.signif")
    print(p)
    return(p)
})
dev.off()

# Tumor region
Braun.CancerCell.2021.tumor <- AddModuleScore(Braun.CancerCell.2021.tumor, target.genes.list, name = TF.names)
colnames(Braun.CancerCell.2021.tumor@meta.data)[(ncol(Braun.CancerCell.2021.tumor@meta.data)-4):ncol(Braun.CancerCell.2021.tumor@meta.data)] <- paste0(TF.names, '_targets')
Braun.CancerCell.2021.tumor$cellType <- factor(Braun.CancerCell.2021.tumor$cellType, levels = c(setdiff(unique(Braun.CancerCell.2021.tumor$cellType), "Tumor"), "Tumor"))

signature.scores <- Braun.CancerCell.2021.tumor@meta.data
pdf("Braun.CancerCell.2021.ccRCCRegion.targetScore.pdf", height = 4)
p <- lapply(TF.names, function(x){
    p <- ggviolin(signature.scores, y=paste0(x, '_targets'), x='cellType', pt.size=0, fill = 'cellType', add = "mean_sd", error.plot = "crossbar") +
    xlab('') + ylab(paste0(x, ' target score')) + ggtitle('') +
    theme(plot.margin = unit(c(0, 0, 0, 0.1), "in"), axis.title.y=element_text(face='bold'), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    NoLegend() + stat_compare_means(ref.group = "Tumor", label = "p.signif")
    print(p)
    return(p)
})
dev.off()

############################################4.CORCES_2018_Science_ATACseq_TCGA###############################################
#### motif TFs
cluster.TFs <- read.table("/data/active_data/lzl/RenalTumor-20200713/Data/CORCES_2018_Science_ATACseq_TCGA/DataS4.TF.enriched.in.clusteredPeak.txt", header = T, stringsAsFactors = F, sep = "\t")
rownames(cluster.TFs) <- cluster.TFs$TF_Motif_Name

tumor.cluster.TFs <- cluster.TFs[which(cluster.TFs$TF_Name %in% tumor.TFs$Name),]
four.cluster.TFs <- cluster.TFs[which(cluster.TFs$TF_Name %in% four.TFs),]
tumor.cluster.TFs <- tumor.cluster.TFs[, -c(1:4)]
four.cluster.TFs <- four.cluster.TFs[, -c(1:4)]

col_fun <- colorRamp2(c(0, 250), c("white", "red"))
pdf("CORCES_2018_Science_ATACseq_TCGA.tumorTFs.heatmap.pdf")
Heatmap(tumor.cluster.TFs, name = "-log10(P value)",  col = col_fun,
        width = unit(10, "cm"), show_row_dend = F, show_column_dend = F,
        column_names_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 5))
dev.off()

pdf("CORCES_2018_Science_ATACseq_TCGA.fourTFs.heatmap.pdf")
Heatmap(four.cluster.TFs, name = "-log10(P value)",  col = col_fun,
        width = unit(10, "cm"), height = unit(5, "cm"), show_row_dend = F, show_column_dend = F,
        column_names_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 10))
dev.off()

# adjust the show 
idx <- c("ISL1_617", "HOXC5_847", "HOXC5_805", "VENTX_748", "OTP_804")
index <- match(idx, cluster.TFs$TF_Motif_Name)
four.cluster.TFs <- cluster.TFs[index,]
rownames(four.cluster.TFs) <- paste0(four.cluster.TFs$TF_Name, " (", four.cluster.TFs$CIS.BP_ID, ")")
four.cluster.TFs <- four.cluster.TFs[, -c(1:4)]
pdf("CORCES_2018_Science_ATACseq_TCGA.fourTFs.heatmap.fined.pdf")
Heatmap(four.cluster.TFs, name = "-log10(P value)",  col = col_fun, 
        column_split = colnames(four.cluster.TFs), row_split = rownames(four.cluster.TFs),
        row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE, border_gp = gpar(col = "grey"),
        width = unit(10, "cm"), height = unit(4, "cm"), show_row_dend = F, show_column_dend = F,
        column_names_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 10))
dev.off()

#### Peak
library(Signac)
library(GenomicRanges)
panCancer.peaks <- read.table("/data/active_data/lzl/RenalTumor-20200713/Data/CORCES_2018_Science_ATACseq_TCGA/TCGA-ATAC_PanCancer_PeakSet.txt", header = T, stringsAsFactors = F, sep = "\t")
# 562709 peaks
panCancer.peaks.gr <- makeGRangesFromDataFrame(panCancer.peaks)
panCancer.peaks.gr$score <- panCancer.peaks$score
panCancer.peaks.gr$annotation <- panCancer.peaks$annotation
panCancer.peaks.gr$percentGC <- panCancer.peaks$percentGC
panCancer.peaks.gr$name <- panCancer.peaks$name
panCancer.peaks.totalNumber <- as.data.frame(table(gsub("_.*", "", panCancer.peaks.gr$name)))

cellType.sig.pos.DARs <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/4.Peak/cellType.sig.pos.DARs.rds")
cellType.gr.overlap <- lapply(as.character(unique(cellType.sig.pos.DARs$cellType)), function(x){
    idx <- which(cellType.sig.pos.DARs$cellType == x)
    cellType.DARs <- cellType.sig.pos.DARs[idx,]
    gr <- StringToGRanges(cellType.DARs$genomicRegion)
    gr$gene <- cellType.DARs$gene
    gr$gene_biotype <- cellType.DARs$gene_biotype
    gr$p_val_adj <- cellType.DARs$p_val_adj
    gr$avg_log2FC <- cellType.DARs$avg_log2FC
    gr$name <- cellType.DARs$cellType
    overlap <- findOverlaps(panCancer.peaks.gr, gr)
    return(list(gr = gr, overlap = overlap))
})
names(cellType.gr.overlap) <- as.character(unique(cellType.sig.pos.DARs$cellType))

percentage.ratio <- sapply(cellType.gr.overlap, function(x){
    panCancer.ratio <- panCancer.peaks.gr[unique(x$overlap@from),]
    panCancer.cellType.totalNumber <- as.data.frame(table(gsub("_.*", "", panCancer.ratio$name)))
    panCancer.cellType.totalNumber$ratio <- panCancer.cellType.totalNumber$Freq/panCancer.peaks.totalNumber$Freq
    return(panCancer.cellType.totalNumber$ratio)
})
rownames(percentage.ratio) <- panCancer.peaks.totalNumber$Var1
# re-order the column
percentage.ratio <- percentage.ratio[, c(7,11,8,5,6,3,9,12,10,1,2,4)]
col_fun <- colorRamp2(c(0, 0.12), c("white", "blue"))
pdf("CORCES_2018_Science_ATACseq_TCGA.peak.heatmap.pdf")
Heatmap(percentage.ratio, name = "Percentage",  col = col_fun, 
        column_split = factor(colnames(percentage.ratio), levels = colnames(percentage.ratio)), row_split = rownames(percentage.ratio),
        row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE, border_gp = gpar(col = "grey"),
        width = unit(6, "cm"), height = unit(8, "cm"), show_row_dend = F, show_column_dend = F,
        cluster_rows = F, cluster_columns = F,
        column_names_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 10))
dev.off()

# permutation test
library(parallel)
peak.number <- as.data.frame(table(cellType.sig.pos.DARs$cellType))
cellType.peak.number <- peak.number$Freq
names(cellType.peak.number) <- peak.number$Var1
permutated.overlap <- function(N, percentage.ratio, cellType.gr.overlap, panCancer.peaks.totalNumber, cellType.sig.pos.DARs, cellType.peak.number){
    
    random.res <- lapply(N, function(x){
        res <- sapply(colnames(percentage.ratio), function(cellType){
            cat("Calculating cell type:", cellType, "\n")
            celltype.overlap <- percentage.ratio[,cellType]
            peak.num <- cellType.peak.number[cellType]
            new.peaks <- cellType.sig.pos.DARs[sample(1:nrow(cellType.sig.pos.DARs), size = peak.num),]

            gr <- StringToGRanges(new.peaks$genomicRegion)
            gr$name <- cellType
            overlap <- findOverlaps(panCancer.peaks.gr, gr)

            panCancer.ratio <- panCancer.peaks.gr[unique(overlap@from),]
            panCancer.cellType.totalNumber <- as.data.frame(table(gsub("_.*", "", panCancer.ratio$name)))
            panCancer.cellType.totalNumber$ratio <- panCancer.cellType.totalNumber$Freq/panCancer.peaks.totalNumber$Freq
            return(panCancer.cellType.totalNumber$ratio)
        })
        rownames(res) <- panCancer.peaks.totalNumber$Var1
        return(res)
    })
    
    if(class(random.res) == "try-error"){
        return("failed random!")
    }else{
        return(random.res)
    }
}

library(parallel)
mc <- getOption("mc.cores", 48)
set.seed(101)
permutated.results <- mclapply(X = 1:1000, 
                FUN = permutated.overlap, 
                percentage.ratio = percentage.ratio, 
                cellType.gr.overlap = cellType.gr.overlap, 
                panCancer.peaks.totalNumber = panCancer.peaks.totalNumber, 
                cellType.sig.pos.DARs = cellType.sig.pos.DARs, 
                cellType.peak.number = cellType.peak.number, 
                mc.cores = mc)
saveRDS(permutated.results, file = "/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/7.TF.analysis/ExtraData/CORCES_2018_Science_ATACseq_TCGA.permutated.results.rds")

# test
p.matrix <- matrix(0, ncol = ncol(percentage.ratio), nrow = nrow(percentage.ratio))
for(i in 1:nrow(percentage.ratio)){
    for(j in 1:ncol(percentage.ratio)){
        actual.value <- percentage.ratio[i,j]
        random.value  <- sapply(permutated.results, function(x){
            return(x[[1]][i,j])
        })
        n <- length(which(random.value<=actual.value)) # H0: μ<=μ0
        res <- binom.test(x = n, n = 1000, p = 0.5, alternative = "greater") # H1:μ>μ0
        p.matrix[i,j] <- res$p.value
    }
}
adj.p.matrix <- matrix(p.adjust(p.matrix, method = "fdr"), ncol = ncol(percentage.ratio), nrow = nrow(percentage.ratio))
rownames(adj.p.matrix) <- rownames(percentage.ratio)
colnames(adj.p.matrix) <- colnames(percentage.ratio)
adj.p.matrix <- -log10(adj.p.matrix)
col_fun <- colorRamp2(c(0, 300), c("white", "red"))
pdf("CORCES_2018_Science_ATACseq_TCGA.peak.heatmap.fdr.pdf")
Heatmap(adj.p.matrix, name = "-log10(FDR)", col = col_fun,
        column_split = factor(colnames(adj.p.matrix), levels = colnames(adj.p.matrix)), row_split = rownames(adj.p.matrix),
        row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE, border_gp = gpar(col = "grey"),
        width = unit(6, "cm"), height = unit(8, "cm"), show_row_dend = F, show_column_dend = F,
        cluster_rows = F, cluster_columns = F,
        column_names_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 10))
dev.off()
