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

# extract ICB and non-ICB
Bi.CancerCell.2021.ICB <- subset(Bi.CancerCell.2021, subset = ICB_Exposed == "ICB") # 14971 cells
Bi.CancerCell.2021.NoICB <- subset(Bi.CancerCell.2021, subset = ICB_Exposed == "NoICB") # 16628 cells

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

# ICB
Bi.CancerCell.2021.ICB <- AddModuleScore(Bi.CancerCell.2021.ICB, target.genes.list, name = TF.names)
colnames(Bi.CancerCell.2021.ICB@meta.data)[(ncol(Bi.CancerCell.2021.ICB@meta.data)-4):ncol(Bi.CancerCell.2021.ICB@meta.data)] <- paste0(TF.names, '_targets')
Bi.CancerCell.2021.ICB$cellType <- factor(Bi.CancerCell.2021.ICB$cellType, levels = c(setdiff(unique(Bi.CancerCell.2021.ICB$cellType), c("Tumor")), c("Tumor")))
signature.scores <- Bi.CancerCell.2021.ICB@meta.data
pdf("Bi.CancerCell.2021.ICB.targetScore.pdf", height = 4)
p <- lapply(TF.names, function(x){
    p <- ggviolin(signature.scores, y=paste0(x, '_targets'), x='cellType', pt.size=0, fill = 'cellType', add = "mean_sd", error.plot = "crossbar") +
    xlab('') + ylab(paste0(x, ' target score')) + ggtitle('') +
    theme(plot.margin = unit(c(0, 0, 0, 0.1), "in"), axis.title.y=element_text(face='bold'), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    NoLegend() + stat_compare_means(ref.group = "Tumor", label = "p.signif")
    print(p)
    return(p)
})
dev.off()

# No ICB
Bi.CancerCell.2021.NoICB <- AddModuleScore(Bi.CancerCell.2021.NoICB, target.genes.list, name = TF.names)
colnames(Bi.CancerCell.2021.NoICB@meta.data)[(ncol(Bi.CancerCell.2021.NoICB@meta.data)-4):ncol(Bi.CancerCell.2021.NoICB@meta.data)] <- paste0(TF.names, '_targets')
Bi.CancerCell.2021.NoICB$cellType <- factor(Bi.CancerCell.2021.NoICB$cellType, levels = c(setdiff(unique(Bi.CancerCell.2021.NoICB$cellType), c("Tumor")), c("Tumor")))
signature.scores <- Bi.CancerCell.2021.NoICB@meta.data
pdf("Bi.CancerCell.2021.NoICB.targetScore.pdf", height = 4)
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

############################################4.Krishna.CancerCell.2021###############################################
Krishna.CancerCell.2021 <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/Krishna_2021_CancerCell_scRNA_ccRCC/ccRCC_6pat_Seurat.rds")
Krishna.CancerCell.2021 <- UpdateSeuratObject(Krishna.CancerCell.2021)
Krishna.CancerCell.2021 <- NormalizeData(Krishna.CancerCell.2021, verbose = FALSE)
# change the ENSG to gene name
source("/home/longzhilin/Analysis_Code/IDConvert.R")
geneName <- IDConvert(genes = rownames(Krishna.CancerCell.2021@assays$RNA@data), method = "clusterProfiler", fromType = "ENSEMBL", toType = "SYMBOL")
idx <- match(geneName$ENSEMBL, rownames(Krishna.CancerCell.2021@assays$RNA@data))
Krishna.CancerCell.2021@assays$RNA@data <- Krishna.CancerCell.2021@assays$RNA@data[idx,]
rownames(Krishna.CancerCell.2021@assays$RNA@data) <- geneName$SYMBOL
Krishna.CancerCell.2021@assays$RNA@counts <- Krishna.CancerCell.2021@assays$RNA@counts[idx,]
rownames(Krishna.CancerCell.2021@assays$RNA@counts) <- geneName$SYMBOL

meta.data <- read.table("/data/active_data/lzl/RenalTumor-20200713/Data/Krishna_2021_CancerCell_scRNA_ccRCC/ccRCC_6pat_cell_annotations.txt", header = T, stringsAsFactors = F, sep = "\t")
idx <- match(rownames(Krishna.CancerCell.2021@meta.data), meta.data$cell)
Krishna.CancerCell.2021@meta.data$cellType <- meta.data$cluster_name[idx]
Krishna.CancerCell.2021@meta.data$region <- meta.data$region[idx]
Krishna.CancerCell.2021@meta.data$Sample <- meta.data$Sample[idx]
Krishna.CancerCell.2021@meta.data$Sample2 <- meta.data$Sample2[idx]

# cell type
refinedCellType <- Krishna.CancerCell.2021@meta.data$cellType
refinedCellType <- gsub("^TAM HLAhi|TAM HLAint|TAM ISGhi|TAM ISGint$", "TAM", refinedCellType)
refinedCellType <- gsub("^CD14\\+ Monocyte|CD14\\+/CD16\\+ Monocyte$", "Monocyte", refinedCellType)
refinedCellType <- gsub("^CD4\\+ Activated IEG|CD4\\+ Effector|CD4\\+ Naive|CD4\\+ Proliferating|CD4\\+ Treg$", "CD4+ T cell", refinedCellType)
refinedCellType <- gsub("^CD8A\\+ Exhausted|CD8A\\+ Exhausted IEG|CD8A\\+ Proliferating|CD8A\\+ Tissue-resident$", "CD8+ T cell", refinedCellType)
refinedCellType <- gsub("^CD8A\\+ NK-like$", "NK-like T cell", refinedCellType)
refinedCellType <- gsub("^cDC1|cDC2|pDC", "DC", refinedCellType)
refinedCellType <- gsub("^NK HSP\\+|Conventional NK", "NK cell", refinedCellType)
Krishna.CancerCell.2021 <- AddMetaData(Krishna.CancerCell.2021, refinedCellType, "refinedCellType")
# remove Ambiguous
Krishna.CancerCell.2021 <- subset(Krishna.CancerCell.2021, subset = refinedCellType %in% setdiff(unique(Krishna.CancerCell.2021$refinedCellType), c("TAM/TCR (Ambiguos)", "Ambiguous", "Ambiguous/Dead"))) # 153501 cells

# only extract the tumor region
Krishna.CancerCell.2021.tumor <- subset(Krishna.CancerCell.2021, subset = region %in% c("Center", "LowerMedial")) # 27007 cells

pdf("Krishna.CancerCell.2021.Four.TF.expression.pdf", width = 5)
DotPlot(Krishna.CancerCell.2021, cols = c("#1e90ff", "#ff5a36"), features = four.TFs, group.by = "refinedCellType", dot.scale = 4) + RotatedAxis() + theme(axis.text.x = element_text(size = 8)) + theme(axis.text.y = element_text(size = 8))
dev.off()
pdf("Krishna.CancerCell.2021.Four.TF.expression.Vlnplot.pdf", width = 5)
VlnPlot(Krishna.CancerCell.2021, features = four.TFs, group.by="refinedCellType", same.y.lims=T,flip = T, stack = T) & xlab("") & ylab("Log-normalized expression") & theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.5))
dev.off()
pdf("Krishna.CancerCell.2021.tumor.TF.expression.pdf")
DotPlot(Krishna.CancerCell.2021, cols = c("#1e90ff", "#ff5a36"), features = tumor.TFs$Name, group.by = "refinedCellType", dot.scale = 4) + theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.text.y = element_text(size = 8))
dev.off()

Krishna.CancerCell.2021 <- AddModuleScore(Krishna.CancerCell.2021, target.genes.list, name = TF.names)
colnames(Krishna.CancerCell.2021@meta.data)[(ncol(Krishna.CancerCell.2021@meta.data)-4):ncol(Krishna.CancerCell.2021@meta.data)] <- paste0(TF.names, '_targets')
Krishna.CancerCell.2021$refinedCellType <- factor(Krishna.CancerCell.2021$refinedCellType, levels = c(setdiff(names(table(Krishna.CancerCell.2021$refinedCellType)), c("CD45- Myofibroblast", "CD45- Vascular Endothelium", "CD45- PAX8+ renal epithelium", "CD45- ccRCC CA9+")), c("CD45- Myofibroblast", "CD45- PAX8+ renal epithelium", "CD45- Vascular Endothelium", "CD45- PAX8+ renal epithelium", "CD45- ccRCC CA9+")))
signature.scores <- Krishna.CancerCell.2021@meta.data
pdf("Krishna.CancerCell.2021.ccRCCSample.targetScore.pdf", height = 4)
p <- lapply(TF.names, function(x){
    p <- ggviolin(signature.scores, y=paste0(x, '_targets'), x='refinedCellType', pt.size=0, fill = 'refinedCellType', add = "mean_sd", error.plot = "crossbar") +
    xlab('') + ylab(paste0(x, ' target score')) + ggtitle('') +
    theme(plot.margin = unit(c(0, 0, 0, 0.1), "in"), axis.title.y=element_text(face='bold'), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    NoLegend() + stat_compare_means(ref.group = "CD45- ccRCC CA9+", label = "p.signif")
    print(p)
    return(p)
})
dev.off()

# only analysis the tumor region
pdf("Krishna.CancerCell.2021.tumor.FourTF.expression.pdf", width = 5)
DotPlot(Krishna.CancerCell.2021.tumor, cols = c("#1e90ff", "#ff5a36"), features = four.TFs, group.by = "refinedCellType", dot.scale = 4) + RotatedAxis() + theme(axis.text.x = element_text(size = 8)) + theme(axis.text.y = element_text(size = 8))
dev.off()
pdf("Krishna.CancerCell.2021.tumor.FourTF.expression.Vlnplot.pdf", width = 5)
VlnPlot(Krishna.CancerCell.2021.tumor, features = four.TFs, group.by="refinedCellType", same.y.lims=T,flip = T, stack = T) & xlab("") & ylab("Log-normalized expression") & theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.5))
dev.off()
pdf("Krishna.CancerCell.2021.tumorTFs.expression.pdf")
DotPlot(Krishna.CancerCell.2021.tumor, cols = c("#1e90ff", "#ff5a36"), features = tumor.TFs$Name, group.by = "refinedCellType", dot.scale = 4) + theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.text.y = element_text(size = 8))
dev.off()

Krishna.CancerCell.2021.tumor <- AddModuleScore(Krishna.CancerCell.2021.tumor, target.genes.list, name = TF.names)
colnames(Krishna.CancerCell.2021.tumor@meta.data)[(ncol(Krishna.CancerCell.2021.tumor@meta.data)-4):ncol(Krishna.CancerCell.2021.tumor@meta.data)] <- paste0(TF.names, '_targets')
Krishna.CancerCell.2021.tumor$refinedCellType <- factor(Krishna.CancerCell.2021.tumor$refinedCellType, levels = c(setdiff(names(table(Krishna.CancerCell.2021.tumor$refinedCellType)), c("CD45- Myofibroblast", "CD45- Vascular Endothelium", "CD45- PAX8+ renal epithelium", "CD45- ccRCC CA9+")), c("CD45- Myofibroblast", "CD45- Vascular Endothelium", "CD45- PAX8+ renal epithelium", "CD45- ccRCC CA9+")))
signature.scores <- Krishna.CancerCell.2021.tumor@meta.data
pdf("Krishna.CancerCell.2021.ccRCCRegion.targetScore.pdf", height = 4)
p <- lapply(TF.names, function(x){
    p <- ggviolin(signature.scores, y=paste0(x, '_targets'), x='refinedCellType', pt.size=0, fill = 'refinedCellType', add = "mean_sd", error.plot = "crossbar") +
    xlab('') + ylab(paste0(x, ' target score')) + ggtitle('') +
    theme(plot.margin = unit(c(0, 0, 0, 0.1), "in"), axis.title.y=element_text(face='bold'), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    NoLegend() + stat_compare_means(ref.group = "CD45- ccRCC CA9+", label = "p.signif")
    print(p)
    return(p)
})
dev.off()

############################################5.Beuselinck_2015_ClinCancerRes_bulk_ccRCC###############################################
expMatrix <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/Beuselinck_2015_ClinCancerRes_bulk_ccRCC/norm_exprs.avg.rds")
clin.info <- read.table("/data/active_data/lzl/RenalTumor-20200713/Data/Beuselinck_2015_ClinCancerRes_bulk_ccRCC/E-MTAB-3267.sdrf.txt", header = T, stringsAsFactors = F, sep = "\t")
clin.info <- clin.info[which(clin.info$Characteristics.disease. == "Tumor"),]
clin.info$Sample <- paste0(clin.info$Assay.Name, ".CEL")
idx <- match(colnames(expMatrix), clin.info$Sample)
clin.info <- clin.info[idx,]

# survival analysis
source("/home/longzhilin/Analysis_Code/SurvivalAnalysis/plot.surv.R")
pdf("Beuselinck_2015_ClinCancerRes.PFS.pdf")
survival.res <- sapply(four.TFs, function(x){
    exp.values <- expMatrix[x,]
    med.value <- median(exp.values)
    group.label <- rep("Low group", length(exp.values))
    idx <- which(exp.values>=med.value)
    group.label[idx] <- "High group"
    PFS.data <- data.frame(Patient_ID = clin.info$Sample, event = clin.info$Characteristics.progression., time = clin.info$Characteristics.progression.free.survival., sample.label = group.label)
    p1 <- plot.surv(PFS.data, risk.table = T, HR = T, ylab = "Progression Free Survival", main = x, surv.median.line = "hv", xlab = "Time (Month)")
    print(p1)
    return(PFS.data)
})
dev.off()

# tumor TF signature & target genes
pdf("Beuselinck_2015_ClinCancerRes.tumorTFs&targetGenes.PFS.pdf")
res <- lapply(target.genes.list, function(x){
    idx <- match(x, rownames(expMatrix))
    interest.matrix <- expMatrix[na.omit(idx),]
    signature.score <- colSums(interest.matrix)/nrow(interest.matrix)
    median.score <- median(signature.score)
    high.group <- which(signature.score >= median.score)
    sample.label <- rep("Low group", length(signature.score))
    sample.label[high.group] <- "High group"
    PFS.data <- data.frame(Patient_ID = clin.info$Sample, event = clin.info$Characteristics.progression., time = clin.info$Characteristics.progression.free.survival., sample.label = sample.label)
    p1 <- plot.surv(PFS.data, risk.table = T, HR = T, ylab = "Progression Free Survival", main = "", surv.median.line = "hv", xlab = "Time (Month)")
    print(p1)
})
dev.off()

############################################6.ICB therapy###############################################
normalized_expression <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/ICB.therapy/normalized_expression.log2.rds")
clin.info <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/ICB.therapy/patient.info.RNA.rds")

# four TF
pdf("ICB.four.TFs.survival.pdf")
survival.res <- lapply(four.TFs, function(x){
    # all cohort
    exp.values <- as.numeric(normalized_expression[x,])
    med.value <- median(exp.values)
    group.label <- rep("Low group", length(exp.values))
    idx <- which(exp.values >= med.value)
    group.label[idx] <- "High group"
    OS.data <- data.frame(Patient_ID = clin.info$SUBJID, event = clin.info$OS_CNSR, time = clin.info$OS, sample.label = group.label)
    p <- plot.surv(OS.data, risk.table = T, HR = T, ylab = "Overall Survival", main = x, surv.median.line = "hv", xlab = "Time (Month)")
    print(p)

    PFS.data <- data.frame(Patient_ID = clin.info$SUBJID, event = clin.info$PFS_CNSR, time = clin.info$PFS, sample.label = group.label)
    p <- plot.surv(PFS.data, risk.table = T, HR = T, ylab = "Progression Free Survival", main = x, surv.median.line = "hv", xlab = "Time (Month)")
    print(p)

    # correlated with immune phenotype
    phenotypes <- c("Angio", "Teff", "Myeloid", "Javelin", "Merck18")
    res <- sapply(phenotypes, function(phenotype){
        a <- data.frame(expression = exp.values, phenotype = clin.info[,phenotype])
        colnames(a) <- c(x, phenotype)
        p <- ggscatter(a, x = x, y = phenotype, add = "reg.line", conf.int = TRUE, xlab = "normalized expression", title = x,
                       add.params = list(color = "blue", fill = "lightgray")) + stat_cor(method = "pearson", label.x = mean(a[,1]), label.y = (max(a[,2])-5))
        print(p)
    })

    # split by cohort
    cohort <- names(table(clin.info$Cohort))
    res <- lapply(cohort, function(y){
        index <- which(clin.info$Cohort == y)
        exp.values <- as.numeric(normalized_expression[x,index])
        med.value <- median(exp.values)
        group.label <- rep("Low group", length(exp.values))
        idx <- which(exp.values >= med.value)
        group.label[idx] <- "High group"

        phenotypes <- c("Angio", "Teff", "Myeloid", "Javelin", "Merck18")
        res <- sapply(phenotypes, function(phenotype){
            a <- data.frame(expression = exp.values, phenotype = clin.info[index,phenotype])
            colnames(a) <- c(x, phenotype)
            p <- ggscatter(a, x = x, y = phenotype, add = "reg.line", conf.int = TRUE, xlab = "normalized expression", title = paste0(x," in ", y),
                        add.params = list(color = "blue", fill = "lightgray")) + stat_cor(method = "pearson", label.x = mean(a[,1]), label.y = (max(a[,2])-5))
            print(p)
        })

        if(length(idx) > length(exp.values)*0.9){
            cat("less sample across each group!\n")
            return("error")
        }else{
            OS.data <- data.frame(Patient_ID = clin.info$SUBJID[index], event = clin.info$OS_CNSR[index], time = clin.info$OS[index], sample.label = group.label)
            p <- plot.surv(OS.data, risk.table = T, HR = T, ylab = "Overall Survival", main = paste0(x," in ", y), surv.median.line = "hv", xlab = "Time (Month)")
            print(p)

            PFS.data <- data.frame(Patient_ID = clin.info$SUBJID[index], event = clin.info$PFS_CNSR[index], time = clin.info$PFS[index], sample.label = group.label)
            p <- plot.surv(PFS.data, risk.table = T, HR = T, ylab = "Progression Free Survival", main = paste0(x," in ", y), surv.median.line = "hv", xlab = "Time (Month)")
            print(p)
            return("ok")
        }
    })
    return(PFS.data)
})
dev.off()

# tumor TF signature & target genes
pdf("ICB.tumorTFs&targetGenes.survival.pdf")
res <- lapply(names(target.genes.list), function(x){
    idx <- match(target.genes.list[[x]], rownames(normalized_expression))
    interest.matrix <- normalized_expression[na.omit(idx),]
    signature.score <- colSums(interest.matrix)/nrow(interest.matrix)
    median.score <- median(signature.score)
    high.group <- which(signature.score >= median.score)
    sample.label <- rep("Low group", length(signature.score))
    sample.label[high.group] <- "High group"

    OS.data <- data.frame(Patient_ID = clin.info$SUBJID, event = clin.info$OS_CNSR, time = clin.info$OS, sample.label = sample.label)
    p <- plot.surv(OS.data, risk.table = T, HR = T, ylab = "Overall Survival", main = x, surv.median.line = "hv", xlab = "Time (Month)")
    print(p)

    PFS.data <- data.frame(Patient_ID = clin.info$SUBJID, event = clin.info$PFS_CNSR, time = clin.info$PFS, sample.label = sample.label)
    p <- plot.surv(PFS.data, risk.table = T, HR = T, ylab = "Progression Free Survival", main = x, surv.median.line = "hv", xlab = "Time (Month)")
    print(p)

    # correlated with immune phenotype
    phenotypes <- c("Angio", "Teff", "Myeloid", "Javelin", "Merck18")
    res <- sapply(phenotypes, function(phenotype){
        a <- data.frame(expression = signature.score, phenotype = clin.info[,phenotype])
        colnames(a) <- c(x, phenotype)
        p <- ggscatter(a, x = x, y = phenotype, add = "reg.line", conf.int = TRUE, xlab = "Targeted gene score", title = x,
                       add.params = list(color = "blue", fill = "lightgray")) + stat_cor(method = "pearson", label.x = mean(a[,1]), label.y = (max(a[,2])-5))
        print(p)
    })

    # split by cohort
    cohort <- names(table(clin.info$Cohort))
    res <- lapply(cohort, function(y){
        index <- which(clin.info$Cohort == y)
        exp.values <- interest.matrix[,index]
        signature.score <- colSums(exp.values)/nrow(exp.values)
        median.score <- median(signature.score)
        high.group <- which(signature.score >= median.score)
        sample.label <- rep("Low group", length(signature.score))
        sample.label[high.group] <- "High group"
        
        res <- sapply(phenotypes, function(phenotype){
            a <- data.frame(expression = signature.score, phenotype = clin.info[index,phenotype])
            colnames(a) <- c(x, phenotype)
            p <- ggscatter(a, x = x, y = phenotype, add = "reg.line", conf.int = TRUE, xlab = "Targeted gene score", title = paste0(x," in ", y),
                        add.params = list(color = "blue", fill = "lightgray")) + stat_cor(method = "pearson", label.x = mean(a[,1]), label.y = (max(a[,2])-5))
            print(p)
        }) 

        if(length(high.group) > length(signature.score)*0.9){
            cat("less sample across each group!\n")
            return("error")
        }else{
            OS.data <- data.frame(Patient_ID = clin.info$SUBJID[index], event = clin.info$OS_CNSR[index], time = clin.info$OS[index], sample.label = sample.label)
            p <- plot.surv(OS.data, risk.table = T, HR = T, ylab = "Overall Survival", main = paste0(x," in ", y), surv.median.line = "hv", xlab = "Time (Month)")
            print(p)

            PFS.data <- data.frame(Patient_ID = clin.info$SUBJID[index], event = clin.info$PFS_CNSR[index], time = clin.info$PFS[index], sample.label = sample.label)
            p <- plot.surv(PFS.data, risk.table = T, HR = T, ylab = "Progression Free Survival", main = paste0(x," in ", y), surv.median.line = "hv", xlab = "Time (Month)")
            print(p)
            return("ok")
        }
       
    })
    return(PFS.data)    
})
dev.off()

############################################7.Obradovic_2021_Cell_scRNA_ccRCC###############################################
#### reference code:https://github.com/Aleksobrad/single-cell-rcc-pipeline/blob/master/bulkRNASeq_validation_cohort_analysis.R
## load Brian Rini dataset
rini_metadata <- read.csv("/data/active_data/lzl/RenalTumor-20200713/Data/Obradovic_2021_Cell_scRNA_ccRCC/validation_dataset_metadata.csv")
rini_data <- read.table("/data/active_data/lzl/RenalTumor-20200713/Data/Obradovic_2021_Cell_scRNA_ccRCC/validation_dataset_bulkRNASeq_counts.gct",sep="\t",header=T)
rini_data <- rini_data[which(rini_data[,2] %in% names(which(table(rini_data[,2])==1))),]
rownames(rini_data) <- rini_data[,2]
rini_data <- rini_data[,3:ncol(rini_data)]
colnames(rini_data) <- unlist(lapply(strsplit(colnames(rini_data),"\\."),function(x){x[3]}))
rini_metadata <- rini_metadata[which(rini_metadata$rna.seq != "#N/A" & rini_metadata$rna.seq != "0"),]
rini_metadata$rna.seq <- unlist(lapply(strsplit(as.character(rini_metadata$rna.seq),"-"),function(x){x[3]}))
rini_metadata <- rini_metadata[!grepl("DUPLICATE",rini_metadata$Batch.ID),]
rownames(rini_metadata) <- rini_metadata$rna.seq
rini_data <- rini_data[,rownames(rini_metadata)]

# dead due to cancer
rini_metadata$Cause.of.death <- as.character(rini_metadata$Cause.of.death)
idx <- which(rini_metadata$Cause.of.death %in% c("rcc", "RCC", "metastatic renal cell carcinoma", "Unknown, likely RCC"))
rini_metadata$dead.due.RCC <- rep("No", nrow(rini_metadata))
rini_metadata$dead.due.RCC[idx] <- "Yes"

#DESeq2 analysis
data <- rini_data
data <- data[which(rowSums(data)>0),]
group <- rini_metadata$dead.due.RCC
data_normalized <- log10(apply(data,2,function(x){x/sum(x)*1000000})+1)
data_normalized_scaled <- data_normalized[which(apply(data_normalized,1,sd)>0),]
data_normalized_scaled <- t(apply(data_normalized_scaled,1,function(x){(x-mean(x))/sd(x)}))
pca <- prcomp(t(data_normalized))
plot(pca$x[,1],pca$x[,2],pch=19,cex=1,col=as.numeric(group)+1,xlab="PC1",ylab="PC2")
text(pca$x[,1], pca$x[,2]-.5, colnames(data_normalized),cex=0.8)

rini_metadata <- rini_metadata[setdiff(rownames(rini_metadata),c("1","62","63","104","28")),]
rini_data <- rini_data[,rownames(rini_metadata)]
data <- rini_data
data <- data[which(rowSums(data)>0),]
group <- rini_metadata$dead.due.RCC
data_normalized <- log10(apply(data,2,function(x){x/sum(x)*1000000})+1)
data_normalized_scaled <- data_normalized[which(apply(data_normalized,1,sd)>0),]
data_normalized_scaled <- t(apply(data_normalized_scaled,1,function(x){(x-mean(x))/sd(x)}))
pca <- prcomp(t(data_normalized))
plot(pca$x[,1],pca$x[,2],pch=19,cex=1,col=as.numeric(group)+1,xlab="PC1",ylab="PC2")
text(pca$x[,1], pca$x[,2]-.5, colnames(data_normalized),cex=0.8)

library(DESeq2)
colData <- matrix(NA,nrow=ncol(data),ncol=2)
rownames(colData) <- colnames(data)
colnames(colData) <- c("condition","type")
colData[,"condition"]=as.character(rini_metadata$dead.due.RCC)
colData[,"type"] <- "paired-end"
colData <- as.data.frame(colData)
dds <- DESeqDataSetFromMatrix(countData = data,colData = colData,design= ~ condition)
dds <- DESeq(dds)
res <- results(dds, name="condition_Yes_vs_No")
resOrdered <- res[order(res$pvalue),]

library(atools)
set.seed(1234)
res <- res[which(!is.na(res$pvalue)),]
x <- res$log2FoldChange
#x=res$pvalue
names(x) <- rownames(res)
x <- x[which(!is.na(x))]
g <- rep(100, length(tumor.TFs$Name))
names(g) <- tumor.TFs$Name
set.seed(1234)
ledge <- gsea(signature = x, geneset = g, twoTails = F, pout = TRUE, per = 1000,colSig = c(.15,0.45, 0.3, 1), colHit = c(.05,0.58, 0.1, 2))$ledge

means <- apply(data_normalized, 1, mean)
sds <- apply(data_normalized, 1, sd)
data_normalized_scaled <- (data_normalized-means)/sds
data_normalized_scaled <- data_normalized_scaled[which(sds>0),]
set.seed(1234)
#gseaBySample <- apply(data_normalized_scaled,2,function(x){gsea(signature = x, geneset =g, twoTails = F, pout = T, per = 500)$nes})

clinical_dat <- rini_metadata
clinical_dat$TumorTFsGSEA <- gseaBySample
p=ggplot(clinical_dat, aes(x=dead.due.RCC, y=TumorTFsGSEA,fill=dead.due.RCC)) +
  geom_boxplot() + 
  #geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5,fill="black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Tumor TFs Enrichment: RCC vs No RCC")+ylab("NES")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=10, face="bold"),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),legend.position = "none")

#### overall suvival
rini_metadata$OS.event <- as.character(rini_metadata$Last.contact)
rini_metadata$OS.event <- gsub("Death notification", "1", rini_metadata$OS.event)
rini_metadata$OS.event <- gsub("Appointment", "0", rini_metadata$OS.event)
rini_metadata$OS.event <- as.numeric(rini_metadata$OS.event)
rini_metadata$Overall.Survival <- as.numeric(rini_metadata$Overall.Survival)

pdf("Obradovic_2021_Cell_bulk_ccRCC.four.TFs.survival.pdf")
survival.res <- lapply(four.TFs, function(x){
    # all cohort
    exp.values <- as.numeric(data_normalized_scaled[x,])
    med.value <- median(exp.values)
    group.label <- rep("Low group", length(exp.values))
    idx <- which(exp.values >= med.value)
    group.label[idx] <- "High group"
    if(length(idx) > length(exp.values)*0.9){
        cat("less sample across each group!\n")
    }else{
        OS.data <- data.frame(Patient_ID = rini_metadata$rna.seq, event = rini_metadata$OS.event, time = rini_metadata$Overall.Survival, sample.label = group.label)
        p <- plot.surv(OS.data, risk.table = T, HR = T, ylab = "Overall Survival", main = x, surv.median.line = "hv", xlab = "Time (Month)")
        print(p)
    }
})
dev.off()

pdf("Obradovic_2021_Cell_bulk_ccRCC.tumorTFs&targetGenes.survival.pdf")
res <- lapply(names(target.genes.list), function(x){
    idx <- match(target.genes.list[[x]], rownames(data_normalized_scaled))
    interest.matrix <- data_normalized_scaled[na.omit(idx),]
    signature.score <- colSums(interest.matrix)/nrow(interest.matrix)
    median.score <- median(signature.score)
    high.group <- which(signature.score >= median.score)
    sample.label <- rep("Low group", length(signature.score))
    sample.label[high.group] <- "High group"

    if(length(high.group) > length(signature.score)*0.9){
        cat("less sample across each group!\n")
    }else{
        OS.data <- data.frame(Patient_ID = rini_metadata$rna.seq, event = rini_metadata$OS.event, time = rini_metadata$Overall.Survival, sample.label = sample.label)
        p <- plot.surv(OS.data, risk.table = T, HR = T, ylab = "Overall Survival", main = x, surv.median.line = "hv", xlab = "Time (Month)")
        print(p)
    }
})
dev.off()

############################################8.CORCES_2018_Science_ATACseq_TCGA###############################################
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