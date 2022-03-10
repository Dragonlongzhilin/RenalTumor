#' @description: identify the tumor-specific TFs

library(Signac)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
set.seed(101)
library(ggpubr)
library(openxlsx)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(matrixStats)
library(future)
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 50000 * 1024^2) #50G

source(file = "/home/longzhilin/Analysis_Code/Visualization/colorPalettes.R")
setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC")

scATAC.data <- readRDS("scATAC.data.pro.rds")
motif.info <- data.frame(originName = names(scATAC.data@assays$Peaks@motifs@motif.names), TF = unlist(scATAC.data@assays$Peaks@motifs@motif.names))
rownames(motif.info) <- NULL
motif.info$originName <- gsub("_", "-", motif.info$originName)

#####基于chromVAR数据分析
cellType.motifs.chromVAR <- readRDS("5.Motif/cellType.motifs.chromVAR.human_pwms_v2.rds")

scRNA.data <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/data.merge.pro.rds")
###################1.plot motifs deviation score heatmap and number ratioplot
top5 <- sapply(cellType.motifs.chromVAR, function(x){
    x <- arrange(x, desc(avg_log2FC))
    return(x$gene[1:3])
})
top5 <- unique(as.character(top5))
top5 <- c(top5, "EOMES", "TBX10")
####Plot --- cellType specific motif heatmap
sig.motifs <- cellType.motifs.chromVAR %>% 
                                bind_rows %>%
                                filter(avg_log2FC > 1 & p_val_adj < 0.05) %>%
                                select(gene) %>%
                                distinct()
idx <- match(sig.motifs$gene, motif.info$TF)
sig.motifNames <- motif.info$originName[idx]

motifs.avgExp <- AverageExpression(scATAC.data, features = sig.motifNames, assays = "chromvar")
motifs.avgExp <- motifs.avgExp$chromvar
rownames(motifs.avgExp) <- sig.motifs$gene
zScore <- function(x){(x - mean(x)) /sd(x)}
motifs.avgExp.scale <- apply(motifs.avgExp, 1, zScore) %>% t() # row: TF; column: cell type
mark.idx <- match(top5, rownames(motifs.avgExp.scale))
pdf("5.Motif/Analysis/chromVAR.cellType.sig.motifs.heatmap.pdf")
ha <- rowAnnotation(link = anno_mark(at = mark.idx, labels = top5, link_width = unit(2, "mm"), labels_gp = gpar(fontsize = 5), padding = unit(1, "mm")))
Heatmap(motifs.avgExp.scale, name = "Deviation score", 
        width = unit(6, "cm"), height = unit(8, "cm"), right_annotation = ha,
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
        show_row_dend = F, show_column_dend = F, show_column_names = T, show_row_names = F, 
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 8), by_row = T))
Heatmap(motifs.avgExp.scale, name = "Deviation score", 
        width = unit(6, "cm"), height = unit(8, "cm"), left_annotation = ha, 
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
        show_row_dend = F, show_column_dend = F, show_column_names = T, show_row_names = F, 
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 8), by_row = T))      
dev.off()

####Plot --- top5 cell-type specific TFs
cellType.sig.motifs <- lapply(names(cellType.motifs.chromVAR), function(x){
    sig <- cellType.motifs.chromVAR[[x]] %>% filter(avg_log2FC > 1 & p_val_adj < 0.05) %>% mutate(celltype = x)
    return(sig)
})
names(cellType.sig.motifs) <- names(cellType.motifs.chromVAR)
sig.motifs.num <- cellType.sig.motifs %>% bind_rows()
sig.motifs.num <- as.data.frame(table(sig.motifs.num$celltype))
pdf("5.Motif/Analysis/chromVAR.cellType.sig.motifs.barplot.pdf")
ggbarplot(sig.motifs.num, x="Var1", y="Freq", fill = "Var1", color = "Var1",
          sort.by.groups=FALSE, sort.val = "desc", #不按组排序
          label = T, xlab = "", ylab = "Cell Number") + theme(legend.position="none") + rotate_x_text(30)
dev.off()

###################2.Screen cell-specific enriched motifs
source(file = "/home/longzhilin/Analysis_Code/code/analysis.diff.survival.TCGA.R")
DESeq2.normalized_counts <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/DESeq2.normalized_counts.rds")
DESeq2.normalized_counts <- log2(DESeq2.normalized_counts+1)
DESeq2.result <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/DESeq2.result.rds")
clin.data <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/clin.data.rds")

# extreact the avg_log2FC and FDR value
order.TFs <- cellType.motifs.chromVAR[[1]]$gene
motifs.FC <- sapply(cellType.motifs.chromVAR, function(x){
  index <- match(order.TFs, x$gene)
  return(round(x$avg_log2FC[index], 2))
})
rownames(motifs.FC) <- order.TFs

motifs.fdr <- sapply(cellType.motifs.chromVAR, function(x){
  index <- match(order.TFs, x$gene)
  return(x$p_val_adj[index])
})
rownames(motifs.fdr) <- order.TFs
a <- motifs.fdr
a[which(a==0)] <- 2
motifs.fdr[which(motifs.fdr==0)] <- min(a)*0.001 ###4.940656e-324, Multiply the minimum value by 0.01, instead of 0
motifs.fdr <- round(-log10(motifs.fdr), 2)

#3.deviation score
all.motifs.avgExp <- AverageExpression(scATAC.data, features = motif.info$originName, assays = "chromvar")
all.motifs.avgExp <- all.motifs.avgExp$chromvar
rownames(all.motifs.avgExp) <- motif.info$TF

################################## screening Strategy
#1.deviation score: sd > median(sd) + 4*mad
#2.Tumor cells: fdr> 0.0001 & avg_log2FC>=4
source(file = "/home/longzhilin/Analysis_Code/DataScience/MAD.R")
avg.sd <- rowSds(all.motifs.avgExp)
names(avg.sd) <- rownames(all.motifs.avgExp)
sd.mad <- DoubleMAD(avg.sd)
mad.threshold <- median(avg.sd) + 4*sd.mad[2] #Take the right

#Only in cell type: FDR < 0.0001 & log2FC>=4; log2FC<1 in other cell types
sig.label <- sapply(order.TFs, function(x){
     if(avg.sd[x] > mad.threshold){
        fdr <- motifs.fdr[x,]
        FC <- motifs.FC[x,]
        sig.fdr <- which(fdr > -log10(0.0001))
        sig.fc <- which(FC >= 4)  
        others <- which(FC >= 1)

        sig.idx <- intersect(sig.fdr, sig.fc)
        common.len <- length(sig.idx)

        if(common.len==1 & length(others)==1){
        #if(common.len==1){
            return("specific")
        }else{
            return("common")
        }
    }else{
        return("low variation")
    }
})
cellType.high.specific.TF <- names(which(sig.label=="specific")) # 49
heatmapEM.fdr <- motifs.fdr[cellType.high.specific.TF,]
heatmapEM.FC <- motifs.FC[cellType.high.specific.TF,]
cellType.sd <- avg.sd[cellType.high.specific.TF]
write.xlsx(list(FDR = heatmapEM.fdr, FC = heatmapEM.FC), file = "5.Motif/Analysis/cellType.specific.TFs.xlsx", sheetName = c("FDR", "Log2FC"), rowNames = T)

#### Tumor specific TFs
idx <- which(colnames(heatmapEM.FC) == "Tumor")
index <- which(heatmapEM.FC[, idx] >= 4)
tumor.specific.TFs <- data.frame(Name = rownames(heatmapEM.fdr)[index], FDR = heatmapEM.fdr[index, idx], avg_log2FC = heatmapEM.FC[index, idx], sd = cellType.sd[index])
heatmapEM.fdr.tumor <- heatmapEM.fdr[index,]
heatmapEM.FC.tumor <- heatmapEM.FC[index,]
saveRDS(tumor.specific.TFs, file = "5.Motif/Analysis/tumor.specific.TFs.rds")

## Filter the TF result
source(file = "/home/longzhilin/Analysis_Code/Visualization/Plot.EnhancedVolcano.R")
##plot --- heatmap
gene.labels <- c("HNF1A", "HNF1B", "HNF4G", "HNF4A", 
                 "HOXC5", "OTP", "ISL1", "VENTX", 
                 "HOXB5", "HOXB4", "HOXA4", "HOXB7", "HOXB3", "HOXB2",
                 "BARX2", "POU6F2", "JUNB", "JUN", "JUND", "BATF")
mark.idx <- match(gene.labels, rownames(heatmapEM.FC.tumor))
pdf("5.Motif/Analysis/tumor.specific.TFs.heatmap.pdf")
#Plot.EnhancedVolcano(TF.DESeq2, x = "log2FoldChange", y = "padj", select.num = nrow(TF.DESeq2), drawConnectors = F)
ha = rowAnnotation(FDR = anno_barplot(heatmapEM.fdr.tumor[,idx], border = F, gp = gpar(fill = "#e5e4e2"), direction = "reverse"))
col_fun <- colorRamp2(c(min(heatmapEM.FC.tumor), 0, max(heatmapEM.FC.tumor)), c("blue", "white", "red"))
ht2 <- Heatmap(heatmapEM.FC.tumor, left_annotation = ha, col = col_fun, 
               show_row_dend = F, show_column_dend = F, 
               column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 8), 
               name = "avg_log2FC", width = unit(7.5, "cm"), height = unit(10, "cm"), 
               show_heatmap_legend = T, 
               heatmap_legend_param = list(col_fun = col_fun, title = "avg_log2FC", direction = "horizontal", title_position = "topcenter")) +
rowAnnotation(link = anno_mark(at = mark.idx, labels = gene.labels, link_width = unit(3, "mm"), labels_gp = gpar(fontsize = 8), padding = unit(1, "mm")))
draw(ht2, heatmap_legend_side = "top")
dev.off()

# plot HNF1A, HNF1B, HNF4G
source("/home/longzhilin/Analysis_Code/SingleCell/FindRegion.R")
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

pdf("5.Motif/Analysis/HNF.coverage.pdf", width = unit(3, "inches"), height = unit(3, "inches"))
#HNF1A <- FindRegion(object = scATAC.data, region = "HNF1A", assay = "Peaks", extend.upstream = 0, extend.downstream = 0)
HNF1A <- StringToGRanges(c("chr12-120977001-120980000"))
regions <- StringToGRanges(c("chr12-120977001-120981801"))
p <- CoveragePlot(
  object = scATAC.data,
  region = HNF1A,
  links = F,
  peaks = F,
  extend.upstream = 0,
  extend.downstream = 0,
  region.highlight = regions)
print(p)

#HNF1B <- FindRegion(object = scATAC.data, region = "HNF1B", assay = "Peaks", extend.upstream = 0, extend.downstream = 0)
HNF1B <- StringToGRanges(c("chr17-37720000-37747000"))
regions <- StringToGRanges(c("chr17-37741201-37747000"))
p <- CoveragePlot(
  object = scATAC.data,
  region = HNF1B,
  links = F,
  peaks = F,
  extend.upstream = 0,
  extend.downstream = 0,
  region.highlight = regions)
print(p)

HNF4G <- FindRegion(object = scATAC.data, region = "HNF4G", assay = "Peaks", extend.upstream = 0, extend.downstream = 0)
HNF4G <- StringToGRanges(c("chr8-75403400-75425000"))
regions <- StringToGRanges(c("chr8-75403401-75410000"))
p <- CoveragePlot(
  object = scATAC.data,
  region = HNF4G,
  links = F,
  peaks = F,
  extend.upstream = 0,
  extend.downstream = 0,
  region.highlight = regions)
print(p)

dev.off()


# survival in TCGA-KIRC data
write.xlsx(tumor.specific.TFs, file = "5.Motif/Analysis/tumor.specific.TFs.xlsx", sheetName = c("Tumor specific TFs"), rowNames = T)
pdf("5.Motif/Analysis/tumor.specific.TFs.survival.pdf")
tumor.TFs.signature.res <- analysis.diff.survival.TCGA(interest.gene = tumor.specific.TFs$Name, diff.gene.pro = DESeq2.result, exp.data.process = DESeq2.normalized_counts, clin.data = clin.data, EnhancedVolcano.plot = F, Box.plot = F, main = "tumor.specific.TFs", meta.signature = T, single.signature = T)
dev.off()

# surcivial in ICB data
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
single <- as.list(tumor.specific.TFs$Name)
names(single) <- tumor.specific.TFs$Name
TF.list <- c(list(All = tumor.specific.TFs$Name), single)
pdf("5.Motif/Analysis/tumor.specific.TFs.ICB.pdf")
cox.res <- RCC.icb.analysis(signature.list = TF.list , expresionMatrix = normalized_expression, clincal.info = patient.info.RNA)
dev.off()


####footprint plot
DefaultAssay(scATAC.data) <- "Peaks"
scATAC.data <- Footprint(
  object = scATAC.data,
  motif.name = c("OTP", "ISL1", "VENTX", "HOXC5", "HNF1A", "HNF1B", "HNF4G"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)
pdf("5.Motif/Analysis/Footprint.pdf")
PlotFootprint(scATAC.data, features = "OTP")
PlotFootprint(scATAC.data, features = "ISL1")
PlotFootprint(scATAC.data, features = "VENTX")
PlotFootprint(scATAC.data, features = "HOXC5")
PlotFootprint(scATAC.data, features = "HNF1A")
PlotFootprint(scATAC.data, features = "HNF1B")
PlotFootprint(scATAC.data, features = "HNF4G")
dev.off()

library(ggplot2)
library(ggseqlogo)
PWMatrixToProbMatrix <- function(x){
        if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
        m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
        m <- t(t(m)/colSums(m))
        m
}
library(chromVARmotifs)
data("human_pwms_v2")
index <- grep("OTP|ISL1|VENTX|HOXC5", names(human_pwms_v2))
PPM.list <- lapply(index, function(y){
    PPM <- PWMatrixToProbMatrix(human_pwms_v2[[y]])
})
pdf("5.Motif/Analysis/tumor.TF.logo.pdf")
names(PPM.list) <- as.character(unlist(scATAC.data@assays$Peaks@motifs@motif.names[index]))
p <- ggseqlogo(PPM.list) + theme(plot.title = element_text(hjust = 0.5))
print(p)
dev.off()

Idents(scRNA.data) <- scRNA.data$cellType_low
tumor.TFs.avgExp <- AverageExpression(scRNA.data, features = tumor.specific.TFs$Name, assays = "RNA", slot = "data")
tumor.TFs.avgExp <- scale(t(tumor.TFs.avgExp$RNA))

pdf("5.Motif/Analysis/tumor.TFs.expression.pdf")
Heatmap(tumor.TFs.avgExp, 
        show_row_dend = F, show_column_dend = F, 
        column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 8), 
        name = "Expression", height = unit(6, "cm"), 
        show_heatmap_legend = T, 
        heatmap_legend_param = list(title = "Expression", direction = "horizontal", title_position = "topcenter"))
dev.off()