#' @description: cellphoneDB analysis

library(psych)
library(qgraph)
library(igraph)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(dplyr)
library(vegan)
library(Seurat)

source(file = "/home/longzhilin/Analysis_Code/code/cellphoneDB.dotplot.R")
source(file = "/home/longzhilin/Analysis_Code/code/Plot.CircularChart.R")

# TCGA data
source(file = "/home/longzhilin/Analysis_Code/code/analysis.diff.survival.TCGA.R")
source(file = "/home/longzhilin/Analysis_Code/code/survival.combined.R")
DESeq2.normalized_counts <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/DESeq2.normalized_counts.rds")
DESeq2.normalized_counts <- log2(DESeq2.normalized_counts+1)
DESeq2.result <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/DESeq2.result.rds")
clin.data <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/clin.data.rds")

# ICB data
source(file = "/home/longzhilin/Analysis_Code/SurvivalAnalysis/Cox.function.R")
source(file = "/home/longzhilin/Analysis_Code/code/RCC.ICB.analysis.R")
normalized_expression <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/ICB.therapy/normalized_expression.log2.rds")
patient.info.RNA <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/ICB.therapy/patient.info.RNA.rds")

patient.info.RNA$Sex <- gsub("^F|FEMALE|Female$", 0, patient.info.RNA$Sex)
patient.info.RNA$Sex <- gsub("^M|Male|MALE$", 1, patient.info.RNA$Sex)
patient.info.RNA$Tumor_Sample_Primary_or_Metastasis <- gsub("PRIMARY", 0, patient.info.RNA$Tumor_Sample_Primary_or_Metastasis)
patient.info.RNA$Tumor_Sample_Primary_or_Metastasis <- gsub("METASTASIS", 1, patient.info.RNA$Tumor_Sample_Primary_or_Metastasis)

data.merge <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/data.merge.pro.rds")
DefaultAssay(data.merge) <- "RNA"
Idents(data.merge) <- data.merge$cellType_low
Macrophage <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/5.Immune/Macrophage/sub.scRNA.harmony.pro.rds") # 3 cluster
CD8.T <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/5.Immune/CD8T/sub.scRNA.harmony.pro.rds")
idx1 <- match(colnames(Macrophage), colnames(data.merge)) # 4319 cells
idx2 <- match(colnames(CD8.T), colnames(data.merge)) # 3728 cells
cell.Type <- data.merge$cellType_low
cell.Type[idx1] <- as.character(Macrophage$cellType3)
cell.Type[idx2] <- paste0("CD8+ T-", CD8.T$cellType3)
data.merge$cellType2 <- cell.Type
cell.Type.new <- setdiff(unique(cell.Type), c("CD8+ T cell", "Macrophage"))
data.merge.new <- subset(data.merge, subset = cellType2 %in% cell.Type.new)
DefaultAssay(data.merge.new) <- "RNA"

setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/6.CrossTalk/CellPhoneDB")
mean.threshold <- 1
pval.threshold <- 0.05
colors <- c("#00A087", "#4DBBD5", "#E64B35", "#3C5488")

## Filter custom ligand receptor, partner_a is ligand, partner_b is receptor
# 3511*7
ligand.receptor.data <- read.table("/data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/ligand_receptor/ligand_receptor.data.csv", header = T, sep = ",", stringsAsFactors = F)
ligand_receptor <- paste0(ligand.receptor.data$partner_a, "_", ligand.receptor.data$partner_b)

All.interaction.number <- read.table("interaction_count.txt", header = T, sep = "\t", stringsAsFactors = F)
cell.types <- All.interaction.number$X
order.value <- c("CD4+ T cell", "Treg", paste0("CD8+ T-", levels(CD8.T$cellType3)), "Proliferative T cell", "NK/NKT cell", "B cell",
                 levels(Macrophage$cellType3), "Monocyte", "Dendritic cell", "Neutrophil", "Mast cell",
                 "Endothelium VCAM1+", "Endothelium VCAM1-", "Mesangial cell", "Tumor")
data.merge.new$cellType2 <- factor(data.merge.new$cellType2, levels = order.value)
Idents(data.merge.new) <- data.merge.new$cellType2
order.value <- gsub("CD8\\+ T-Exhausted IEG", "CD8\\+ T-IEG", order.value)
####################################Overall Cell Interaction Mode######## ##########################
#Only keep user_curated, remove the interaction pair in cellphoneDB (including complex)
mypvals <- read.delim("Renal/pvalues.txt", check.names = FALSE) # 3479 *411
mymeans <- read.delim("Renal/means.txt", check.names = FALSE)
sig.means <- read.delim("Renal/significant_means.txt", check.names = FALSE)

mypvals.usr <- mypvals[which(mypvals$annotation_strategy == "user_curated"),] # 2930 * 411
mymeans.usr <- mymeans[which(mymeans$annotation_strategy == "user_curated"),]
sig.means.usr <- sig.means[which(sig.means$annotation_strategy == "user_curated"),]

# Identification ligand and receptor, used to identify the first gene
mypvals.usr$ligand <- rep("Ligand", nrow(mypvals.usr))
mymeans.usr$ligand <- rep("Ligand", nrow(mymeans.usr))
sig.means.usr$ligand <- rep("Ligand", nrow(sig.means.usr))

res <- for(i in 1:nrow(mypvals.usr)){
    index1 <- match(mypvals.usr$interacting_pair[i], ligand_receptor)
    if(is.na(index1)){
        mypvals.usr$ligand[i] <- "Receptor"
        mymeans.usr$ligand[i] <- "Receptor"
    }

    index2 <- match(sig.means.usr$interacting_pair[i], ligand_receptor)
    if(is.na(index2)){
        sig.means.usr$ligand[i] <- "Receptor"
    }
}
# Move the position of the ligand column, position the 11th column
mypvals.usr <- mypvals.usr[,c(1:10, 412, 11:411)]
mymeans.usr <- mymeans.usr[,c(1:10, 412, 11:411)]
sig.means.usr <- sig.means.usr[,c(1:10, 413, 11:412)]

# 2930 ligand-receptors
saveRDS(mypvals.usr, file = "mypvals.usr.rds")
saveRDS(mymeans.usr, file = "mymeans.usr.rds")
saveRDS(sig.means.usr, file = "sig.means.usr.rds")
write.table(mypvals.usr, file = "Renal/mypvals.usr.txt", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(mymeans.usr, file = "Renal/mymeans.usr.txt", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(sig.means.usr, file = "Renal/sig.means.usr.txt", row.names = F, col.names = T, sep = "\t", quote = F)

####Statistics of the ligand-receptor interactions of each cell type
#ligand
sig.interaction.cellType <- lapply(cell.types, function(x){
    mymeans.usr %>% dplyr::select("interacting_pair", "ligand", starts_with(x), ends_with(x)) %>% reshape2::melt() -> meansdf
    colnames(meansdf)<- c("interacting_pair", "ligand", "CC","means")
    mypvals.usr %>% dplyr::select("interacting_pair", "ligand", starts_with(x), ends_with(x)) %>% reshape2::melt()-> pvalsdf
    colnames(pvalsdf)<- c("interacting_pair", "ligand", "CC","pvals")
    pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair, "_", pvalsdf$CC)
    meansdf$joinlab<- paste0(meansdf$interacting_pair, "_", meansdf$CC)
    pldf <- merge(pvalsdf, meansdf, by = "joinlab")
    sig.interaction <- pldf %>% filter(means >= mean.threshold & pvals < pval.threshold)
    sig.interaction <- sig.interaction[, -c(6,7,8)]
    colnames(sig.interaction) <- c("joinlab", "interacting_pair", "ligand", "CC", "pvals", "means")

    ####Split the interaction between ligand and receptor
    a <- apply(sig.interaction, 1, function(y){
        CC.interaction <- unlist(strsplit(as.character(y[4]), "\\|"))
        if((y[3] == "Ligand" & CC.interaction[1] == x) | (y[3] == "Receptor" & CC.interaction[2] == x)){
            return("ok")
        }else{
            return("no")
        }
    })
    ligand.interactions <- sig.interaction[which(a=="ok"),]
    receptor.interactions <- sig.interaction[which(a=="no"),]

    sig.interaction <- sig.interaction[order(sig.interaction$means, decreasing = T),]
    ligand.interactions <- ligand.interactions[order(ligand.interactions$means, decreasing = T),]
    receptor.interactions <- receptor.interactions[order(receptor.interactions$means, decreasing = T),]      

    return(list(All = sig.interaction, Ligand = ligand.interactions, Receptor = receptor.interactions))
})
names(sig.interaction.cellType) <- cell.types
saveRDS(sig.interaction.cellType, file = "sig.interaction.cellType.rds")

a <- lapply(sig.interaction.cellType, function(x){
    return(x$All)
})
library(openxlsx)
write.xlsx(a, file = "sig.interaction.cellType.xlsx", sheetName = names(sig.interaction.cellType), rowNames = F)

####----Drawing a heat map of interaction

interaction.counts <- sapply(names(sig.interaction.cellType), function(x){
    interaction.pairs <- sig.interaction.cellType[[x]]$All
    cellType.counts <- sapply(names(sig.interaction.cellType), function(y){
        a <- t(as.data.frame(strsplit(x = as.character(interaction.pairs$CC), split = "\\|")))
        idx1 <- which(a[,1]==y) 
        idx2 <- which(a[,2]==y)
        return(length(unique(c(idx1, idx2))))
    })
    idx <- which(cellType.counts == max(cellType.counts))

    self.pairs <- length(which(interaction.pairs$CC == paste0(x, "|", x)))
    cellType.counts[idx] <- self.pairs
    return(cellType.counts)
})
saveRDS(interaction.counts, file = "sig.interaction.counts.rds")
pdf("sig.interaction.number.heatmap.pdf")
col_fun <- c("#124F8B", "#FED9B8", "#8B0B50")
Heatmap(as.matrix(interaction.counts), name = "number of interactions", col = col_fun, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), width = unit(8, "cm"), height = unit(8, "cm"))
Heatmap(log2(as.matrix(interaction.counts)), name = "number of interactions", col = col_fun, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), width = unit(8, "cm"), height = unit(8, "cm"))
#只提取immune cell
idx <- which(rownames(interaction.counts) %in% c("Endothelium VCAM1-", "Endothelium VCAM1+", "Mesangial cell"))
Heatmap(as.matrix(interaction.counts[-idx, -idx]), name = "number of interactions", col = col_fun, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), width = unit(8, "cm"), height = unit(8, "cm"))
Heatmap(log2(as.matrix(interaction.counts[-idx, -idx])), name = "number of interactions", col = col_fun, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), width = unit(8, "cm"), height = unit(8, "cm"))
dev.off()

####-----count the interaction number
interaction.type.number <- sapply(sig.interaction.cellType, function(x){
    return(c(nrow(x$Ligand), nrow(x$Receptor)))
})
rownames(interaction.type.number) <- c("Ligand", "Receptor")

##1.all cell types
interaction.type.number <- data.frame(interaction.type.number, check.names = F)
interaction.type.number$Type <- c("Ligand", "Receptor")
library(tidyr)
res <- gather(interaction.type.number, cellType, count, -Type)
saveRDS(res, file = "Ligand.Receptor.interaction.rds")

res$cellType <- factor(res$cellType, levels = order.value)
res.ligand <- res[which(res$Type == "Ligand"),]
res.receptor <- res[which(res$Type == "Receptor"),]
type <- rep("Lymphoid", nrow(res.ligand))
type[c(1, 8, 13:16, 19)] <- "Myeloid"
type[10] <- "Tumor"
type[c(4, 6, 12)] <- "Other"
res.ligand$type <- factor(type, levels = c("Lymphoid", "Myeloid", "Tumor", "Other"))
res.receptor$type <- factor(type, levels = c("Lymphoid", "Myeloid", "Tumor", "Other"))
pdf("Ligand.Receptor.interaction.pdf")
ggbarplot(res, x="cellType", y="count", color = "Type", xlab = "", ylab = "Number of interactions", fill = "Type", palette = c("#E18727", "#0072B5"), x.text.angle = 60, position = position_dodge())
ggbarplot(res.ligand, x="cellType", y="count", color = "type", xlab = "", ylab = "Number of interactions", fill = "type", palette = colors, sort.val = "desc", sort.by.groups = FALSE, x.text.angle = 60)
ggbarplot(res.receptor, x="cellType", y="count", color = "type", xlab = "", ylab = "Number of interactions", fill = "type", palette = colors, sort.val = "desc", sort.by.groups = FALSE, x.text.angle = 60)

#2.focus on the immune cell
other <- c("Endothelium VCAM1-", "Endothelium VCAM1+", "Mesangial cell", "Tumor")
index1 <- which(res$cellType %in% other)

res.immune <- res[-index1,]
ggbarplot(res.immune, x="cellType", y="count", color = "Type", xlab = "", ylab = "Number of interactions", fill = "Type", palette = c("#E18727", "#0072B5"), x.text.angle = 60, position = position_dodge())
index1 <- which(res.ligand$cellType %in% other)
ggbarplot(res.ligand[-index1,], x="cellType", y="count", color = "type", xlab = "", ylab = "Number of interactions", fill = "type", palette = colors, sort.val = "desc", sort.by.groups = FALSE, x.text.angle = 60) 
ggbarplot(res.receptor[-index1,], x="cellType", y="count", color = "type", xlab = "", ylab = "Number of interactions", fill = "type", palette = colors, sort.val = "desc", sort.by.groups = FALSE, x.text.angle = 60)

# heatmap展示ligand和receptor情况
ligand.receptor.data <- data.frame(Ligand = res.ligand[,3], Receptor = res.receptor[,3])
rownames(ligand.receptor.data) <- res.ligand$cellType
row_split <- res.ligand$type
Heatmap(ligand.receptor.data, name = "Number of Interaction", col = colorRamp2(c(0,1200), c("white", "red")), show_row_dend = F, show_column_dend = F, row_split = row_split, 
        width = unit(2, "cm"), height = unit(10, "cm"), cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.0f", ligand.receptor.data[i, j]), x, y, gp = gpar(fontsize = 7))}
        )
dev.off()

# heatmap shows the conditions of ligand and receptor
All.interaction.number <- sapply(sig.interaction.cellType, function(x){
    return(nrow(x$All))
})
All.interaction.number <- data.frame(cellType = names(All.interaction.number), interaction.num = All.interaction.number)

class.type <- rep("Lymphoid", nrow(All.interaction.number))
class.type[c(1, 8, 13:16, 19)] <- "Myeloid"
class.type[which(All.interaction.number$cellType %in% c("Endothelium VCAM1+", "Endothelium VCAM1-", "Mesangial cell"))] <- "Other"
class.type[which(All.interaction.number$cellType %in% c("Tumor"))] <- "Tumor"
All.interaction.number$class <- factor(class.type, levels = c("Lymphoid", "Myeloid", "Tumor", "Other"))
pdf("sig.interaction.number.pdf")
ggbarplot(All.interaction.number, x="cellType", y="interaction.num", xlab = "", ylab = "Number of interactions", fill = "class", color = "white", palette = colors,sort.val = "desc",sort.by.groups = FALSE,x.text.angle = 60)
dev.off()
saveRDS(All.interaction.number, "All.interaction.number.rds")

#### Prepare circular plots for all cell types
track.matrix <- matrix(c(rep(0, nrow(All.interaction.number)), All.interaction.number$interaction.num), ncol = 2)
rownames(track.matrix) <- All.interaction.number[,1]
track.matrix <- track.matrix[order.value,]

#function --- plot circular
plot.circular <- function(track.matrix, gene.pair, mypvals, mymeans, data.merge = data.merge, idx = c(1, 3:10, 12), means.threshold = 1, pval.threshold = 0.01, risk.table = T, High.low = F, colors = c("#D95319", "#F39B7F", "#3B6793","#4285F4")){
    gene.a <- gene.pair[1]
    gene.b <- gene.pair[2]
    mymeans %>% dplyr::filter(gene_a == gene.a & gene_b == gene.b) %>% reshape2::melt() -> meansdf
    meansdf <- meansdf[, -idx]
    colnames(meansdf)<- c("interacting_pair", "ligand", "CC","means")
    mypvals %>% dplyr::filter(gene_a == gene.a & gene_b == gene.b) %>% reshape2::melt() -> pvalsdf
    pvalsdf <- pvalsdf[, -idx]
    colnames(pvalsdf)<- c("interacting_pair", "ligand", "CC","means")
    pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair, "_", pvalsdf$CC)
    meansdf$joinlab<- paste0(meansdf$interacting_pair, "_", meansdf$CC)
    pldf <- merge(pvalsdf, meansdf, by = "joinlab")
    pldf <- pldf[,-c(1,6:8)]
    colnames(pldf) <- c("interacting_pair", "ligand", "CC", "pvals", "means")
    sig.interaction <- pldf%>% filter(means >= means.threshold & pvals < pval.threshold)

    interaction.matrix <- apply(sig.interaction, 1, function(x){
        if(x[2] == "Receptor"){
            a <- unlist(strsplit(x[3], "\\|"))
            b <- c(a[2], a[1], x[5])
        }else{
            a <- unlist(strsplit(x[3], "\\|"))
            b <- c(a[1], a[2], x[5])
        }
        return(b)
    })
    interaction.matrix <- as.data.frame(t(interaction.matrix))
    interaction.matrix$means <- as.numeric(interaction.matrix$means)
    p <- Plot.CircularChart(track.matrix = track.matrix, interaction.matrix = interaction.matrix, 
                            bg.col = NULL, cex = 0.8,
                            link.col = NULL, arr.lwd = 0.3)
    print(p)

    survival.res <- analysis.diff.survival.TCGA(interest.gene = c(gene.a, gene.b), diff.gene.pro = DESeq2.result, exp.data.process = DESeq2.normalized_counts, clin.data = clin.data, EnhancedVolcano.plot = F, Box.plot = T, main = paste0(gene.a, "-", gene.b), meta.signature = F)
    combined.res <- survival.combined(geneA = gene.a, geneB = gene.b, clin.data = clin.data, DESeq2.result = DESeq2.result, DESeq2.normalized_counts = DESeq2.normalized_counts, risk.table = risk.table, High.low = High.low, colors = colors)

    p <- FeaturePlot(data.merge, features = gene.a, cols = c("lightgrey", "red"))
    print(p)
    p <- FeaturePlot(data.merge, features = gene.b, cols = c("lightgrey", "red"))
    print(p)
    p <- VlnPlot(data.merge, features = gene.a, pt.size = 0.5) + NoLegend()
    print(p)
    p <- VlnPlot(data.merge, features = gene.b, pt.size = 0.5) + NoLegend()
    print(p)
    return(interaction.matrix)
}

#########################################################analysis tumor cell################################################
tumor <- sig.interaction.cellType[["Tumor"]]
#Statistics and tumor interaction cell type distribution
#ligand
ligand.interact.tumor <- apply(tumor$Ligand, 1, function(x){
    a <- unlist(strsplit(x[4], "\\|"))
    return(a)
})
ligand.interact.tumor.freq <- data.frame(table(as.character(ligand.interact.tumor)))

idx <- which(ligand.interact.tumor.freq$Var1 == "Tumor")
ligand.interact.tumor.freq[idx, 2] <- ligand.interact.tumor.freq[idx, 2]-nrow(tumor$Ligand)

receptor.interact.tumor <- apply(tumor$Receptor, 1, function(x){
    a <- unlist(strsplit(x[4], "\\|"))
    return(a)
})
receptor.interact.tumor.freq <- data.frame(table(as.character(receptor.interact.tumor)))
idx <- which(receptor.interact.tumor.freq$Var1 == "Tumor")
receptor.interact.tumor.freq[idx, 2] <- receptor.interact.tumor.freq[idx, 2]-nrow(tumor$Receptor)

type <- rep("Lymphoid", nrow(receptor.interact.tumor.freq))
type[c(7,10,12:13,16:18)] <- "Myeloid"
type[20] <- "Tumor"
type[c(8:9, 11)] <- "Other"
receptor.interact.tumor.freq$Type <- factor(type, levels = c("Lymphoid", "Myeloid", "Tumor", "Other"))
ligand.interact.tumor.freq$Type <- factor(type, levels = c("Lymphoid", "Myeloid", "Tumor", "Other"))

ligand.interact.tumor.freq$Var1 <- factor(ligand.interact.tumor.freq$Var1, levels = order.value)
receptor.interact.tumor.freq$Var1 <- factor(receptor.interact.tumor.freq$Var1, levels = order.value)
pdf("Tumor/tumor.ligand.receptor.pdf")
ggbarplot(ligand.interact.tumor.freq, x="Var1", y="Freq", color = "Type", xlab = "", palette = colors, ylab = "Number of Ligand", fill = "Type", sort.val = "desc", sort.by.groups = FALSE, x.text.angle = 60) 
ggbarplot(receptor.interact.tumor.freq, x="Var1", y="Freq", color = "Type", xlab = "", palette = colors, ylab = "Number of Receptor", fill = "Type", sort.val = "desc", sort.by.groups = FALSE, x.text.angle = 60)
ggbarplot(ligand.interact.tumor.freq[-c(8,9,11),], x="Var1", y="Freq", color = "Type", xlab = "", palette = colors[1:3], ylab = "Number of Ligand", fill = "Type", sort.val = "desc", sort.by.groups = FALSE, x.text.angle = 60) 
ggbarplot(receptor.interact.tumor.freq[-c(8,9,11),], x="Var1", y="Freq", color = "Type", xlab = "", palette = colors[1:3], ylab = "Number of Receptor", fill = "Type", sort.val = "desc", sort.by.groups = FALSE, x.text.angle = 60)

ligand.receptor.data <- data.frame(Ligand = ligand.interact.tumor.freq[,2], Receptor = receptor.interact.tumor.freq[,2])
rownames(ligand.receptor.data) <- ligand.interact.tumor.freq$Var1
row_split <- factor(ligand.interact.tumor.freq$Type, levels = c("Lymphoid", "Myeloid", "Tumor", "Other"))
Heatmap(ligand.receptor.data, name = "Number of Interaction", col = colorRamp2(c(0,50), c("white", "red")), show_row_dend = F, show_column_dend = F, row_split = row_split, 
        width = unit(2, "cm"), height = unit(8, "cm"), cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.0f", ligand.receptor.data[i, j]), x, y, gp = gpar(fontsize = 7))}
        )

# total number
total.interaction <- data.frame(cellType = rownames(ligand.receptor.data), number = ligand.receptor.data[,1]+ligand.receptor.data[,2], type = ligand.interact.tumor.freq$Type)
total.interaction$type <- factor(total.interaction$type, levels = c("Lymphoid", "Myeloid", "Tumor", "Other"))
ggbarplot(total.interaction, x="cellType", y="number", color = "type", xlab = "", palette = colors, ylab = "Number of interaction", fill = "type", sort.val = "desc", sort.by.groups = FALSE, x.text.angle = 60) 
dev.off()
saveRDS(total.interaction, file = "Tumor/tumor.interaction.number.rds")

#### Interaction between cancer cells and other cells
#Screening significant interaction pairs-ligand
selected_rows <- unique(tumor$All$interacting_pair[which(tumor$All$means>=mean.threshold & tumor$All$pvals<pval.threshold)])

idx <- which(selected_rows %in% c("VEGFA_FLT1", "NRP1_VEGFA", "SPP1_CD44", "C5AR1_RPS19", "HLA-E_KLRD1", "EGFR_VEGFA", "VEGFA_GPC1",
                                  "CXCL8_SDC2", "APOE_TREM2", "APOE_VLDLR", "EGFR_VEGFA",
                                  "IL1B_SIGIRR", "HMGB1_CD163", "HLA-E_KLRK1", "CXCL8_ACKR1", "CXCL2_DPP4", "CD24_SIGLEC10", "CD14-RIPK1",
                                  "APOC1_VLDLR", "CXCL8_KDR", "CXCL8_CD79A", "CD27_CD70", "CD2_CD48", "HMMR_CALM1", "APP_NOTCH2",
                                  "CD44_TIMP3", "TLR4_HMGB1", "CD74_MIF", "MIF_CXCR4", "VEGFA_KDR", "VEGFA_EPHB2"))
selected_rows <- selected_rows[idx]

selected_columns <- sort(unique(as.character(tumor$All$CC)))
pdf("Tumor/tumor.sig.interaction.pdf")
p1 <- dot_plot(selected_rows = selected_rows, selected_columns = selected_columns, 
               breaks = c(-2, -1, 0, 1, 2), limit = c(-2, 2), legend.key.size = unit(0.1, "inches"),
               size = 8, max.dot.size = 3, means_path = "Renal/mymeans.usr.txt", pvalues_path = "Renal/mypvals.usr.txt")
print(p1)
dev.off()

pdf("Tumor/VEGFA.expression.pdf", width = 3, height = 3)
p <- DotPlot(data.merge.new, features = c("VEGFA", "NRP1", "KDR", "GPC1"), cols = c("#1e90ff", "#F15F30"), dot.scale = 3.5) + xlab("") + ylab("") + theme(axis.text.x = element_text(size = 8, angle = 90)) + theme(axis.text.y = element_text(size = 8)) + NoLegend()
print(p)
dev.off()
pdf("Tumor/VEGFA.expression-test.pdf", width = 3, height = 3)
p <- DotPlot(data.merge.new, features = c("VEGFA", "NRP1", "KDR", "GPC1"), cols = c("#1e90ff", "#F15F30"), dot.scale = 3.5) + xlab("") + ylab("") + theme(axis.text.x = element_text(size = 8, angle = 90)) + theme(axis.text.y = element_text(size = 8))
print(p)
dev.off()
pdf("Tumor/interested.interaction.pairs.pdf")
p <- VlnPlot(data.merge.new, features = c("MIF", "CD74", "CXCR4"), pt.size = 0, ncol = 2) + NoLegend()
print(p)
p <- VlnPlot(data.merge.new, features = c("VEGFA", "NRP1", "KDR", "GPC1"), pt.size = 0, ncol = 2) + NoLegend()
print(p)
p <- FeaturePlot(data.merge, features = c("VEGFA", "NRP1", "KDR", "GPC1"), cols = c("lightgrey", "red"))
print(p)
p <- VlnPlot(data.merge.new, features = c("CXCL2", "DPP4", "CD24", "SIGLEC10"), pt.size = 0, ncol = 2)
print(p)
p <- VlnPlot(data.merge.new, features = c("APOE", "APOC1", "TREM2", "VLDLR"), pt.size = 0, ncol = 2)
print(p)
p <- VlnPlot(data.merge.new, features = c("CD70", "CD27", "CALM1", "HMMR"), pt.size = 0, ncol = 2)
print(p)
dev.off()

features <- c("MIF",  "VEGFA", "GPC1", "CD24",
              "DPP4", "CD70", "CALM1", "APOE", "VLDLR") 
pdf("Tumor/tumor.sig.ligand.receptor.survival.pdf")
res <- analysis.diff.survival.TCGA(interest.gene = features, diff.gene.pro = DESeq2.result, exp.data.process = DESeq2.normalized_counts, clin.data = clin.data, EnhancedVolcano.plot = F, Box.plot = T, main = "cancer.sig.ligand.receptor", meta.signature = F)
features <- c("MIF",  "VEGFA", "GPC1", 
              "DPP4", "CD70", "CALM1", "APOE", "VLDLR")
features.list <- as.list(features)
names(features.list) <- features
cox.res <- RCC.icb.analysis(signature.list = features.list, expresionMatrix = normalized_expression, clincal.info = patient.info.RNA)
dev.off()


#### Prepare circular plots 
track.matrix <- matrix(c(rep(0, nrow(total.interaction)), total.interaction$number), ncol = 2)
rownames(track.matrix) <- total.interaction[,1]
track.matrix <- track.matrix[order.value,]

pdf("Tumor/tumor.interaction.circular.pdf")
plot.circular(track.matrix = track.matrix, gene.pair = c("CD24", "SIGLEC10"), mypvals = mypvals.usr, mymeans = mymeans.usr)
plot.circular(track.matrix = track.matrix, gene.pair = c("APOE", "TREM2"), mypvals = mypvals.usr, mymeans = mymeans.usr)
plot.circular(track.matrix = track.matrix, gene.pair = c("APOC1", "VLDLR"), mypvals = mypvals.usr, mymeans = mymeans.usr)
plot.circular(track.matrix = track.matrix, gene.pair = c("CXCL2", "DPP4"), mypvals = mypvals.usr, mymeans = mymeans.usr)
dev.off()


###########################Tumor VS TAM
tumor_TAM <- tumor$All[grep("TAM", tumor$All$CC),]
selected_rows <- unique(tumor_TAM$interacting_pair[which(tumor_TAM$means>=mean.threshold & tumor_TAM$pvals<pval.threshold)])

# features <- as.character(unlist(strsplit(selected_rows, "_")))
# pdf("Tumor/Myeloid/tumor.sig.ligand.receptor.survival.pdf")
# res <- analysis.diff.survival.TCGA(interest.gene = features, diff.gene.pro = DESeq2.result, exp.data.process = DESeq2.normalized_counts, clin.data = clin.data, EnhancedVolcano.plot = F, Box.plot = T, main = "Tumor.TAM", meta.signature = F)
# dev.off()
# write.csv(res,file = "Tumor/Myeloid/tumor.sig.ligand.receptor.survival.csv", row.names = F)

TAM.DEGs <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/5.Immune/Macrophage/cellType.DE.pro.rds")
Tumor.DEGs <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/2.Cluster/AnnotateCellType/cellType.sig.pos.DEGs.rds")
Tumor.DEGs <- Tumor.DEGs %>% filter(cluster == "Tumor") %>% filter(p_val_adj<0.05 & avg_log2FC>0.25)
TAM.DEGs <- TAM.DEGs %>% filter(p_val_adj<0.05 & avg_log2FC>0.25)
DEGs <- unique(c(Tumor.DEGs$gene, TAM.DEGs$gene))
interaction.genes <- unlist(str_split(selected_rows, pattern="_"))
sig.DEGs <- intersect(interaction.genes, DEGs)

interaction.info <- sapply(selected_rows, function(x){
    a <- unlist(str_split(x, "_"))
    first <- match(a[1], sig.DEGs)
    second <- match(a[2], sig.DEGs)
    if(!is.na(first)&!is.na(second)){
        return("ok")
    }else{
        return("no")
    }
})

idx <- which(interaction.info == "ok")
selected_rows <- selected_rows[idx]
idx <- which(selected_rows %in% c("ITGB2_VCAM1", "ITGB1_TIMP2", "CD99_CD81","TIMP1_CD63", "ITGB1_SPP1", "ITGA5_SPP1", "ITGB1_VEGFA", "CP_SLC40A1", "EGFR_ANXA1",
                                  "EGFR_MIF", "EGFR_S100A4", "EGFR_GSTP1", "IL1B_SIGIRR", "C3_CD81", "ITGB1_CD14", "ITGB2_CD14", "APOE_SDC2", "APOE_SORL1",
                                  "TNFRSF1B_GRN", "VEGFA_ITGA9", "SPP1_ITGA4", "SPP1_PTGER4", "S100A8_CD68", "S100A9_CD68", "HSPA8_LRP2",
                                  "TLR4_HSPA1A", "PKM_CD44", "SPP1_ITGA9", "ITGB1_LGALS1", "EGFR_GRN", "FN1_PLAUR", "CXCL14_CXCR4", "CXCL8_SDC2"))
selected_rows <- selected_rows[-idx]
selected_rows <- unique(c(selected_rows, "CD24_SIGLEC10", "C5AR1_RPS19", "HLA-F_B2M", "CXCL2_DPP4", "C1QB_LRP1", "PSAP_LRP1", "CD74_APP",
                   "INSR_NAMPT", "CXCL8_SDC2", "VEGFA_EPHB2", "SPP1_CD44", "EGFR_HSP90AA1", "ITGB1_CD14", "CD14_ITGA4", "APOE_LRP2", "APOE_LSR", "APOE_LRP1"))
a <- apply(tumor_TAM, 1, function(x){
    b <- unlist(str_split(as.character(x[4]), "\\|"))
    if((x[3]=="Ligand" & b[1]=="Tumor")|(x[3]=="Receptor" & b[2]=="Tumor")){
        return("ok")
    }else{
        return("error")
    }
})
selected_rows1 <- selected_rows[which(selected_rows %in% tumor_TAM$interacting_pair[which(a=="error")])]
selected_columns1 <- as.character(tumor_TAM$CC[which(tumor_TAM$interacting_pair %in% selected_rows1)])
selected_rows2 <- selected_rows[which(selected_rows %in% tumor_TAM$interacting_pair[which(a=="ok")])]
selected_columns2 <- as.character(tumor_TAM$CC[which(tumor_TAM$interacting_pair %in% selected_rows2)])

pdf("Tumor/Myeloid/tumor.TAM.interaction.receptor.pdf", width = 3.8, height = 4.2)
p1 <- dot_plot(selected_rows = selected_rows1, selected_columns = unique(selected_columns1), 
               breaks = c(-2, -1, 0, 1, 2), limit = c(-2, 2), legend.key.size = unit(0.1, "inches"), order.column = c("TAM-C1QB|Tumor", "TAM-RGCC|Tumor", "TAM-LGALS3|Tumor"),
               size = 8, max.dot.size = 4, means_path = "Renal/mymeans.usr.txt", pvalues_path = "Renal/mypvals.usr.txt")
print(p1)
dev.off()
pdf("Tumor/Myeloid/tumor.TAM.interaction.ligand.pdf", width = 3.8, height = 4.2)
p1 <- dot_plot(selected_rows = selected_rows2, selected_columns = unique(selected_columns2), 
               breaks = c(-2, -1, 0, 1, 2), limit = c(-2, 2), legend.key.size = unit(0.1, "inches"), order.column = c("Tumor|TAM-C1QB", "Tumor|TAM-RGCC", "Tumor|TAM-LGALS3"),
               size = 8, max.dot.size = 4, means_path = "Renal/mymeans.usr.txt", pvalues_path = "Renal/mypvals.usr.txt")
print(p1)
dev.off()

pdf("Tumor/Myeloid/SPP1.expression.pdf", width = 3, height = 3)
p <- DotPlot(data.merge.new, features = c("SPP1", "PTGER4", "ITGA9", "ITGA4"), cols = c("#1e90ff", "#F15F30"), dot.scale = 3.5) + xlab("") + ylab("") + theme(axis.text.x = element_text(size = 8, angle = 90)) + theme(axis.text.y = element_text(size = 8)) + NoLegend()
print(p)
dev.off()
pdf("Tumor/Myeloid/SPP1.expression-test.pdf", width = 3, height = 3)
p <- DotPlot(data.merge.new, features = c("SPP1", "PTGER4", "ITGA9", "ITGA4"), cols = c("#1e90ff", "#F15F30"), dot.scale = 3.5) + xlab("") + ylab("") + theme(axis.text.x = element_text(size = 8, angle = 90)) + theme(axis.text.y = element_text(size = 8))
print(p)
dev.off()
pdf("Tumor/Myeloid/interested.interaction.pairs.pdf")
p <- VlnPlot(data.merge.new, features = c("MIF", "CD74", "CD44", "SPP1"), pt.size = 0, ncol = 2) + NoLegend()
print(p)
p <- VlnPlot(data.merge.new, features = c("CD24", "SIGLEC10", "APOE", "TREM2"), pt.size = 0, ncol = 2)
print(p)
p <- VlnPlot(data.merge.new, features = c("VIM", "CD44", "RPS19", "C5AR1"), pt.size = 0, ncol = 2)
print(p)
p <- VlnPlot(data.merge.new, features = c("PKM", "CD44", "CXCL8", "SDC2"), pt.size = 0, ncol = 2)
print(p)
p <- VlnPlot(data.merge.new, features = c("FN1", "PLAUR", "CXCL8", "SDC2"), pt.size = 0, ncol = 2)
print(p)
dev.off()

features <- c("MIF",  "VEGFA", "GPC1", 
              "DPP4", "CD70", "CALM1", "APOE", "VLDLR") 
features.list <- as.list(features)
names(features.list) <- features
cox.res <- RCC.icb.analysis(signature.list = features.list, expresionMatrix = normalized_expression, clincal.info = patient.info.RNA)
dev.off()

pdf("Tumor/Myeloid/tumor-TAM.interaction.circular.pdf")
plot.circular(track.matrix = track.matrix, gene.pair = c("FN1", "PLAUR"), 
              mypvals = mypvals.usr, mymeans = mymeans.usr, 
              means.threshold = mean.threshold, pval.threshold,
              High.low = T, colors = c("#D95319", "#3B6793"))
plot.circular(track.matrix = track.matrix, gene.pair = c("VIM", "CD44"), 
              mypvals = mypvals.usr, mymeans = mymeans.usr, 
              means.threshold = mean.threshold, pval.threshold,
              High.low = T, colors = c("#D95319", "#3B6793"))
plot.circular(track.matrix = track.matrix, gene.pair = c("C5AR1", "RPS19"), 
              mypvals = mypvals.usr, mymeans = mymeans.usr, 
              means.threshold = mean.threshold, pval.threshold,
              High.low = T, colors = c("#D95319", "#3B6793"))              
dev.off()

###########################Tumor VS T cell
tumor_lymphoid <- tumor$All[grep("Proliferative T cell|Treg|Exhausted", tumor$All$CC),]
selected_rows <- unique(tumor_lymphoid$interacting_pair[which(tumor_lymphoid$means>=mean.threshold & tumor_lymphoid$pvals<pval.threshold)])

DEGs <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/5.Immune/CD8T/cellType.DE.rds")
T.DEGs <- DEGs[["Exhausted"]]  %>% filter(p_val_adj<0.05 & avg_log2FC>0.25)
all.DEGs <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/2.Cluster/AnnotateCellType/cellType.sig.pos.DEGs.rds")
other.DEGs <- all.DEGs %>% filter(cluster %in% c("Proliferative T cell", "Treg", "Tumor")) %>% filter(p_val_adj<0.05 & avg_log2FC>0.25)
DEGs <- unique(c(T.DEGs$gene, other.DEGs$gene))
interaction.genes <- unlist(str_split(selected_rows, pattern="_"))
sig.DEGs <- intersect(interaction.genes, DEGs)

interaction.info <- sapply(selected_rows, function(x){
    a <- unlist(str_split(x, "_"))
    first <- match(a[1], sig.DEGs)
    second <- match(a[2], sig.DEGs)
    if(!is.na(first)&!is.na(second)){
        return("ok")
    }else{
        return("no")
    }
})
idx <- which(interaction.info == "ok")
selected_rows <- selected_rows[idx]

idx <- which(selected_rows %in% c("ITGAV_VEGFA", "ITGB1_VCAM1", "ITGB2_VCAM1", "TIMP1_CD63", "SPP1_PTGER4", "ITGA5_SPP1", "ARPC5_LRP2", "CD3G_B2M", "SPP1_ITGA4",
                                  "EGFR_S100A4", "ERBB2_S100A4", "SPP1_ITGA4", "S100A8_CD69", "PKM_CD44", "CD74_MIF", "PTPRC_LGALS1",
                                  "PTPRC_LGALS1", "RPSA_LAMB2", "HSP90AA1_LRP1", "CD3D_B2M", "CD3G_B2M", "TFRC_B2M", "CD247_B2M", "EGFR_CALM2",
                                  "HLA-F_B2M", "MIF_CXCR4", "MIF_CD44", "CD2_CD58", "EGFR_MIF", "VIM_CD44", "HSPA8_LRP2", "CALR_HLA-F", "HMMR_CALM1",
                                  "CALM1_FAS", "CALM2_AQP1", "EGFR_HSP90AA1", "ERBB2_HSP90AA1", "TGFB1_CXCR4", "EGFR_CALM1", "EZR_VCAM1"))
selected_rows <- selected_rows[-idx]
selected_rows <- unique(c(selected_rows, "CD27_CD70", "CD2_CD59", "CD2_CD48", "NRP1_VEGFA", "SPP1_CD44", "HLA-F_B2M", "EGFR_MIF",
                  "CD74_APP", "VEGFA_FLT1", "LTBR_LTB", "ITGB1_SPP1", "TFRC_B2M", "EGFR_HSP90AA1"))

a <- apply(tumor_lymphoid, 1, function(x){
    b <- unlist(str_split(as.character(x[4]), "\\|"))
    if((x[3]=="Ligand" & b[1]=="Tumor")|(x[3]=="Receptor" & b[2]=="Tumor")){
        return("ok")
    }else{
        return("error")
    }
})
selected_rows1 <- selected_rows[which(selected_rows %in% tumor_lymphoid$interacting_pair[which(a=="error")])]
selected_columns1 <- as.character(tumor_lymphoid$CC[which(tumor_lymphoid$interacting_pair %in% selected_rows1)])
selected_rows2 <- selected_rows[which(selected_rows %in% tumor_lymphoid$interacting_pair[which(a=="ok")])]
selected_columns2 <- as.character(tumor_lymphoid$CC[which(tumor_lymphoid$interacting_pair %in% selected_rows2)])

pdf("Tumor/Lymphoid/tumor.Lymphoid.interaction.receptor.pdf", width = 3.8, height = 4)
p1 <- dot_plot(selected_rows = selected_rows1, selected_columns = unique(selected_columns1), 
               breaks = c(-2, -1, 0, 1, 2), limit = c(-2, 2), legend.key.size = unit(0.1, "inches"), order.column = c("CD8+ T-Exhausted|Tumor", "Proliferative T cell|Tumor", "Treg|Tumor"),
               size = 8, max.dot.size = 4, means_path = "Renal/mymeans.usr.txt", pvalues_path = "Renal/mypvals.usr.txt")
print(p1)
dev.off()
pdf("Tumor/Lymphoid/tumor.Lymphoid.interaction.ligand.pdf", width = 3.65, height = 4)
p1 <- dot_plot(selected_rows = selected_rows2, selected_columns = unique(selected_columns2), 
               breaks = c(-2, -1, 0, 1, 2), limit = c(-2, 2), legend.key.size = unit(0.1, "inches"), order.column = c("Tumor|CD8+ T-Exhausted", "Tumor|Proliferative T cell", "Tumor|Treg"),
               size = 8, max.dot.size = 4, means_path = "Renal/mymeans.usr.txt", pvalues_path = "Renal/mypvals.usr.txt")
print(p1)
dev.off()

pdf("Tumor/Lymphoid/tumor-lymphoid.interaction.circular.pdf")
plot.circular(track.matrix = track.matrix, gene.pair = c("TGFB1", "CXCR4"), 
              mypvals = mypvals.usr, mymeans = mymeans.usr, 
              means.threshold = mean.threshold, pval.threshold = pval.threshold,
              High.low = T, colors = c("#D95319", "#3B6793"))
plot.circular(track.matrix = track.matrix, gene.pair = c("CCL5", "SDC4"), 
              mypvals = mypvals.usr, mymeans = mymeans.usr, 
              means.threshold = mean.threshold, pval.threshold = pval.threshold,
              High.low = T, colors = c("#D95319", "#3B6793"))
plot.circular(track.matrix = track.matrix, gene.pair = c("CD27", "CD70"), 
              mypvals = mypvals.usr, mymeans = mymeans.usr, 
              means.threshold = mean.threshold, pval.threshold,
              High.low = T, colors = c("#D95319", "#3B6793"))
plot.circular(track.matrix = track.matrix, gene.pair = c("LTBR", "LTB"), 
              mypvals = mypvals.usr, mymeans = mymeans.usr, 
              means.threshold = mean.threshold, pval.threshold,
              High.low = T, colors = c("#D95319", "#3B6793"))                                                
dev.off()

pdf("Tumor/Lymphoid/interested.interaction.pairs.pdf")
p <- VlnPlot(data.merge.new, features = c("LTB", "LTBR", "CD27", "CD70"), pt.size = 0, ncol = 2) + NoLegend()
print(p)
dev.off()