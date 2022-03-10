#' @description: Endothelium cells

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
source(file = "/home/longzhilin/Analysis_Code/Combined.P.FC.R")
source(file = "/home/longzhilin/Analysis_Code/Visualization/colorPalettes.R")
source(file = "/home/longzhilin/Analysis_Code/code/ratio.plot.R")

scRNA.data <- readRDS("data.merge.pro.rds")

###########################1.differential gene#################################
Idents(scRNA.data) <- scRNA.data$cellType_low
idents <- as.character(levels(scRNA.data))
Endothelium.markers <- FindMarkers(scRNA.data, 
                                    ident.1 = "Endothelium (VCAM1+)",
                                    ident.2 = "Endothelium (VCAM1-)",
                                    logfc.threshold = 0, 
                                    min.pct = 0.1, 
                                    test.use = "MAST", 
                                    latent.vars = "orig.ident")
Endothelium.markers <- arrange(Endothelium.markers, desc(avg_log2FC))
Endothelium.markers$gene <- rownames(Endothelium.markers)
Endothelium.DEGs <- Endothelium.markers %>% filter(abs(avg_log2FC) >= 0.25 & p_val_adj < 0.05)
write.xlsx(Endothelium.markers, file = "2.Cluster/AnnotateCellType/Endothelium.all.DEGs.xlsx", rowNames = F, overwrite = T)
saveRDS(Endothelium.markers, file = "2.Cluster/AnnotateCellType/Endothelium.all.DEGs.rds")

#### volcano map
source("/home/longzhilin/Analysis_Code/Visualization/Plot.EnhancedVolcano.R")
pdf("2.Cluster/AnnotateCellType/Endothelium.EnhancedVolcano.pdf")
Plot.EnhancedVolcano(Endothelium.DEGs, x = "avg_log2FC", y = "p_val_adj", selected.showType = "Both", select.num = 15, title = "Endothelium (VCAM1+) VS Endothelium (VCAM1-)")
dev.off()

#### survival analysis
source(file = "/home/longzhilin/Analysis_Code/code/analysis.diff.survival.TCGA.R")
DESeq2.normalized_counts <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/DESeq2.normalized_counts.rds")
DESeq2.normalized_counts <- log2(DESeq2.normalized_counts+1)
DESeq2.result <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/DESeq2.result.rds")
clin.data <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/clin.data.rds")
pdf("2.Cluster/AnnotateCellType/Endothelium.signature.survival.pdf")
Endothelium.VCAM1.positive <- analysis.diff.survival.TCGA(interest.gene = Endothelium.DEGs$gene[which(Endothelium.DEGs$avg_log2FC > 0.5 & Endothelium.DEGs$p_val_adj < 0.05)], diff.gene.pro = DESeq2.result, exp.data.process = DESeq2.normalized_counts, clin.data = clin.data, EnhancedVolcano.plot = F, main = "Endothelium (VCAM1+)", Box.plot = F, meta.signature = T, single.signature = F)
Endothelium.VCAM1.negative <- analysis.diff.survival.TCGA(interest.gene = Endothelium.DEGs$gene[which(Endothelium.DEGs$avg_log2FC < -0.5 & Endothelium.DEGs$p_val_adj < 0.05)], diff.gene.pro = DESeq2.result, exp.data.process = DESeq2.normalized_counts, clin.data = clin.data, EnhancedVolcano.plot = F, main = "Endothelium (VCAM1-)", Box.plot = F, meta.signature = T, single.signature = F)
dev.off()

#### pathway enrichment
source("/home/longzhilin/Analysis_Code/PathwayEnrichment/clusterProfiler.enricher.R")
up.DEGs <- filter(Endothelium.DEGs, avg_log2FC > 1 & p_val_adj < 0.05)
down.DEGs <- filter(Endothelium.DEGs, avg_log2FC< -1 & p_val_adj < 0.05)
DEGs <- list(up = up.DEGs$gene, down = down.DEGs$gene)
pdf("2.Cluster/AnnotateCellType/Endothelium.program.pathway.enrichment.pdf")
res <- lapply(names(DEGs), function(x){
    y <- DEGs[[x]]
    res <- cluterProfiler.enricher(gene = y, geneType = "SYMBOL", db.type = "MsigDB", saveDir = paste0(getwd(),"/2.Cluster/AnnotateCellType"),
                            title = x, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
    # gene ratio
    res <- res$em.res.genesymbol@result %>% filter(p.adjust<0.05) #fdr adjust
    pathway.geneNum <- unlist(strsplit(res$BgRatio, "/"))[seq(1, nrow(res),by=2)]
    gene.ratio <- as.numeric(res$Count)/as.numeric(pathway.geneNum)
    res$gene.ratio <- gene.ratio
    return(res)
})
dev.off()

####cluster profiler enrichment result display
library(ggthemes)
#top10
enrich1 <- read.csv("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/2.Cluster/AnnotateCellType/up_enricherResult.csv", header = T)
enrich2 <- read.csv("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/2.Cluster/AnnotateCellType/down_enricherResult.csv", header = T)
enrich.list <- list(up = enrich1, down = enrich2)
enrich <- list(up = enrich1[c(1,3,8,9,16,17,18,21,22,23,25,26,27,28,42,49,50,72), c("ID")],
               down = enrich2[c(2:5,8:10,11,12,19,21,22), c("ID")])
enrich.path <- Reduce(function(x,y) c(x,y), enrich)
enrich.path <- unique(enrich.path)
enrich.path <- enrich.path[-c(3,6,26,28,29,4,13,11,5,20,23,27,30)]
enrich.list.select <- lapply(enrich.list, function(x){
  idx <- which(x$ID %in% enrich.path)
  return(x[idx, c("ID", "Count", "pvalue", "p.adjust")])
})
enrich.res <- Reduce(function(x,y) rbind(x,y), enrich.list.select)
a <- sapply(enrich.list.select, function(x){nrow(x)})
enrich.res$Type <- c(rep("Endothelium (VCAM1+) ", as.numeric(a[1])), rep("Endothelium (VCAM1-) ", as.numeric(a[2])))
enrich.res$p.adjust <- -log10(enrich.res$p.adjust)

#处理名字
enrich.res$ID <- gsub("_", " ", enrich.res$ID)
# factor ID
idx1 <- na.omit(match(enrich$up, enrich.path))
idx2 <- na.omit(match(enrich$down, enrich.path))
enrich.path <- enrich.path[c(idx1, idx2)]
enrich.res$ID <- factor(enrich.res$ID, levels = gsub("_", " ", enrich.path))
#source(file = "/home/longzhilin/Analysis_Code/Visualization/ggplot.bubble.plot.R")
pdf("2.Cluster/AnnotateCellType/Endothelium.enrichment.bubble.pdf", width = 5, height = 5)
p1 <- ggplot(enrich.res, 
            aes(x = Type, 
                y = ID, 
                size = Count, 
                fill = p.adjust))
p2 <- p1 + guides(color=FALSE) + geom_point(shape = 21, alpha = 0.7) + theme_few() + scale_size_continuous(breaks = c(3,5,7,9,11), labels = c(3,5,7,9,11)) + scale_fill_gradientn(colours = c("#fed71b", "#fcc41c", "#f27620", "#e92925"), breaks = c(1.69897, 4, 6, 8), labels = c(0.02, 0.0001, 0.000001, 0.00000001))
p2 <- p2 + theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1, size = 8), axis.text.y = element_text(size = 8)) + xlab("") + ylab("")
print(p2)
dev.off()

###########################2.interaction with immune cell and tumor cell#################################
source(file = "/home/longzhilin/Analysis_Code/code/cellphoneDB.dotplot.R")
source(file = "/home/longzhilin/Analysis_Code/code/Plot.CircularChart.R")
setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/6.CrossTalk/CellPhoneDB")
mean.threshold <- 1
pval.threshold <- 0.05
colors <- c("#00A087", "#4DBBD5", "#E64B35", "#3C5488")
mypvals.usr <- readRDS(file = "mypvals.usr.rds")
mymeans.usr <- readRDS(file = "mymeans.usr.rds")
sig.means.usr <- readRDS(file = "sig.means.usr.rds")
sig.interaction.cellType <- readRDS(file = "sig.interaction.cellType.rds")
interaction.counts <- readRDS("sig.interaction.counts.rds")
All.interaction.number <- readRDS("All.interaction.number.rds")

order.value <- c("CD4+ T cell", "Treg", "CD8+ T-Tissue-resident.C1", "CD8+ T-Tissue-resident.C2", "CD8+ T-IEG", "CD8+ T-Exhausted", "Proliferative T cell", "NK cell", "B cell",
                 "TAM-C1QB", "TAM-RGCC", "TAM-LGALS3", "Monocyte", "Dendritic cell", "Neutrophil", "Mast cell",
                 "Endothelium VCAM1+", "Endothelium VCAM1-", "Mesangial cell", "Tumor")

#### ----Endothelium VCAM1+
Endothelium.VCAM1.positive <- sig.interaction.cellType[["Endothelium VCAM1+"]]
#ligand
ligand.interact.Endothelium.VCAM1.positive <- apply(Endothelium.VCAM1.positive$Ligand, 1, function(x){
    a <- unlist(strsplit(x[4], "\\|"))
    return(a)
})
ligand.interact.Endothelium.VCAM1.positive.freq <- data.frame(table(as.character(ligand.interact.Endothelium.VCAM1.positive)))

idx <- which(ligand.interact.Endothelium.VCAM1.positive.freq$Var1 == "Endothelium VCAM1+")
ligand.interact.Endothelium.VCAM1.positive.freq[idx, 2] <- ligand.interact.Endothelium.VCAM1.positive.freq[idx, 2]-nrow(Endothelium.VCAM1.positive$Ligand)

receptor.interact.Endothelium.VCAM1.positive <- apply(Endothelium.VCAM1.positive$Receptor, 1, function(x){
    a <- unlist(strsplit(x[4], "\\|"))
    return(a)
})
receptor.interact.Endothelium.VCAM1.positive.freq <- data.frame(table(as.character(receptor.interact.Endothelium.VCAM1.positive)))
idx <- which(receptor.interact.Endothelium.VCAM1.positive.freq$Var1 == "Endothelium VCAM1+")
receptor.interact.Endothelium.VCAM1.positive.freq[idx, 2] <- receptor.interact.Endothelium.VCAM1.positive.freq[idx, 2]-nrow(Endothelium.VCAM1.positive$Receptor)

type <- rep("Lymphoid", nrow(receptor.interact.Endothelium.VCAM1.positive.freq))
type[c(7,10,12:13,16:18)] <- "Myeloid"
type[20] <- "Tumor"
type[c(8:9, 11)] <- "Other"
receptor.interact.Endothelium.VCAM1.positive.freq$Type <- factor(type, levels = c("Lymphoid", "Myeloid", "Tumor", "Other"))
ligand.interact.Endothelium.VCAM1.positive.freq$Type <- factor(type, levels = c("Lymphoid", "Myeloid", "Tumor", "Other"))

ligand.interact.Endothelium.VCAM1.positive.freq$Var1 <- factor(ligand.interact.Endothelium.VCAM1.positive.freq$Var1, levels = order.value)
receptor.interact.Endothelium.VCAM1.positive.freq$Var1 <- factor(receptor.interact.Endothelium.VCAM1.positive.freq$Var1, levels = order.value)
pdf("Endothelium/VCAM1.positive/Endothelium.VCAM1.positive.ligand.receptor.pdf")
ggbarplot(ligand.interact.Endothelium.VCAM1.positive.freq, x="Var1", y="Freq", color = "Type", xlab = "", palette = colors, ylab = "Number of Ligand", fill = "Type", sort.val = "desc", sort.by.groups = FALSE, x.text.angle = 60) 
ggbarplot(receptor.interact.Endothelium.VCAM1.positive.freq, x="Var1", y="Freq", color = "Type", xlab = "", palette = colors, ylab = "Number of Receptor", fill = "Type", sort.val = "desc", sort.by.groups = FALSE, x.text.angle = 60)
ggbarplot(ligand.interact.Endothelium.VCAM1.positive.freq[-c(8,9,11),], x="Var1", y="Freq", color = "Type", xlab = "", palette = colors[1:3], ylab = "Number of Ligand", fill = "Type", sort.val = "desc", sort.by.groups = FALSE, x.text.angle = 60) 
ggbarplot(receptor.interact.Endothelium.VCAM1.positive.freq[-c(8,9,11),], x="Var1", y="Freq", color = "Type", xlab = "", palette = colors[1:3], ylab = "Number of Receptor", fill = "Type", sort.val = "desc", sort.by.groups = FALSE, x.text.angle = 60)

ligand.receptor.data <- data.frame(Ligand = ligand.interact.Endothelium.VCAM1.positive.freq[,2], Receptor = receptor.interact.Endothelium.VCAM1.positive.freq[,2])
rownames(ligand.receptor.data) <- ligand.interact.Endothelium.VCAM1.positive.freq$Var1
row_split <- factor(ligand.interact.Endothelium.VCAM1.positive.freq$Type, levels = c("Lymphoid", "Myeloid", "Tumor", "Other"))
Heatmap(ligand.receptor.data, name = "Number of Interaction", col = colorRamp2(c(0,50), c("white", "red")), show_row_dend = F, show_column_dend = F, row_split = row_split, 
        width = unit(2, "cm"), height = unit(8, "cm"), cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.0f", ligand.receptor.data[i, j]), x, y, gp = gpar(fontsize = 7))}
        )

# total number
total.interaction <- data.frame(cellType = rownames(ligand.receptor.data), number = ligand.receptor.data[,1]+ligand.receptor.data[,2], type = ligand.interact.Endothelium.VCAM1.positive.freq$Type)
total.interaction$type <- factor(total.interaction$type, levels = c("Lymphoid", "Myeloid", "Tumor", "Other"))
ggbarplot(total.interaction, x="cellType", y="number", color = "type", xlab = "", palette = colors, ylab = "Number of interaction", fill = "type", sort.val = "desc", sort.by.groups = FALSE, x.text.angle = 60) 
dev.off()

#### Interaction between cancer cells and other cells
Endothelium.VCAM1.positive.partial <- Endothelium.VCAM1.positive$All[grep("Proliferative T cell|Treg|Exhausted|TAM|Tumor", Endothelium.VCAM1.positive$All$CC),]
selected_rows <- unique(Endothelium.VCAM1.positive.partial$interacting_pair[which(Endothelium.VCAM1.positive.partial$means>=mean.threshold & Endothelium.VCAM1.positive.partial$pvals<pval.threshold)])
a <- apply(Endothelium.VCAM1.positive.partial, 1, function(x){
    b <- unlist(str_split(as.character(x[4]), "\\|"))
    if((x[3]=="Ligand" & b[1]=="Endothelium VCAM1+")|(x[3]=="Receptor" & b[2]=="Endothelium VCAM1+")){
        return("ok")
    }else{
        return("error")
    }
})
# receptor
selected_rows1 <- selected_rows[which(selected_rows %in% Endothelium.VCAM1.positive.partial$interacting_pair[which(a=="error")])]
selected_columns1 <- as.character(Endothelium.VCAM1.positive.partial$CC[which(Endothelium.VCAM1.positive.partial$interacting_pair %in% selected_rows1)])
# ligand
selected_rows2 <- selected_rows[which(selected_rows %in% Endothelium.VCAM1.positive.partial$interacting_pair[which(a=="ok")])]
selected_columns2 <- as.character(Endothelium.VCAM1.positive.partial$CC[which(Endothelium.VCAM1.positive.partial$interacting_pair %in% selected_rows2)])

# filtering
selected_rows1 <- selected_rows1[grep("ACKR|VWF|VCAM1|CALM1|VEGFA|FN1|KDR",selected_rows1)]
selected_rows2 <- selected_rows2[grep("ACKR|VWF|VCAM1|CALM1|VEGFA|FN1",selected_rows2)]

pdf("Endothelium/VCAM1.positive/Endothelium.VCAM1.positive.interaction.receptor.pdf")
p1 <- dot_plot(selected_rows = selected_rows1[1:45], selected_columns = unique(selected_columns1), 
               breaks = c(-2, -1, 0, 1, 2), limit = c(-2, 2), legend.key.size = unit(0.1, "inches"), order.column = c("TAM-C1QB|Endothelium VCAM1+", "TAM-RGCC|Endothelium VCAM1+", "TAM-LGALS3|Endothelium VCAM1+", "CD8+ T-Exhausted|Endothelium VCAM1+", "Proliferative T cell|Endothelium VCAM1+", "Treg|Endothelium VCAM1+", "Tumor|Endothelium VCAM1+"),
               size = 8, max.dot.size = 4, means_path = "Renal/mymeans.usr.txt", pvalues_path = "Renal/mypvals.usr.txt")
print(p1)
p1 <- dot_plot(selected_rows = selected_rows1[46:90], selected_columns = unique(selected_columns1), 
               breaks = c(-2, -1, 0, 1, 2), limit = c(-2, 2), legend.key.size = unit(0.1, "inches"), order.column = c("TAM-C1QB|Endothelium VCAM1+", "TAM-RGCC|Endothelium VCAM1+", "TAM-LGALS3|Endothelium VCAM1+", "CD8+ T-Exhausted|Endothelium VCAM1+", "Proliferative T cell|Endothelium VCAM1+", "Treg|Endothelium VCAM1+", "Tumor|Endothelium VCAM1+"),
               size = 8, max.dot.size = 4, means_path = "Renal/mymeans.usr.txt", pvalues_path = "Renal/mypvals.usr.txt")
print(p1)
p1 <- dot_plot(selected_rows = selected_rows1[91:135], selected_columns = unique(selected_columns1), 
               breaks = c(-2, -1, 0, 1, 2), limit = c(-2, 2), legend.key.size = unit(0.1, "inches"), order.column = c("TAM-C1QB|Endothelium VCAM1+", "TAM-RGCC|Endothelium VCAM1+", "TAM-LGALS3|Endothelium VCAM1+", "CD8+ T-Exhausted|Endothelium VCAM1+", "Proliferative T cell|Endothelium VCAM1+", "Treg|Endothelium VCAM1+", "Tumor|Endothelium VCAM1+"),
               size = 8, max.dot.size = 4, means_path = "Renal/mymeans.usr.txt", pvalues_path = "Renal/mypvals.usr.txt")
print(p1)
dev.off()
pdf("Endothelium/VCAM1.positive/Endothelium.VCAM1.positive.interaction.ligand.pdf")
p1 <- dot_plot(selected_rows = selected_rows2[1:43], selected_columns = unique(selected_columns2), 
               breaks = c(-2, -1, 0, 1, 2), limit = c(-2, 2), legend.key.size = unit(0.1, "inches"), order.column = c("Endothelium VCAM1+|TAM-C1QB", "Endothelium VCAM1+|TAM-RGCC", "Endothelium VCAM1+|TAM-LGALS3", "Endothelium VCAM1+|CD8+ T-Exhausted", "Endothelium VCAM1+|Proliferative T cell", "Endothelium VCAM1+|Treg", "Endothelium VCAM1+|Tumor"),
               size = 8, max.dot.size = 4, means_path = "Renal/mymeans.usr.txt", pvalues_path = "Renal/mypvals.usr.txt")
print(p1)
p1 <- dot_plot(selected_rows = selected_rows2[44:86], selected_columns = unique(selected_columns2), 
               breaks = c(-2, -1, 0, 1, 2), limit = c(-2, 2), legend.key.size = unit(0.1, "inches"), order.column = c("Endothelium VCAM1+|TAM-C1QB", "Endothelium VCAM1+|TAM-RGCC", "Endothelium VCAM1+|TAM-LGALS3", "Endothelium VCAM1+|CD8+ T-Exhausted", "Endothelium VCAM1+|Proliferative T cell", "Endothelium VCAM1+|Treg", "Endothelium VCAM1+|Tumor"),
               size = 8, max.dot.size = 4, means_path = "Renal/mymeans.usr.txt", pvalues_path = "Renal/mypvals.usr.txt")
print(p1)
p1 <- dot_plot(selected_rows = selected_rows2[89:129], selected_columns = unique(selected_columns2), 
               breaks = c(-2, -1, 0, 1, 2), limit = c(-2, 2), legend.key.size = unit(0.1, "inches"), order.column = c("Endothelium VCAM1+|TAM-C1QB", "Endothelium VCAM1+|TAM-RGCC", "Endothelium VCAM1+|TAM-LGALS3", "Endothelium VCAM1+|CD8+ T-Exhausted", "Endothelium VCAM1+|Proliferative T cell", "Endothelium VCAM1+|Treg", "Endothelium VCAM1+|Tumor"),
               size = 8, max.dot.size = 4, means_path = "Renal/mymeans.usr.txt", pvalues_path = "Renal/mypvals.usr.txt")
print(p1)
dev.off()

#### ----Endothelium VCAM1-
Endothelium.VCAM1.negative <- sig.interaction.cellType[["Endothelium VCAM1-"]]
#ligand
ligand.interact.Endothelium.VCAM1.negative <- apply(Endothelium.VCAM1.negative$Ligand, 1, function(x){
    a <- unlist(strsplit(x[4], "\\|"))
    return(a)
})
ligand.interact.Endothelium.VCAM1.negative.freq <- data.frame(table(as.character(ligand.interact.Endothelium.VCAM1.negative)))

idx <- which(ligand.interact.Endothelium.VCAM1.negative.freq$Var1 == "Endothelium VCAM1-")
ligand.interact.Endothelium.VCAM1.negative.freq[idx, 2] <- ligand.interact.Endothelium.VCAM1.negative.freq[idx, 2]-nrow(Endothelium.VCAM1.negative$Ligand)

receptor.interact.Endothelium.VCAM1.negative <- apply(Endothelium.VCAM1.negative$Receptor, 1, function(x){
    a <- unlist(strsplit(x[4], "\\|"))
    return(a)
})
receptor.interact.Endothelium.VCAM1.negative.freq <- data.frame(table(as.character(receptor.interact.Endothelium.VCAM1.negative)))
idx <- which(receptor.interact.Endothelium.VCAM1.negative.freq$Var1 == "Endothelium VCAM1-")
receptor.interact.Endothelium.VCAM1.negative.freq[idx, 2] <- receptor.interact.Endothelium.VCAM1.negative.freq[idx, 2]-nrow(Endothelium.VCAM1.negative$Receptor)

type <- rep("Lymphoid", nrow(receptor.interact.Endothelium.VCAM1.negative.freq))
type[c(7,10,12:13,16:18)] <- "Myeloid"
type[20] <- "Tumor"
type[c(8:9, 11)] <- "Other"
receptor.interact.Endothelium.VCAM1.negative.freq$Type <- factor(type, levels = c("Lymphoid", "Myeloid", "Tumor", "Other"))
ligand.interact.Endothelium.VCAM1.negative.freq$Type <- factor(type, levels = c("Lymphoid", "Myeloid", "Tumor", "Other"))

ligand.interact.Endothelium.VCAM1.negative.freq$Var1 <- factor(ligand.interact.Endothelium.VCAM1.negative.freq$Var1, levels = order.value)
receptor.interact.Endothelium.VCAM1.negative.freq$Var1 <- factor(receptor.interact.Endothelium.VCAM1.negative.freq$Var1, levels = order.value)
pdf("Endothelium/VCAM1.negative/Endothelium.VCAM1.negative.ligand.receptor.pdf")
ggbarplot(ligand.interact.Endothelium.VCAM1.negative.freq, x="Var1", y="Freq", color = "Type", xlab = "", palette = colors, ylab = "Number of Ligand", fill = "Type", sort.val = "desc", sort.by.groups = FALSE, x.text.angle = 60) 
ggbarplot(receptor.interact.Endothelium.VCAM1.negative.freq, x="Var1", y="Freq", color = "Type", xlab = "", palette = colors, ylab = "Number of Receptor", fill = "Type", sort.val = "desc", sort.by.groups = FALSE, x.text.angle = 60)
ggbarplot(ligand.interact.Endothelium.VCAM1.negative.freq[-c(8,9,11),], x="Var1", y="Freq", color = "Type", xlab = "", palette = colors[1:3], ylab = "Number of Ligand", fill = "Type", sort.val = "desc", sort.by.groups = FALSE, x.text.angle = 60) 
ggbarplot(receptor.interact.Endothelium.VCAM1.negative.freq[-c(8,9,11),], x="Var1", y="Freq", color = "Type", xlab = "", palette = colors[1:3], ylab = "Number of Receptor", fill = "Type", sort.val = "desc", sort.by.groups = FALSE, x.text.angle = 60)

ligand.receptor.data <- data.frame(Ligand = ligand.interact.Endothelium.VCAM1.negative.freq[,2], Receptor = receptor.interact.Endothelium.VCAM1.negative.freq[,2])
rownames(ligand.receptor.data) <- ligand.interact.Endothelium.VCAM1.negative.freq$Var1
row_split <- factor(ligand.interact.Endothelium.VCAM1.negative.freq$Type, levels = c("Lymphoid", "Myeloid", "Tumor", "Other"))
Heatmap(ligand.receptor.data, name = "Number of Interaction", col = colorRamp2(c(0,50), c("white", "red")), show_row_dend = F, show_column_dend = F, row_split = row_split, 
        width = unit(2, "cm"), height = unit(8, "cm"), cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.0f", ligand.receptor.data[i, j]), x, y, gp = gpar(fontsize = 7))}
        )

# total number
total.interaction <- data.frame(cellType = rownames(ligand.receptor.data), number = ligand.receptor.data[,1]+ligand.receptor.data[,2], type = ligand.interact.Endothelium.VCAM1.negative.freq$Type)
total.interaction$type <- factor(total.interaction$type, levels = c("Lymphoid", "Myeloid", "Tumor", "Other"))
ggbarplot(total.interaction, x="cellType", y="number", color = "type", xlab = "", palette = colors, ylab = "Number of interaction", fill = "type", sort.val = "desc", sort.by.groups = FALSE, x.text.angle = 60) 
dev.off()

#### Interaction between cancer cells and other cells
Endothelium.VCAM1.negative.partial <- Endothelium.VCAM1.negative$All[grep("Proliferative T cell|Treg|Exhausted|TAM|Tumor", Endothelium.VCAM1.negative$All$CC),]
selected_rows <- unique(Endothelium.VCAM1.negative.partial$interacting_pair[which(Endothelium.VCAM1.negative.partial$means>=mean.threshold & Endothelium.VCAM1.negative.partial$pvals<pval.threshold)])
a <- apply(Endothelium.VCAM1.negative.partial, 1, function(x){
    b <- unlist(str_split(as.character(x[4]), "\\|"))
    if((x[3]=="Ligand" & b[1]=="Endothelium VCAM1-")|(x[3]=="Receptor" & b[2]=="Endothelium VCAM1-")){
        return("ok")
    }else{
        return("error")
    }
})
# receptor
selected_rows1 <- selected_rows[which(selected_rows %in% Endothelium.VCAM1.negative.partial$interacting_pair[which(a=="error")])]
selected_columns1 <- as.character(Endothelium.VCAM1.negative.partial$CC[which(Endothelium.VCAM1.negative.partial$interacting_pair %in% selected_rows1)])
# ligand
selected_rows2 <- selected_rows[which(selected_rows %in% Endothelium.VCAM1.negative.partial$interacting_pair[which(a=="ok")])]
selected_columns2 <- as.character(Endothelium.VCAM1.negative.partial$CC[which(Endothelium.VCAM1.negative.partial$interacting_pair %in% selected_rows2)])

pdf("Endothelium/VCAM1.negative/Endothelium.VCAM1.negative.interaction.receptor.pdf")
p1 <- dot_plot(selected_rows = selected_rows1[1:45], selected_columns = unique(selected_columns1), 
               breaks = c(-2, -1, 0, 1, 2), limit = c(-2, 2), legend.key.size = unit(0.1, "inches"), order.column = c("TAM-C1QB|Endothelium VCAM1-", "TAM-RGCC|Endothelium VCAM1-", "TAM-LGALS3|Endothelium VCAM1-", "CD8+ T-Exhausted|Endothelium VCAM1-", "Proliferative T cell|Endothelium VCAM1-", "Treg|Endothelium VCAM1-", "Tumor|Endothelium VCAM1-"),
               size = 8, max.dot.size = 4, means_path = "Renal/mymeans.usr.txt", pvalues_path = "Renal/mypvals.usr.txt")
print(p1)
p1 <- dot_plot(selected_rows = selected_rows1[46:90], selected_columns = unique(selected_columns1), 
               breaks = c(-2, -1, 0, 1, 2), limit = c(-2, 2), legend.key.size = unit(0.1, "inches"), order.column = c("TAM-C1QB|Endothelium VCAM1-", "TAM-RGCC|Endothelium VCAM1-", "TAM-LGALS3|Endothelium VCAM1-", "CD8+ T-Exhausted|Endothelium VCAM1-", "Proliferative T cell|Endothelium VCAM1-", "Treg|Endothelium VCAM1-", "Tumor|Endothelium VCAM1-"),
               size = 8, max.dot.size = 4, means_path = "Renal/mymeans.usr.txt", pvalues_path = "Renal/mypvals.usr.txt")
print(p1)
p1 <- dot_plot(selected_rows = selected_rows1[91:135], selected_columns = unique(selected_columns1), 
               breaks = c(-2, -1, 0, 1, 2), limit = c(-2, 2), legend.key.size = unit(0.1, "inches"), order.column = c("TAM-C1QB|Endothelium VCAM1-", "TAM-RGCC|Endothelium VCAM1-", "TAM-LGALS3|Endothelium VCAM1-", "CD8+ T-Exhausted|Endothelium VCAM1-", "Proliferative T cell|Endothelium VCAM1-", "Treg|Endothelium VCAM1-", "Tumor|Endothelium VCAM1-"),
               size = 8, max.dot.size = 4, means_path = "Renal/mymeans.usr.txt", pvalues_path = "Renal/mypvals.usr.txt")
print(p1)
dev.off()
pdf("Endothelium/VCAM1.negative/Endothelium.VCAM1.negative.interaction.ligand.pdf")
p1 <- dot_plot(selected_rows = selected_rows2[1:40], selected_columns = unique(selected_columns2), 
               breaks = c(-2, -1, 0, 1, 2), limit = c(-2, 2), legend.key.size = unit(0.1, "inches"), order.column = c("Endothelium VCAM1-|TAM-C1QB", "Endothelium VCAM1-|TAM-RGCC", "Endothelium VCAM1-|TAM-LGALS3", "Endothelium VCAM1-|CD8+ T-Exhausted", "Endothelium VCAM1-|Proliferative T cell", "Endothelium VCAM1-|Treg", "Endothelium VCAM1-|Tumor"),
               size = 8, max.dot.size = 4, means_path = "Renal/mymeans.usr.txt", pvalues_path = "Renal/mypvals.usr.txt")
print(p1)
p1 <- dot_plot(selected_rows = selected_rows2[41:80], selected_columns = unique(selected_columns2), 
               breaks = c(-2, -1, 0, 1, 2), limit = c(-2, 2), legend.key.size = unit(0.1, "inches"), order.column = c("Endothelium VCAM1-|TAM-C1QB", "Endothelium VCAM1-|TAM-RGCC", "Endothelium VCAM1-|TAM-LGALS3", "Endothelium VCAM1-|CD8+ T-Exhausted", "Endothelium VCAM1-|Proliferative T cell", "Endothelium VCAM1-|Treg", "Endothelium VCAM1-|Tumor"),
               size = 8, max.dot.size = 4, means_path = "Renal/mymeans.usr.txt", pvalues_path = "Renal/mypvals.usr.txt")
print(p1)
p1 <- dot_plot(selected_rows = selected_rows2[81:122], selected_columns = unique(selected_columns2), 
               breaks = c(-2, -1, 0, 1, 2), limit = c(-2, 2), legend.key.size = unit(0.1, "inches"), order.column = c("Endothelium VCAM1-|TAM-C1QB", "Endothelium VCAM1-|TAM-RGCC", "Endothelium VCAM1-|TAM-LGALS3", "Endothelium VCAM1-|CD8+ T-Exhausted", "Endothelium VCAM1-|Proliferative T cell", "Endothelium VCAM1-|Treg", "Endothelium VCAM1-|Tumor"),
               size = 8, max.dot.size = 4, means_path = "Renal/mymeans.usr.txt", pvalues_path = "Renal/mypvals.usr.txt")
print(p1)
dev.off()