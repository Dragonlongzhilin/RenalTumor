#' @description: plot the drug screening result

#### Lincs data
# /data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/Drug/drug.LINCS.ipynb
setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/Drug/LINCS")
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
drug.data <- read.table("drug.data.txt", header = T, stringsAsFactors = F, sep = "\t")
heatmap.data <- pivot_wider(drug.data, names_from = c(Drug.name, cellLine), values_from = ZSCORE)
heatmap.data <- as.data.frame(heatmap.data)
rownames(heatmap.data) <- heatmap.data[,1]
heatmap.data <- as.matrix(heatmap.data[,-1])
cols <- colorRamp2(c(-3, -1), c("blue", "white"))
col_split <- c(rep("HA1E", 11), "A549")
pdf("drug.heatmap.pdf")
Heatmap(heatmap.data, name = "Zscore", border_gp = gpar(col = "#fafafa"), column_split = col_split, cluster_rows = F, col = cols, cluster_columns = F, width = unit(6, "cm"), height = unit(2, "cm"))
dev.off()


#### other cell types
drug.data <- read.table("lincs.TF.negative.drug.txt", header = T, stringsAsFactors = F, sep = "\t")
#remove RBD chemical
idx <- grep("BRD|BG-|BDR-", drug.data$Drug.name)
drug.data <- drug.data[-idx,]
#Tumor cell line
cellLine.data <- read.table("/data/ExtraDisk/sdd/longzhilin/Data/drugData/CMAP_LINCS_2020/cellLine.data.txt", header = T, stringsAsFactors = F, sep = "\t")
cellLine.data <- cellLine.data[which(cellLine.data$cell_type == "tumor"), ]
idx <- match(drug.data$cellLine, cellLine.data$cell_iname)
tumor.cellLines <- cellLine.data[idx,]
drug.data <- drug.data[which(!is.na(idx)),]

heatmap.data <- pivot_wider(drug.data, names_from = c(Drug.name, cellLine), values_from = ZSCORE)
heatmap.data <- as.data.frame(heatmap.data)
rownames(heatmap.data) <- heatmap.data[,1]
heatmap.data <- as.matrix(heatmap.data[,-1])
cols <- colorRamp2(c(-5, -1), c("blue", "white"))
pdf("drug.allCellLine.heatmap.pdf")
Heatmap(t(heatmap.data), name = "Zscore", border_gp = gpar(col = "#fafafa"), row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6), cluster_rows = F, col = cols, cluster_columns = F, width = unit(2, "cm"))
dev.off()