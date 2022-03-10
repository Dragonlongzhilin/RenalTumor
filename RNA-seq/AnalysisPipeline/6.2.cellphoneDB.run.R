#' @description: run cellphoneDB
python -m venv cpdb-venv
source cpdb-venv/bin/activate
pip install cellphonedb

library(Seurat)
setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA")
####RNA count matrix
data.merge <- readRDS("data.merge.pro.rds")
Macrophage <- readRDS("5.Immune/Macrophage/sub.scRNA.harmony.pro.rds") # 3 cluster
CD8.T <- readRDS("5.Immune/CD8T/sub.scRNA.harmony.pro.rds")
idx1 <- match(colnames(Macrophage), colnames(data.merge)) # 4319 cells
idx2 <- match(colnames(CD8.T), colnames(data.merge)) # 3728 cells
cell.Type <- data.merge$cellType_low
cell.Type[idx1] <- as.character(Macrophage$cellType3)
cell.Type[idx2] <- paste0("CD8+ T-", CD8.T$cellType3)
data.merge$cellType2 <- cell.Type

cell.Type.new <- setdiff(unique(cell.Type), c("CD8+ T cell", "Macrophage"))
data.merge.new <- subset(data.merge, subset = cellType2 %in% cell.Type.new)
DefaultAssay(data.merge.new) <- "RNA"
#In order to avoid the trouble of handling CD8+ T-Exhausted and CD8+ T-Exhausted IEG, change CD8+ T-Exhausted IEG --- CD8+ T-IEG
data.merge.new$cellType2 <- gsub("CD8\\+ T-Exhausted IEG", "CD8\\+ T-IEG", data.merge.new$cellType2)

RNA.normalizedCount <- GetAssayData(data.merge.new, slot = "data") #normalised counts
write.table(as.matrix(RNA.normalizedCount), file = "6.CrossTalk/CellPhoneDB/RNA.normailzedCount.txt", sep = "\t", quote = F)
meta_data <- cbind(rownames(data.merge.new@meta.data), data.merge.new@meta.data[,'cellType2', drop=F])  
meta_data <- as.matrix(meta_data)
write.table(meta_data, '6.CrossTalk/CellPhoneDB/RNA.cell_meta.txt', sep='\t', quote=F, row.names=F)

#python code
cd /data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/6.CrossTalk/CellPhoneDB
cellphonedb method statistical_analysis --counts-data=gene_name --project-name=Renal --threshold 0.1 --database /data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/ligand_receptor/cellphoneDB/out/cellphonedb_user_2021-08-02-20_16.db --output-path=/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/6.CrossTalk/CellPhoneDB --threads=48 RNA.cell_meta.txt RNA.normailzedCount.txt
cellphonedb plot dot_plot --means-path=./Renal/means.txt --pvalues-path=./Renal/pvalues.txt --output-path=/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/6.CrossTalk/CellPhoneDB
cellphonedb plot heatmap_plot --pvalues-path=./Renal/pvalues.txt --output-path=/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/6.CrossTalk/CellPhoneDB RNA.cell_meta.txt

