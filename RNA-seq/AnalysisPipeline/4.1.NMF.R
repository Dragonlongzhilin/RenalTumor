#' @description: NMF analysis for tumor cell

seed.data <- 101
set.seed(seed.data)
library(NMF) #0.23.0
library(ggpubr)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(corrplot)
library(cowplot)
library(dendsort)
library(openxlsx)
setwd(dir = "/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/4.NMF/SCT")

#### data preparation
library(future)
plan("multiprocess", workers = 5)
options(future.globals.maxSize = 50000 * 1024^2) # set 50G RAM
# Since the similarity of the cells needs to be evaluated in the overall scale expression, the overall scale needs to be performed
data.merge <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/data.merge.pro.rds")
DefaultAssay(data.merge) <- "RNA"
CancerCell <- subset(data.merge, subset = cellType == "Tumor")
source(file = "/data/active_data/lzl/RenalTumor-20200713/code/Filter.gene.R")
CancerCell <- Filter.gene(CancerCell) #Remove genes less than 3 cells expressed
# scale.data: do.scale = T, do.centre = T
CancerCell <- SCTransform(CancerCell, vars.to.regress = c("nCount_RNA", "percent.mt"), return.only.var.genes = F, verbose = FALSE)
saveRDS(CancerCell, file = "CancerCell.rds")

stallion = c("1"="#EF7F48","2"="#D69100","3"="#C69900","4"="#83AD00","5"="#00BE6B", "6"="#00C0BA","7"="#00B9E1","8"="#00B3F1","19"="#E6C2DC", "10"="#8794FF", "11"="#B086FF","12"="#E46DF6","13"="#F365E6","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#00ABFD","20"="#3D3D3D")
## NMF, only analyze the cancer cell cluster, calculate in a single sample
## NMF requires that the input data has been standardized
DefaultAssay(data.merge) <- "RNA"
T1.cancer.data <- subset(data.merge, subset = cellType == "Tumor" & orig.ident == "T1") # 606 cells
T1.cancer.data <- SCTransform(T1.cancer.data, vars.to.regress = c("nCount_RNA", "percent.mt"), return.only.var.genes = F, verbose = FALSE)

T2.cancer.data <- subset(data.merge, subset = cellType == "Tumor" & orig.ident == "T2") # 2382 cells
T2.cancer.data <- SCTransform(T2.cancer.data, vars.to.regress = c("nCount_RNA", "percent.mt"), return.only.var.genes = F, verbose = FALSE)

T3.cancer.data <- subset(data.merge, subset = cellType == "Tumor" & orig.ident == "T3") # 320 cells
T3.cancer.data <- SCTransform(T3.cancer.data, vars.to.regress = c("nCount_RNA", "percent.mt"), return.only.var.genes = F, verbose = FALSE)

T4.cancer.data <- subset(data.merge, subset = cellType == "Tumor" & orig.ident == "T4") # 256 cells
T4.cancer.data <- SCTransform(T4.cancer.data, vars.to.regress = c("nCount_RNA", "percent.mt"), return.only.var.genes = F, verbose = FALSE)

saveRDS(T1.cancer.data, file = "T1.cancer.data.rds")
saveRDS(T2.cancer.data, file = "T2.cancer.data.rds")
saveRDS(T3.cancer.data, file = "T3.cancer.data.rds")
saveRDS(T4.cancer.data, file = "T4.cancer.data.rds")

#######################################################1.NMF calculate###################################################
T1.cancer.data <- readRDS("T1.cancer.data.rds") # 14896*606
T2.cancer.data <- readRDS("T2.cancer.data.rds") # 17153*2382
T3.cancer.data <- readRDS("T3.cancer.data.rds") # 14685*320
T4.cancer.data <- readRDS("T4.cancer.data.rds") # 12468*256
# Convert the negative value in the scale.data data to a value of 0, and it is required to have expression in at least 10 cells
Filter.gene <- function(data, min.cell = 10){
    exp.num <- apply(data, 1, function(x){
        return(length(which(x>0)))
    })
    data  <- data[which(exp.num>=min.cell),]
    return(data)
}
T1.cancer.scale.data <- GetAssayData(T1.cancer.data, slot = "scale.data")
T1.cancer.scale.data[which(T1.cancer.scale.data < 0)] <- 0
T1.cancer.scale.data <- Filter.gene(T1.cancer.scale.data) # 14634*606
T2.cancer.scale.data <- GetAssayData(T2.cancer.data, slot = "scale.data") 
T2.cancer.scale.data[which(T2.cancer.scale.data < 0)] <- 0 
T2.cancer.scale.data <- Filter.gene(T2.cancer.scale.data) # 16831*2382
T3.cancer.scale.data <- GetAssayData(T3.cancer.data, slot = "scale.data") 
T3.cancer.scale.data[which(T3.cancer.scale.data < 0)] <- 0
T3.cancer.scale.data <- Filter.gene(T3.cancer.scale.data) # 14390*320
T4.cancer.scale.data <- GetAssayData(T4.cancer.data, slot = "scale.data") 
T4.cancer.scale.data[which(T4.cancer.scale.data < 0)] <- 0 
T4.cancer.scale.data <- Filter.gene(T4.cancer.scale.data) # 12042*256

####stardard vatiation analysis
T1.sd <- apply(T1.cancer.scale.data, 1, sd)
T2.sd <- apply(T2.cancer.scale.data, 1, sd)
T3.sd <- apply(T3.cancer.scale.data, 1, sd)
T4.sd <- apply(T4.cancer.scale.data, 1, sd)

sd.data <- data.frame(Sample = c(rep("T1", length(T1.sd)), rep("T2", length(T2.sd)),rep("T3", length(T3.sd)), rep("T4", length(T4.sd))),
                      sd = c(T1.sd, T2.sd, T3.sd, T4.sd))
pdf("SD/sample.sd.distribution.pdf")
ggdensity(sd.data, x="sd", rug = TRUE, color = "Sample", fill = "Sample", palette = "npg") + geom_vline(xintercept = 0.5)
dev.off()
T1.cancer.scale.data.pro <- T1.cancer.scale.data[names(T1.sd)[which(T1.sd >= 0.5)],]
T2.cancer.scale.data.pro <- T2.cancer.scale.data[names(T2.sd)[which(T2.sd >= 0.5)],]
T3.cancer.scale.data.pro <- T3.cancer.scale.data[names(T3.sd)[which(T3.sd >= 0.5)],]
T4.cancer.scale.data.pro <- T4.cancer.scale.data[names(T4.sd)[which(T4.sd >= 0.5)],]

#### Determine the best rank
T1.nmf <- nmf(T1.cancer.scale.data.pro, 2:6, seed=seed.data, nrun = 30, .options = "vP48")
T2.nmf <- nmf(T2.cancer.scale.data.pro, 2:6, seed=seed.data, nrun = 30, .options = "vP48")
T3.nmf <- nmf(T3.cancer.scale.data.pro, 2:6, seed=seed.data, nrun = 30, .options = "vP48")
T4.nmf <- nmf(T4.cancer.scale.data.pro, 2:6, seed=seed.data, nrun = 30, .options = "vP48")
saveRDS(T1.nmf, "T1.nmf.rds")
saveRDS(T2.nmf, "T2.nmf.rds")
saveRDS(T3.nmf, "T3.nmf.rds")
saveRDS(T4.nmf, "T4.nmf.rds")
pdf("patient.optiaml.rank.pdf")
plot(T1.nmf)
plot(2:6,T1.nmf$measures$cophenetic, type="b", col="purple")
plot(T2.nmf)
plot(2:6,T2.nmf$measures$cophenetic, type="b", col="purple")
plot(T3.nmf)
plot(2:6,T3.nmf$measures$cophenetic, type="b", col="purple")
plot(T4.nmf)
plot(2:6,T4.nmf$measures$cophenetic, type="b", col="purple")
dev.off()

# function: Determine the best rank
optimal.rank <- function(nmf){
    coph <- nmf$measures$cophenetic
    coph_diff <- NULL
    for (i in 2:length(coph)){
        coph_diff <- c(coph_diff, coph[i-1]-coph[i])
    }
    k.best <- which.max(coph_diff)+1
    return(k.best)
}
optimal.rank(T1.nmf)
optimal.rank(T2.nmf)
optimal.rank(T3.nmf)
optimal.rank(T4.nmf)

#optimal rank
T1.nmf.optimal <- nmf(T1.cancer.scale.data.pro, optimal.rank(T1.nmf), nrun = 30, seed=seed.data, .options = "vP48")
T2.nmf.optimal <- nmf(T2.cancer.scale.data.pro, optimal.rank(T2.nmf), nrun = 30, seed=seed.data, .options = "vP48")
T3.nmf.optimal <- nmf(T3.cancer.scale.data.pro, optimal.rank(T3.nmf), nrun = 30, seed=seed.data, .options = "vP48")
T4.nmf.optimal <- nmf(T4.cancer.scale.data.pro, optimal.rank(T4.nmf), nrun = 30, seed=seed.data, .options = "vP48")

pdf("nmf.map.pdf")
coefmap(T1.nmf.optimal)
basismap(T1.nmf.optimal)
consensusmap(T1.nmf.optimal)
coefmap(T2.nmf.optimal)
basismap(T2.nmf.optimal)
consensusmap(T2.nmf.optimal)
coefmap(T3.nmf.optimal)
basismap(T3.nmf.optimal)
consensusmap(T3.nmf.optimal)
coefmap(T4.nmf.optimal)
basismap(T4.nmf.optimal)
consensusmap(T4.nmf.optimal)
dev.off()

saveRDS(T1.nmf.optimal, "T1.nmf.optimal.rds")
saveRDS(T2.nmf.optimal, "T2.nmf.optimal.rds")
saveRDS(T3.nmf.optimal, "T3.nmf.optimal.rds")
saveRDS(T4.nmf.optimal, "T4.nmf.optimal.rds")

#######################################################2.NMF anlaysis###################################################
#### load data
T1.cancer.data <- readRDS("T1.cancer.data.rds")
T2.cancer.data <- readRDS("T2.cancer.data.rds")
T3.cancer.data <- readRDS("T3.cancer.data.rds")
T4.cancer.data <- readRDS("T4.cancer.data.rds")
CancerCell <- readRDS("CancerCell.rds")

setwd(paste0(getwd(), "/SD"))
T1.nmf.optimal <- readRDS("T1.nmf.optimal.rds")
T2.nmf.optimal <- readRDS("T2.nmf.optimal.rds")
T3.nmf.optimal <- readRDS("T3.nmf.optimal.rds")
T4.nmf.optimal <- readRDS("T4.nmf.optimal.rds")

# Calculate the correlation between each sample and each program
T1.w <- basis(T1.nmf.optimal)
T1.H <- coef(T1.nmf.optimal)
T2.w <- basis(T2.nmf.optimal)
T2.H <- coef(T2.nmf.optimal)
T3.w <- basis(T3.nmf.optimal)
T3.H <- coef(T3.nmf.optimal)
T4.w <- basis(T4.nmf.optimal)
T4.H <- coef(T4.nmf.optimal)

T1.basis.score <- as.numeric(T1.w)
T2.basis.score <- as.numeric(T2.w)
T3.basis.score <- as.numeric(T3.w)
T4.basis.score <- as.numeric(T4.w)
basis.score <- data.frame(patient = c(rep("T1", length(T1.basis.score)),
                                      rep("T2", length(T2.basis.score)),
                                      rep("T3", length(T3.basis.score)),
                                      rep("T4", length(T4.basis.score))), 
                          basis.score = c(T1.basis.score, T2.basis.score, T3.basis.score, T4.basis.score))
pdf("basis.score.pdf")
ggdensity(basis.score, x = "basis.score", add = "mean", rug = TRUE, color = "patient", fill = "patient")
dev.off()

#####################################################2.identify program/meta signature###################################################
#############2.1 program signature
####Extract the signature gene corresponding to each patient program
T1.program.orders <- extractFeatures(T1.nmf.optimal, 30) 
T1.program.signature <- lapply(T1.program.orders, function(x){
    return(rownames(T1.w)[x])
})
T2.program.orders <- extractFeatures(T2.nmf.optimal, 30)
T2.program.signature <- lapply(T2.program.orders, function(x){
    return(rownames(T2.w)[x])
})
T3.program.orders <- extractFeatures(T3.nmf.optimal, 30)
T3.program.signature <- lapply(T3.program.orders, function(x){
    return(rownames(T3.w)[x])
})
T4.program.orders <- extractFeatures(T4.nmf.optimal, 30)
T4.program.signature <- lapply(T4.program.orders, function(x){
    return(rownames(T4.w)[x])
})

colnames(T1.w) <- paste0("T1_Program", 1:ncol(T1.w))
colnames(T2.w) <- paste0("T2_Program", 1:ncol(T2.w))
colnames(T3.w) <- paste0("T3_Program", 1:ncol(T3.w))
colnames(T4.w) <- paste0("T4_Program", 1:ncol(T4.w))
NMF.score <- list(T1 = T1.w, T2 = T2.w, T3 = T3.w, T4 = T4.w)
saveRDS(NMF.score, file = "NMF.score.rds")

T1.data.scale <- GetAssayData(T1.cancer.data, slot = "scale.data")
T2.data.scale <- GetAssayData(T2.cancer.data, slot = "scale.data")
T3.data.scale <- GetAssayData(T3.cancer.data, slot = "scale.data")
T4.data.scale <- GetAssayData(T4.cancer.data, slot = "scale.data")
pdf("program.signatures.heatmap.eachPatient.pdf")
col.list <- list(Phase = c("S" = "#3BBCA8", "G1" = "#6E4B9E", "G2M" = "#F59899"))
row.list <- list(Type = c("Program1"="#208A42", "Program2"="#7700FF","Program3"="#89288F", "Program4"="#e0598b","Program5"="#F47D2B"))
anno.info <- data.frame(Type = rep(paste0("Program", 1:length(T1.program.orders)), sapply(T1.program.orders, function(x){length(x)})))
row.anno <- HeatmapAnnotation(Type = anno.info$Type, which = "row", col = row.list)
top.anno <- HeatmapAnnotation(Phase = T1.cancer.data$Phase, 
                              col = col.list)
p <- Heatmap(T1.data.scale[unlist(T1.program.signature),], name = "Expression", cluster_rows = F, top_annotation = top.anno, left_annotation = row.anno, show_row_dend = F, show_column_dend = F, show_row_names = F, show_column_names = F)
print(p)

anno.info <- data.frame(Type = rep(paste0("Program", 1:length(T2.program.orders)), sapply(T2.program.orders, function(x){length(x)})))
row.anno <- HeatmapAnnotation(Type = anno.info$Type, which = "row", col = row.list)
top.anno <- HeatmapAnnotation(Phase = T2.cancer.data$Phase, 
                              col = col.list)
p <- Heatmap(T2.data.scale[unlist(T2.program.signature),], name = "Expression", cluster_rows = F, top_annotation = top.anno, left_annotation = row.anno, show_row_dend = F, show_column_dend = F, show_row_names = F, show_column_names = F)
print(p)

anno.info <- data.frame(Type = rep(paste0("Program", 1:length(T3.program.orders)), sapply(T3.program.orders, function(x){length(x)})))
row.anno <- HeatmapAnnotation(Type = anno.info$Type, which = "row", col = row.list)
top.anno <- HeatmapAnnotation(Phase = T3.cancer.data$Phase, 
                              col = col.list)
p <- Heatmap(T3.data.scale[unlist(T3.program.signature),], name = "Expression", cluster_rows = F, top_annotation = top.anno, left_annotation = row.anno, show_row_dend = F, show_column_dend = F, show_row_names = F, show_column_names = F)
print(p)

anno.info <- data.frame(Type = rep(paste0("Program", 1:length(T4.program.orders)), sapply(T4.program.orders, function(x){length(x)})))
row.anno <- HeatmapAnnotation(Type = anno.info$Type, which = "row", col = row.list)
top.anno <- HeatmapAnnotation(Phase = T4.cancer.data$Phase, 
                              col = col.list)
p <- Heatmap(T4.data.scale[unlist(T4.program.signature),], name = "Expression", cluster_rows = F, top_annotation = top.anno, left_annotation = row.anno, show_row_dend = F, show_column_dend = F, show_row_names = F, show_column_names = F)
print(p)
dev.off()

#####################################Calculate the correlation of each program
names(T1.program.signature) <- paste0("T1_Program", 1:length(T1.program.signature))
names(T2.program.signature) <- paste0("T2_Program", 1:length(T2.program.signature))
names(T3.program.signature) <- paste0("T3_Program", 1:length(T3.program.signature))
names(T4.program.signature) <- paste0("T4_Program", 1:length(T4.program.signature))
All.program.signature <- c(T1.program.signature, T2.program.signature, T3.program.signature, T4.program.signature)
saveRDS(All.program.signature, file = "All.program.signature.rds")
#out put list
a <- as.data.frame(All.program.signature)
write.csv(a, file = "All.program.signature.csv", row.names = F)

####Calculate the program score in the overall scale.data data and then perform correlation comparison analysis
CancerCell.four.matrix <- GetAssayData(CancerCell, slot = "scale.data")

##1. Calculate the program score of each sample across samples
CancerCell <- AddModuleScore(CancerCell, features = All.program.signature, name = paste0(names(All.program.signature), "_"))
program.score <- CancerCell@meta.data[, grep("Program", colnames(CancerCell@meta.data))]
colnames(program.score) <- gsub("_\\d+", "", colnames(program.score))
saveRDS(program.score, file = "program.score.rds")

##2. Calculation of each program score of a single sample
T1.cancer.data <- AddModuleScore(T1.cancer.data, features = All.program.signature, name = paste0(names(All.program.signature), "_"))
T1.program.score <- T1.cancer.data@meta.data[, grep("Program", colnames(T1.cancer.data@meta.data))]
colnames(T1.program.score) <- gsub("_\\d+", "", colnames(T1.program.score))
T2.cancer.data <- AddModuleScore(T2.cancer.data, features = All.program.signature, name = paste0(names(All.program.signature), "_"))
T2.program.score <- T2.cancer.data@meta.data[, grep("Program", colnames(T2.cancer.data@meta.data))]
colnames(T2.program.score) <- gsub("_\\d+", "", colnames(T2.program.score))
T3.cancer.data <- AddModuleScore(T3.cancer.data, features = All.program.signature, name = paste0(names(All.program.signature), "_"))
T3.program.score <- T3.cancer.data@meta.data[, grep("Program", colnames(T3.cancer.data@meta.data))]
colnames(T3.program.score) <- gsub("_\\d+", "", colnames(T3.program.score))
T4.cancer.data <- AddModuleScore(T4.cancer.data, features = All.program.signature, name = paste0(names(All.program.signature), "_"))
T4.program.score <- T4.cancer.data@meta.data[, grep("Program", colnames(T4.cancer.data@meta.data))]
colnames(T4.program.score) <- gsub("_\\d+", "", colnames(T4.program.score))
single.programScore <- list(T1 = T1.program.score, T2 = T2.program.score, T3 = T3.program.score, T4 = T4.program.score)
saveRDS(single.programScore, file = "single.programScore.rds")

# program score correlation analysis
pdf("program.correlation.pdf")
corrM <- cor(program.score)
corrplot(corrM, order = "hclust", col=colorRampPalette(c("blue","white","#EA4335"))(100), hclust.method = "ward.D", title = "All sample", addCoef.col = "black", number.cex = 0.75, number.font = 3, tl.col = "black", tl.cex = 0.8)
T1.corrM <- cor(T1.program.score)
corrplot(T1.corrM, order = "hclust", col=colorRampPalette(c("blue","white","#EA4335"))(100), hclust.method = "ward.D", title = "T1", addCoef.col = "black", number.cex = 0.75, number.font = 3, tl.col = "black", tl.cex = 0.8)
T2.corrM <- cor(T2.program.score)
corrplot(T2.corrM, order = "hclust", col=colorRampPalette(c("blue","white","#EA4335"))(100), hclust.method = "ward.D", title = "T2", addCoef.col = "black", number.cex = 0.75, number.font = 3, tl.col = "black", tl.cex = 0.8)
T3.corrM <- cor(T3.program.score)
corrplot(T3.corrM, order = "hclust", col=colorRampPalette(c("blue","white","#EA4335"))(100), hclust.method = "ward.D", title = "T3", addCoef.col = "black", number.cex = 0.75, number.font = 3, tl.col = "black", tl.cex = 0.8)
T4.corrM <- cor(T4.program.score)
corrplot(T4.corrM, order = "hclust", col=colorRampPalette(c("blue","white","#EA4335"))(100), hclust.method = "ward.D", title = "T4", addCoef.col = "black", number.cex = 0.75, number.font = 3, tl.col = "black", tl.cex = 0.8)
mean.corrM <- (T1.corrM+T2.corrM+T3.corrM+T4.corrM)/4
corrplot(mean.corrM, order = "hclust", col=colorRampPalette(c("blue","white","#EA4335"))(100), hclust.method = "ward.D", title = "Mean correlation", addCoef.col = "black", number.cex = 0.75, number.font = 3, tl.col = "black", tl.cex = 0.8)
dev.off()

#### Functional enrichment analysis based on these programs
source("/home/longzhilin/Analysis_Code/PathwayEnrichment/clusterProfiler.enricher.R")
pdf("program.pathway.enrichment.pdf")
res <- sapply(names(All.program.signature), function(x){
    y <- All.program.signature[[x]]
    res <- cluterProfiler.enricher(gene = y, geneType = "SYMBOL", db.type = "MsigDB", saveDir = paste0(getwd(),"/cluterProfiler_MsigDB_Program"),
                            title = x, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
    return(res)
})
dev.off()

############2.2 Calculate the similarity between cells
####Cell similarity
#First realize clustering through ComplexHeatmap
#Hierarchical clustering, with one minus Pearson correlation as the distance metric and Ward’s linkage
#The first way: first cluster a single sample in the overall scale matrix, obtain the sort position, and then combine
cluster.order <- sapply(c("T1_", "T2_", "T3_", "T4_"), function(x){
    index <- grep(x, colnames(CancerCell.four.matrix)) 
    patient.matrix <- CancerCell.four.matrix[,index]
    corrMatrix <- (1- cor(patient.matrix, method="pearson"))
    hc <- hclust(as.dist(corrMatrix), method="ward.D")
    row_dend <- dendsort(hc) 
    return(colnames(patient.matrix)[row_dend$order])
})
cluster.order <- unlist(cluster.order)
index <- match(cluster.order, colnames(CancerCell.four.matrix))
CancerCell.four.matrix <- CancerCell.four.matrix[,index]
saveRDS(CancerCell.four.matrix, file = "CancerCell.four.matrix.rds")
cell.similarity <- cor(CancerCell.four.matrix, method="pearson")
saveRDS(cell.similarity, file = "cell.similarity.rds")
diag(cell.similarity) <- 0
diag(cell.similarity) <- max(cell.similarity)
pdf("cell.similarity.pdf")
row_split <- gsub("_\\d+", "", names(cluster.order))
column_split <- gsub("_\\d+", "", names(cluster.order))
Heatmap(cell.similarity, cluster_rows = F, border = TRUE, row_split = row_split, column_split = column_split, row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), width = unit(12, "cm"), height = unit(12, "cm"), cluster_columns = F, show_column_names = F, show_row_names = F, use_raster = T, heatmap_legend_param = list(title = "Correlation coeficient"))
dev.off()

#############2.3 Observe the similarity of each program score based on the order of cell similarity
pdf("programScore.hierarchical.clustering.pdf")
## all patient
a <- t(program.score)
index <- match(cluster.order, colnames(a))
a <- a[,index]
p1 <- Heatmap(a, cluster_rows = T, border = TRUE, row_split = 4, clustering_distance_rows = "pearson", clustering_method_rows = "ward.D", column_split = column_split, row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), cluster_columns = F, show_column_names = F, show_row_names = F, use_raster = T, width = unit(8, "cm"), height = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), heatmap_legend_param = list(title = "Expression score"))

program.matrix <- matrix(0, ncol = 4, nrow = nrow(a))
colnames(program.matrix) <- c("T1", "T2", "T3", "T4")
rownames(program.matrix) <- rownames(a)
for(i in 1:nrow(a)){
    x <- rownames(a)[i]
    x <- gsub("_.*", "", x)
    program.matrix[i, x] <- 1
}
colors = structure(c("white", "#6387C5"), names = c("0", "1"))
p2 <- Heatmap(program.matrix, col = colors, cluster_rows = F, rect_gp = gpar(col = "white", lwd = 2), border = TRUE, row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), cluster_columns = F, show_column_names = T, show_row_names = F, width = unit(1.5, "cm"), row_names_gp = gpar(fontsize = 10))
plot.list <- p1+p2
draw(plot.list, heatmap_height = unit(6, "cm"))

p1 <- Heatmap(a, cluster_rows = T, border = TRUE, row_split = 4, clustering_distance_rows = "pearson", clustering_method_rows = "ward.D", column_split = column_split, row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), cluster_columns = F, show_column_names = F, show_row_names = T, use_raster = T, width = unit(8, "cm"), height = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), heatmap_legend_param = list(title = "Expression score"))
p2 <- Heatmap(program.matrix, col = colors, cluster_rows = F, rect_gp = gpar(col = "white", lwd = 2), border = TRUE, row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), cluster_columns = F, show_column_names = T, show_row_names = T, width = unit(1.5, "cm"), row_names_gp = gpar(fontsize = 10))
plot.list <- p1+p2
draw(plot.list, heatmap_height = unit(6, "cm"))

## single patient
singleScore.cluster <- function(program.score, cluster.order, row_split){
    a <- t(program.score)
    index <- match(cluster.order, colnames(a))
    a <- a[,na.omit(index)]
    p1 <- Heatmap(a, cluster_rows = T, border = TRUE, row_split = row_split, clustering_distance_rows = "pearson", clustering_method_rows = "ward.D", cluster_column = F, row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), cluster_columns = F, show_column_names = F, show_row_names = T, use_raster = T, width = unit(8, "cm"), height = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), heatmap_legend_param = list(title = "Expression score"))
    print(p1)
}
res <- singleScore.cluster(program.score = T1.program.score, cluster.order = cluster.order, row_split = 5)
res <- singleScore.cluster(program.score = T2.program.score, cluster.order = cluster.order, row_split = 5)
res <- singleScore.cluster(program.score = T3.program.score, cluster.order = cluster.order, row_split = 5)
res <- singleScore.cluster(program.score = T4.program.score, cluster.order = cluster.order, row_split = 5)
Heatmap(mean.corrM, cluster_rows = T, cluster_columns = T, border = TRUE, row_split = 5, column_split = 5, row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), show_column_names = T, show_row_names = T, use_raster = T, width = unit(6, "cm"), height = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), heatmap_legend_param = list(title = "Mean correlation"))
Heatmap(mean.corrM, cluster_rows = T, cluster_columns = T, clustering_distance_rows = "pearson", clustering_method_rows = "ward.D", clustering_distance_columns = "pearson", clustering_method_columns = "ward.D", border = TRUE, row_split = 5, column_split = 5, row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), show_column_names = T, show_row_names = T, use_raster = T, width = unit(6, "cm"), height = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), heatmap_legend_param = list(title = "Mean correlation"))
dev.off()

###############2.4 Build metaProgram
#ref:Resolving medulloblastoma cellular architecture by single-cell genomics. 2019 Nature
All.program.signature <- readRDS("All.program.signature.rds")
metaProgram1 <- c("T2_Program2", "T1_Program1", "T3_Program1", "T4_Program1")
metaProgram2 <- c("T4_Program2", "T2_Program1", "T3_Program3")
#metaProgram3 <- c("T2_Program3", "T3_Program2")

#' @note function ---  metaProgram signature gene
#average NMF score top30 gene as the signature
getMetaSignature <- function(NMF.score, program.set){
    # first step：extract NMF score with correlated program
    patient.info <- gsub("_.*", "", program.set)
    program.nmf.score <- lapply(1:length(patient.info), function(x){
        nmf.score <- NMF.score[[patient.info[x]]] # focus patient
        nmf.score <- nmf.score[, program.set[x]] # focus program
        nmf.score <- data.frame(gene = names(nmf.score), score = as.numeric(nmf.score))
        return(nmf.score)
    })
    program.nmf.score.merge <- Reduce(function(x,y) merge(x, y, by="gene", all.x=TRUE, all.y=TRUE), program.nmf.score)
    colnames(program.nmf.score.merge) <- c("gene", program.set)

    # second step: calculate average score
    avg.score <- apply(program.nmf.score.merge, 1, function(x){
        return(sum(na.omit(as.numeric(x[-1])))/length(program.set))
    })
    program.nmf.score.merge$avgScore <- avg.score
    program.nmf.score.merge <- program.nmf.score.merge[order(program.nmf.score.merge$avgScore, decreasing = T), ]
    return(program.nmf.score.merge$gene[1:30]) 
}
meta.signature1 <- getMetaSignature(program.set = metaProgram1, NMF.score = NMF.score)
meta.signature2 <- getMetaSignature(program.set = metaProgram2, NMF.score = NMF.score)
meta.Signature <- list(metaProgram1 = meta.signature1, 
                       metaProgram2 = meta.signature2)
a <- as.data.frame(meta.Signature)
write.csv(a, file = "meta.Signature.csv", row.names = F)
saveRDS(meta.Signature, file = "meta.Signature.rds")

##Functional enrichment analysis
source("/home/longzhilin/Analysis_Code/PathwayEnrichment/clusterProfiler.enricher.R")
pdf("MetaProgram.pathway.enrichment.pdf")
res <- sapply(1:length(meta.Signature), function(x){
    y <- meta.Signature[[x]]
    res <- cluterProfiler.enricher(gene = y, geneType = "SYMBOL", db.type = "MsigDB", saveDir = paste0(getwd(), "/cluterProfiler_MsigDB_MetaProgram"),
                            title = x, qvalueCutoff = 0.05, pvalueCutoff = 0.05, showCategory = 10)
    return(res)
})
dev.off()

Expression.matrix <- CancerCell.four.matrix
meta.genes <- unlist(meta.Signature[c(1:length(meta.Signature))]) #meta gene signature
exp.matrix <- Expression.matrix[meta.genes,]
row_split <- rep(names(meta.Signature)[c(1:length(meta.Signature))],each=30)
column_split <- gsub("_.+", "", colnames(CancerCell.four.matrix))
pdf("MetaSignature.Expression.pdf")
Heatmap(exp.matrix, cluster_rows = T, cluster_columns = T, row_split = row_split, column_split = column_split, width = 12, height = 10, show_column_dend = F, show_row_dend = F, show_column_names = F, show_row_names = F, use_raster = T, row_names_gp = gpar(fontsize = 10), heatmap_legend_param = list(title = "Relative expression"))   
Heatmap(exp.matrix, cluster_rows = T, cluster_columns = T, row_split = row_split, column_split = column_split, width = 12, height = 10, show_column_dend = F, show_row_dend = F, show_column_names = F, show_row_names = T, use_raster = T, row_names_gp = gpar(fontsize = 10), heatmap_legend_param = list(title = "Relative expression"))   
dev.off()

##Analyze the correlation between the two meatPrograms
#Based on the original strategy for programScore, --- mainly consider this,
CancerCell <- AddModuleScore(CancerCell, features = meta.Signature, name =paste0(names(meta.Signature), "_"))
metaSignature.score <- CancerCell@meta.data[, grep("metaProgram", colnames(CancerCell@meta.data))]
colnames(metaSignature.score) <- gsub("_\\d+", "", colnames(metaSignature.score))
metaSignature.score$orig.ident <- CancerCell@meta.data$orig.ident
saveRDS(metaSignature.score, file = "metaSignature.score.rds")
saveRDS(CancerCell, file = "CancerCell.pro.rds")

cell.group <- apply(metaSignature.score[,1:ncol(metaSignature.score)], 1, function(x){
    idx <- which(x == max(x))
    return(idx)
})
cell.group <- gsub(1, "metaProgram1", cell.group)
cell.group <- gsub(2, "metaProgram2", cell.group)
CancerCell$program.group <- cell.group
pdf("program.group.pdf")
DimPlot(CancerCell, group.by = "program.group")
ggscatter(metaSignature.score, x = "metaProgram1", y = "metaProgram2", add = "reg.line", cor.coef = T)
ggscatter(metaSignature.score, x = "metaProgram1", y = "metaProgram2", color = "orig.ident", add = "reg.line", cor.coef = T)
dev.off()

#############2.6 Functional enrichment analysis based on gprofiler2
#hypergeometric test followed by correction for multiple testing
source(file = "/home/longzhilin/Analysis_Code/PathwayEnrichment/pathwayEnrichment.gProfiler.R")
gProfiler.res <- pathwayEnrichment.gProfiler(geneList = meta.Signature, saveDir = paste0(getwd(), "/gProfiler"),
                            filename = "MetaProgram", ordered_query = T, max.term.size = 350, user_threshold = 0.001, dbs = c("GO:BP", "REAC", "KEGG"), log10fdr = F)

####NMF enrich the TF database using 
source(file = "/home/longzhilin/Analysis_Code/PathwayEnrichment/pathwayEnrichment.Enrichr.R")
pdf(paste0(getwd(), "/enrichr/enrichr.result.pdf"))
enrichr.res <- pathwayEnrichment.Enrichr(geneList = meta.Signature, 
                                         dbs = c("ChEA_2016", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "TF_Perturbations_Followed_by_Expression", "ARCHS4_TFs_Coexp"))
dev.off()
write.xlsx(enrichr.res$metaProgram1, file = paste0(getwd(), "/enrichr/metaProgram1.xlsx"), rowNames = F)
write.xlsx(enrichr.res$metaProgram2, file = paste0(getwd(), "/enrichr/metaProgram2.xlsx"), rowNames = F)

#############2.8 TCGA survival
source(file = "/home/longzhilin/Analysis_Code/code/analysis.diff.survival.TCGA.R")
DESeq2.normalized_counts <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/DESeq2.normalized_counts.rds")
DESeq2.normalized_counts <- log2(DESeq2.normalized_counts+1)
DESeq2.result <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/DESeq2.result.rds")
clin.data <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/clin.data.rds")
pdf("meta.signature.survival.pdf")
Metabolic.TCGA <- analysis.diff.survival.TCGA(interest.gene = meta.Signature$metaProgram1, diff.gene.pro = DESeq2.result, exp.data.process = DESeq2.normalized_counts, clin.data = clin.data, EnhancedVolcano.plot = F, main = "Meta program1", Box.plot = F, meta.signature = T, single.signature = F)
Immune.TCGA <- analysis.diff.survival.TCGA(interest.gene = meta.Signature$metaProgram2, diff.gene.pro = DESeq2.result, exp.data.process = DESeq2.normalized_counts, clin.data = clin.data, EnhancedVolcano.plot = F, main = "Meta program2", Box.plot = F, meta.signature = T, single.signature = F)
dev.off()

#############2.9 ICB survival
library(ggpubr)
source(file = "/home/longzhilin/Analysis_Code/SurvivalAnalysis/Cox.function.R")
source(file = "/home/longzhilin/Analysis_Code/code/RCC.ICB.analysis.R")
normalized_expression <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/ICB.therapy/normalized_expression.rds")
patient.info.RNA <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/ICB.therapy/patient.info.RNA.rds")
patient.info.RNA$Sex <- gsub("^F|FEMALE|Female$", 0, patient.info.RNA$Sex)
patient.info.RNA$Sex <- gsub("^M|Male|MALE$", 1, patient.info.RNA$Sex)
patient.info.RNA$Tumor_Sample_Primary_or_Metastasis <- gsub("PRIMARY", 0, patient.info.RNA$Tumor_Sample_Primary_or_Metastasis)
patient.info.RNA$Tumor_Sample_Primary_or_Metastasis <- gsub("METASTASIS", 1, patient.info.RNA$Tumor_Sample_Primary_or_Metastasis)

pdf("meta.signature.ICB.survival.pdf")
cox.res <- RCC.icb.analysis(signature.list = meta.Signature, expresionMatrix = normalized_expression, clincal.info = patient.info.RNA)
dev.off()

####cluster profiler enrichment result display
library(ggthemes)
#top10
enrich1 <- read.csv("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/4.NMF/SCT/SD/cluterProfiler_MsigDB_MetaProgram/1_enricherResult.csv", header = T)
enrich2 <- read.csv("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/4.NMF/SCT/SD/cluterProfiler_MsigDB_MetaProgram/2_enricherResult.csv", header = T)
enrich <- list(MetaProgram1 = enrich1[c(1:3,7,9,12,14,20,21,27,33,77,131,141,165), c("ID", "Count", "p.adjust")],
               MetaProgram2 = enrich2[c(1:4,6,7,11,15,17,18,26,33,39,42,43), c("ID", "Count", "p.adjust")])
enrich.res <- Reduce(function(x,y) rbind(x,y), enrich)
enrich.res$Type <- rep(c("Meta program1", "Meta program2"), each = 15)
enrich.res$ID <- factor(enrich.res$ID, levels = unique(enrich.res$ID))
enrich.res$p.adjust <- -log10(enrich.res$p.adjust)
#处理名字
enrich.res$ID <- gsub("_", " ", enrich.res$ID)
enrich.res$ID <- factor(enrich.res$ID, levels = unique(enrich.res$ID))
#source(file = "/home/longzhilin/Analysis_Code/Visualization/ggplot.bubble.plot.R")
pdf("ggplot.bubble.pathway.pdf", width = 6.6)
p1 <- ggplot(enrich.res, 
            aes(x = Type, 
                y = ID, 
                size = Count, 
                fill = p.adjust))
p2 <- p1 + guides(color=FALSE) + geom_point(shape = 21, alpha = 0.7) + theme_few() + scale_fill_gradientn(colours = c("#fed71b", "#fcc41c", "#f27620", "#e92925"), breaks = c(1.39794, 2, 4, 6), labels = c(0.04, 0.01, 0.0001, 0.000001))
p2 <- p2 + theme(axis.text.x = element_text(angle = 90,vjust = 1,hjust = 1, size = 8), axis.text.y = element_text(size = 8)) + xlab("") + ylab("")
print(p2)
dev.off()

values = c(1, 1.3, 8), breaks = c(1, 1.3, 8), labels = c(0.1, 0.05, 0.00000001)
enri ch <- list(MetaProgram1 = enrich1[c(1:3,7,9,12,14,20,21,27,33,77,131,141,165), c("ID", "Count", "p.adjust", "geneID")], #
               MetaProgram2 = enrich2[c(1:4,6,7,11,15,17,18,26,33,39,42,43), c("ID", "Count", "p.adjust", "geneID")])
pathway.gene <- sapply(enrich, function(x){
    res <- apply(x, 1, function(y){
        a <- unlist(strsplit(y[4], split = "/"))
        return(a)
    })
    return(unlist(res))
})
program1 <- as.data.frame(table(pathway.gene$MetaProgram1))
program2 <- as.data.frame(table(pathway.gene$MetaProgram2))
pdf("pathway.hub.gene.pdf")
ggbarplot(program1, x = "Var1", y = "Freq", sort.val = c("desc")) + xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggbarplot(program2, x = "Var1", y = "Freq", sort.val = c("desc")) + xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()