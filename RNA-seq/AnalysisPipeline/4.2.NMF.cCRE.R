#' @description: Analyze the TFs involved in each NMF program
set.seed(101)
library(Seurat)
library(Signac)
library(openxlsx)
library(dplyr)
library(tibble)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)

setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/4.NMF/SCT/SD")
meta.Signature <- readRDS("meta.Signature.rds")
cCREs <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/6.Co-Accessible/Tumor_peak_gene_correlation.rds")

source(file = "/home/longzhilin/Analysis_Code/code/analysis.diff.survival.TCGA.R")
DESeq2.normalized_counts <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/DESeq2.normalized_counts.rds")
DESeq2.normalized_counts <- log2(DESeq2.normalized_counts+1)
DESeq2.result <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/DESeq2.result.rds")
clin.data <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/clin.data.rds")

scATAC.data <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/scATAC.data.pro.rds")
scRNA.data <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/data.merge.pro.rds")

####---------------------------------------------------1.Observe the expression of metaSiganture in scATAC data
DefaultAssay(scATAC.data) <- "ACTIVITY"
Tumor.scATAC <- subset(scATAC.data, subset = AnnotatedcellType == "Tumor")
exp.matrix <- GetAssayData(Tumor.scATAC, slot = "scale.data")
meta.genes <- data.frame(genes = unlist(meta.Signature[c(1:length(meta.Signature))]), type = c(rep("Meta-program1", 30), rep("Meta-program2", 30)))
rownames(meta.genes) <- NULL
idx <- match(meta.genes$genes, rownames(exp.matrix))
meta.exp.matrix <- exp.matrix[na.omit(idx),]
meta.genes <- meta.genes[which(!is.na(idx)),]
row_split <- meta.genes$type
column_split <- gsub("_.+", "", colnames(exp.matrix))
pdf("TF_cCREs/MetaSignature.Expression.pdf")
Heatmap(as.matrix(meta.exp.matrix), cluster_rows = T, cluster_columns = T, row_split = row_split, column_split = column_split, width = 12, height = 10, show_column_dend = F, show_row_dend = F, show_column_names = F, show_row_names = F, use_raster = T, row_names_gp = gpar(fontsize = 10), heatmap_legend_param = list(title = "Gene activity"))   
Heatmap(as.matrix(meta.exp.matrix), cluster_rows = T, cluster_columns = T, row_split = row_split, column_split = column_split, width = 12, height = 10, show_column_dend = F, show_row_dend = F, show_column_names = F, show_row_names = T, use_raster = T, row_names_gp = gpar(fontsize = 10), heatmap_legend_param = list(title = "Gene activity"))   
dev.off()

####---------------------------------------------------2.The intersection of metaSignature and DEG, DAR
DEGs <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/2.Cluster/AnnotateCellType/cellType.all.DEGs.rds")
Tumor.DEGs <- DEGs %>% filter(cluster == "Tumor") %>% filter(abs(avg_log2FC)>0.25 & p_val_adj < 0.05)
NMF.DEGs <- lapply(meta.Signature, function(x){
    a <- sapply(x, function(y){
        idx <- which(Tumor.DEGs$gene == y)
        if(length(idx)>0){
            return(Tumor.DEGs$avg_log2FC[idx])
        }else{
            return(NA)
        }
    })
    return(a)
})

DARs <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/4.Peak/cellType.all.DARs.rds")
Tumor.DARs <- DARs %>% filter(cellType == "Tumor") %>% filter(abs(avg_log2FC)>0.25 & p_val_adj < 0.05)
NMF.DARs <- lapply(meta.Signature, function(x){
    a <- sapply(x, function(y){
        idx <- which(Tumor.DARs$gene == y)
        if(length(idx)>0){
            return(max(Tumor.DARs$avg_log2FC[idx]))
        }else{
            return(NA)
        }
    })
    return(a)
})
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
pdf("TF_cCREs/MetaSignature.DEG.DAR.pdf")
a <- matrix(NMF.DEGs$metaProgram1)
rownames(a) <- names(NMF.DEGs$metaProgram1)
p1 <- Heatmap(a, name = "log2FC", col = col_fun, width = unit(1, "cm"), height = unit(8, "cm"), na_col = "grey", cluster_rows = F, cluster_columns = F,
        show_column_dend = F, show_row_dend = F, show_column_names = F, show_row_names = T, row_names_side = "left",
        row_names_gp = gpar(fontsize = 8))
a <- matrix(NMF.DARs$metaProgram1)
rownames(a) <- names(NMF.DARs$metaProgram1)
p2 <- Heatmap(a, name = "log2FC", col = col_fun, width = unit(1, "cm"), height = unit(8, "cm"), na_col = "grey", cluster_rows = F, cluster_columns = F,
        show_column_dend = F, show_row_dend = F, show_column_names = F, show_row_names = T, row_names_side = "left",
        row_names_gp = gpar(fontsize = 8))
p1+p2
a <- matrix(NMF.DEGs$metaProgram2)
rownames(a) <- names(NMF.DEGs$metaProgram2)
p1 <- Heatmap(a, name = "log2FC", col = col_fun, width = unit(1, "cm"), height = unit(8, "cm"), na_col = "grey", cluster_rows = F, cluster_columns = F,
        show_column_dend = F, show_row_dend = F, show_column_names = F, show_row_names = T, row_names_side = "left",
        row_names_gp = gpar(fontsize = 8))
a <- matrix(NMF.DARs$metaProgram2)
rownames(a) <- names(NMF.DARs$metaProgram2)
p2 <- Heatmap(a, name = "log2FC", col = col_fun, width = unit(1, "cm"), height = unit(8, "cm"), na_col = "grey", cluster_rows = F, cluster_columns = F,
        show_column_dend = F, show_row_dend = F, show_column_names = F, show_row_names = T, row_names_side = "left",
        row_names_gp = gpar(fontsize = 8))
p1+p2
dev.off()

####---------------------------------------------------3. peak in NMF program gene promoter
# cCREs in NMF program
meta.program.cCREs <- lapply(meta.Signature, function(x){
    idx <- match(x, cCREs$peak1_nearestGene)
    return(cCREs[na.omit(idx),])
})
write.xlsx(meta.program.cCREs, file = "TF_cCREs/meta.program.cCREs.xlsx", sheetName = names(meta.program.cCREs), rowNames = F)

# TF target meta-program genes
All.TF.targets <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/7.TF.analysis/targets/Tumor.All.TF.targets.rds")
TF.metaProgram <- lapply(meta.Signature, function(x){
    idx <- which(All.TF.targets$to %in% x)
    return(All.TF.targets[na.omit(idx),])
})

pdf("TF_cCREs/TF.regulation.Type.pdf", height = 4)
TF.regulate.number <- lapply(names(TF.metaProgram), function(y){
    x <- TF.metaProgram[[y]]
    a <- as.data.frame(table(x[, c("from","type")]), stringsAsFactors = F)
    #取top30
    regulation.num <- tapply(a$Freq, a$from, sum)
    regulation.num <- sort(regulation.num, decreasing = T)
    top30 <- names(regulation.num)[1:30]
    top.data <- a[which(a$from %in% top30),]
    top.data$from <- factor(top.data$from, levels = top30)
    write.xlsx(arrange(top.data, from), file = paste0("TF_cCREs/TF.", y, ".top30.xlsx"), rowNames = F)
    p <- ggbarplot(top.data, x = "from", y = "Freq", fill = "type", color = "type", lab.pos = "in",
                   label = T, xlab = "", ylab = paste0("Target numbers in", y)) + theme(legend.position="top", axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)

    # calculate the regulation number of each TFs
    targetNum <- table(unique(x[, c("from","to")])[,1])
    targetNum <- sort(targetNum, decreasing = T)
    targetNum <- data.frame(TF = names(targetNum), targetNum = as.numeric(targetNum))
    return(targetNum)
})
dev.off()
names(TF.regulate.number) <- names(TF.metaProgram)
# regulation number
pdf("TF_cCREs/TF.regulation.number.pdf" , height = 3)
p <- ggbarplot(TF.regulate.number$metaProgram1[1:30,], x = "TF", y = "targetNum", fill = "grey", color = "grey", lab.pos = "out", sort.val = "desc",
                label = T, xlab = "", ylab = "Target numbers in meta-program1") + theme(legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1))
print(p)
p <- ggbarplot(TF.regulate.number$metaProgram2[1:30,], x = "TF", y = "targetNum", fill = "grey", color = "grey", lab.pos = "out", sort.val = "desc",
                label = T, xlab = "", ylab = "Target numbers in meta-program2") + theme(legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1))
print(p)
dev.off()

### Screen TF-specific regulation of a program
program1 <- TF.regulate.number$metaProgram1[which(TF.regulate.number$metaProgram1$targetNum>=10),]
program2 <- TF.regulate.number$metaProgram2[which(TF.regulate.number$metaProgram2$targetNum>=10),]
spe.TF1 <- setdiff(program1$TF, program2$TF)
spe.TF2 <- setdiff(program2$TF, program1$TF)
spe.TFs <- c(spe.TF1, spe.TF2)
data1 <- TF.regulate.number$metaProgram1[match(spe.TFs, TF.regulate.number$metaProgram1$TF),]
data2 <- TF.regulate.number$metaProgram2[match(spe.TFs, TF.regulate.number$metaProgram2$TF),]
TF.data <- cbind(data1, data2[,2])
colnames(TF.data) <- c("TF", "Meta_program1", "Meta_program2")
TF.data1 <- TF.data %>% filter(Meta_program1>=10 & Meta_program2<5)
TF.data2 <- TF.data %>% filter(Meta_program1<5 & Meta_program2>10)

####### screen the TF associated NMF program
#### 1.combined the top 30 TF in each meta-program
TF.regulations <- Reduce(function(x,y) merge(x, y, by="TF", all.x=TRUE), TF.regulate.number)
colnames(TF.regulations)[2:ncol(TF.regulations)] <- names(TF.metaProgram)
rownames(TF.regulations) <- TF.regulations$TF
top.list <- sapply(TF.regulate.number, function(x){return(x$TF[1:30])})
top.list <- unique(as.character(top.list))
TF.regulation.data <- TF.regulations[top.list,]

#### 2.combined the tumor enriched TF information
Tumor.enriched.TFs <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/5.Motif/Analysis/tumor.specific.TFs.rds")
chromVAR.TFs <- TF.regulations[Tumor.enriched.TFs$Name,]
chromVAR.TFs <- chromVAR.TFs[which(!is.na(chromVAR.TFs$TF)),]
pdf("TF_cCREs/TF.cCREs.chromVAR.pdf")
col_fun <- colorRamp2(c(0, 15, 30), c("#00B3B6", "white", "#FF2610"))
a <- TF.regulation.data[,2:ncol(TF.regulation.data)]
p <- Heatmap(a, show_row_dend = F, show_column_dend = F, col = col_fun,
             width = unit(2, "cm"), name = "Target number",
             row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
             cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(a[i, j], x, y, gp = gpar(fontsize = 7))})
print(p)

col_fun <- colorRamp2(c(0, 7, 14), c("#00B3B6", "white", "#FF2610"))
a <- chromVAR.TFs[,2:ncol(chromVAR.TFs)]
p <- Heatmap(a, show_row_dend = F, show_column_dend = F, col = col_fun, cluster_rows = F, cluster_columns = F,
             width = unit(2, "cm"), height = unit(14, "cm"), name = "Target number",
             row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
             cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(a[i, j], x, y, gp = gpar(fontsize = 7))})
print(p)
dev.off()

#### 3. TF within DEG in tumor
source(file = "/home/longzhilin/Analysis_Code/Visualization/Plot.EnhancedVolcano.R")
scRNA.DEGs <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/2.Cluster/AnnotateCellType/cellType.all.DEGs.rds")
tumor.DEGs <- scRNA.DEGs %>% dplyr::filter(cluster == "Tumor")
pdf("TF_cCREs/TF.metaProgram.DEGs.pdf")
TF.regulate.DEGs <- lapply(TF.regulate.number, function(x){
    idx <- match(x$TF, tumor.DEGs$gene)
    x$avg_log2FC <- tumor.DEGs$avg_log2FC[idx]
    x$p_val_adj <- tumor.DEGs$p_val_adj[idx]
    return(x) 
})
##Draw a volcano map that regulates more than 5 genes
plot.data <- TF.regulate.DEGs$metaProgram1[TF.regulate.DEGs$metaProgram1$targetNum >=5, ]
plot.data$pointSize <- plot.data$targetNum*0.2 # 384
p <- Plot.EnhancedVolcano(plot.data, x = "avg_log2FC", y = "p_val_adj", id.column = "TF", 
                          FC.cutoff = 0.25, selected.showType = c("FC"), select.num = 10, 
                          drawConnectors = T, pointSize = plot.data$pointSize, title = "Meta-program1 associated with TFs")
p <- p + guides(colour = guide_legend(override.aes = list(size= c(1, 2, 3, 4))))
print(p)

plot.data <- TF.regulate.DEGs$metaProgram2[TF.regulate.DEGs$metaProgram2$targetNum >=5, ]
plot.data$pointSize <- plot.data$targetNum*0.2 # 435
p <- Plot.EnhancedVolcano(plot.data, x = "avg_log2FC", y = "p_val_adj", id.column = "TF", 
                            FC.cutoff = 0.25, selected.showType = c("FC"), select.num = 10, 
                            drawConnectors = T, pointSize = plot.data$pointSize, title = "Meta-program2 associated with TFs")
p <- p + guides(colour = guide_legend(override.aes = list(size= c(1, 2, 3, 4))))
print(p)
write.xlsx(TF.regulate.DEGs, file = "TF_cCREs/TF.regulate.DEGs.xlsx", sheetName = names(TF.regulate.DEGs), rowNames = F)
dev.off()

interest.genes <- c("EGR1", "HSF4", "JUND", "ATF5", "ZBTB7A", "SPI1", "KLF2", "NR0B1", "ZNF263", "SP1")
pdf("TF_cCREs/NMF.TF.survival.pdf")
NMF.TCGA <- analysis.diff.survival.TCGA(interest.gene = interest.genes, diff.gene.pro = DESeq2.result, exp.data.process = DESeq2.normalized_counts, clin.data = clin.data, EnhancedVolcano.plot = F, main = "NMF TFs", Box.plot = F, meta.signature = F, single.signature = T)
dev.off()

####---- plot sankey 
# look TF motif
meta.program1.specific <- setdiff(TF.regulate.number$metaProgram1$TF, TF.regulate.number$metaProgram2$TF)
idx <- match(meta.program1.specific, TF.regulate.number$metaProgram1$TF)
metaProgram1.regulate.number <- TF.regulate.number$metaProgram1[idx,]
meta.program2.specific <- setdiff(TF.regulate.number$metaProgram2$TF, TF.regulate.number$metaProgram1$TF)
idx <- match(meta.program2.specific, TF.regulate.number$metaProgram2$TF)
metaProgram2.regulate.number <- TF.regulate.number$metaProgram2[idx,]
metaProgram1.regulate.number <- filter(metaProgram1.regulate.number, targetNum >= 30*0.1)
metaProgram2.regulate.number <- filter(metaProgram2.regulate.number, targetNum >= 30*0.1)
# plot sankey
program1 <- as.data.frame(table(TF.metaProgram$metaProgram1[,c(1,2)]), stringsAsFactors = F) %>% filter(Freq>0) %>% filter(from %in% meta.program1.specific)
program2 <- as.data.frame(table(TF.metaProgram$metaProgram2[,c(1,2)]), stringsAsFactors = F) %>% filter(Freq>0) %>% filter(from %in% meta.program2.specific)
source("/home/longzhilin/Analysis_Code/Visualization/Plot.Sankey.R")
Plot.Sankey(data = program1, saveDir = file.path(getwd(), "TF_cCREs"), filename = "TF.metaProgram1.specific.sankey")
Plot.Sankey(data = program2, saveDir = file.path(getwd(), "TF_cCREs"), filename = "TF.metaProgram2.specific.sankey")

# top30 TF motifs
program1 <- as.data.frame(table(TF.metaProgram$metaProgram1[,c(1,2)])) %>% filter(Freq>0) %>% filter(from %in% TF.regulate.number$metaProgram1$TF[1:30])
program2 <- as.data.frame(table(TF.metaProgram$metaProgram2[,c(1,2)])) %>% filter(Freq>0) %>% filter(from %in% TF.regulate.number$metaProgram2$TF[1:30])
Plot.Sankey(data = program1, saveDir = file.path(getwd(), "TF_cCREs"), filename = "TF.metaProgram1.top30.sankey")
Plot.Sankey(data = program2, saveDir = file.path(getwd(), "TF_cCREs"), filename = "TF.metaProgram2.top30.sankey")

# OTP, VENTX, ISL1, HOXC5
program1 <- as.data.frame(table(TF.metaProgram$metaProgram1[,c(1,2)])) %>% filter(Freq>0) %>% filter(from %in% c("HOXC5", "ISL1", "VENTX", "OTP"))
program2 <- as.data.frame(table(TF.metaProgram$metaProgram2[,c(1,2)])) %>% filter(Freq>0) %>% filter(from %in% c("HOXC5", "ISL1", "VENTX", "OTP"))
Plot.Sankey(data = program1, saveDir = file.path(getwd(), "TF_cCREs"), filename = "TF.metaProgram1.TF4.sankey", nodeWidth = 15, fontSize = 12, width = 600, height = 300)
Plot.Sankey(data = program2, saveDir = file.path(getwd(), "TF_cCREs"), filename = "TF.metaProgram2.TF4.sankey", nodeWidth = 15, fontSize = 12, width = 600, height = 300)

# cancer specific TFs
program1 <- as.data.frame(table(TF.metaProgram$metaProgram1[,c(1,2)])) %>% filter(Freq>0) %>% filter(from %in% Tumor.enriched.TFs$Name)
program2 <- as.data.frame(table(TF.metaProgram$metaProgram2[,c(1,2)])) %>% filter(Freq>0) %>% filter(from %in% Tumor.enriched.TFs$Name)
Plot.Sankey(data = program1, saveDir = file.path(getwd(), "TF_cCREs"), filename = "TF.metaProgram1.tumorTFs.sankey")
Plot.Sankey(data = program2, saveDir = file.path(getwd(), "TF_cCREs"), filename = "TF.metaProgram2.tumorTFs.sankey")

#### ---- Coverage plot: VEGFA & ALDOB, PCK1
source(file = "/home/longzhilin/Analysis_Code/Plot_colorPaletters.R")
source(file = "/home/longzhilin/Analysis_Code/SingleCell/user.CoveragePlot.R")
df <- readRDS(file = "/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/6.Co-Accessible/Tumor_peak_gene_correlation.rds")
cCREs.genes <- na.omit(unique(c(df$peak1_nearestGene, df$peak2_nearestGene))) # 5407
interest.genes <- c("VEGFA", "EGR1", "ALDOB", "PCK1", "CXCL14", "CRYAB")
NMF.cCREs <- intersect(interest.genes, cCREs.genes)
NMF.cCREs.info <- apply(df, 1, function(x){
  gene1 <- as.character(na.omit(x[7]))
  gene2 <- as.character(na.omit(x[10]))
  idx1 <- which(gene1 %in% NMF.cCREs)
  idx2 <- which(gene2 %in% NMF.cCREs)
  if(length(idx1)>0 | length(idx2)>0){
      return(t(data.frame(c(x))))
  }else{
    return(data.frame())
  }
})
NMF.cCREs.info <- Reduce(rbind, NMF.cCREs.info)
rownames(NMF.cCREs.info) <- NULL
write.table(NMF.cCREs.info, file = "TF_cCREs/Tumor.NMF.cCREs.info.csv", sep = ",", quote=FALSE, row.names=F)
saveRDS(NMF.cCREs.info, file = "TF_cCREs/Tumor.NMF.cCREs.info.rds")

gene_model <- readRDS("/data/ExtraDisk/sdd/longzhilin/Data/ReferenceFasta/Annotation/gene_model.rds")
#CRYAB coverage plot
interest.genes <- c("VEGFA", "EGR1", "ALDOB", "CXCL14", "CRYAB")
Tumor.cicero.conns <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/6.Co-Accessible/ccans/tumor.conns.rds")
Tumor.ccans <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC/6.Co-Accessible/ccans/tumor.ccans.rds")
Tumor.ccans$Peak <- gsub("_", "-", Tumor.ccans$Peak)
rownames(Tumor.ccans) <- Tumor.ccans$Peak
DefaultAssay(scATAC.data) <- "Peaks"
pdf("TF_cCREs/Tumor.NMF.coveragePlot.pdf")
res <- user.CoveragePlot(scATAC.object = scATAC.data,
                        interest.genes = interest.genes, 
                        interested.conns = NMF.cCREs.info, 
                        ccans = Tumor.ccans,
                        gene_model = gene_model,
                        coaccess_cutoff = 0.2,
                        cicero = T,
                        idents = "Tumor",
                        heights = c(2,1,1,2))
dev.off()
pdf("TF_cCREs/Tumor.NMF.gene.pdf")
FeaturePlot(scRNA.data, features = interest.genes, cols = c("lightgrey", "red"), reduction = 'tsne')
FeaturePlot(scRNA.data, features = interest.genes, cols = c("lightgrey", "red"), reduction = 'umap') 
VlnPlot(scRNA.data, features = interest.genes, group.by = "cellType_low", ncol = 2, pt.size = 0)
dev.off()