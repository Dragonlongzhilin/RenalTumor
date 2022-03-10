setwd("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result")
source(file = "/data/active_data/lzl/RenalTumor-20200713/code/analysis.diff.survival.TCGA.R")

#### surivial analysis in TCGA
#http://firebrowse.org/?cohort=KIRC&download_dialog=true，数据下载来源
exp.data <- read.table("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/KIRC__RSEM_normalized__data.data.txt", 
                        header = T, stringsAsFactors = F, sep = "\t")
exp.data <- exp.data[-1,]
rownames(exp.data) <- exp.data[,1]
exp.data.process <- exp.data[,-1]
exp.data.process <- apply(exp.data.process, 2, as.numeric)
rownames(exp.data.process) <- rownames(exp.data)
exp.data.process <- log2(exp.data.process+1)
saveRDS(exp.data.process, file = "Log2.expression.rds") #20531*606

#### raw count
exp.matrix <- read.table("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/KIRC_RSEM_genes_data.data.txt", 
                        header = T, stringsAsFactors = F, sep = "\t")
count.index <- seq(1,1818,by=3)+1
exp.matrix.count <- exp.matrix[,count.index]
exp.matrix.count <- exp.matrix.count[-1,]
rownames(exp.matrix.count) <- exp.matrix[2:nrow(exp.matrix),1]

#处理表达谱名字
gene.id <- gsub("\\|\\d+", "", rownames(exp.matrix.count))
index <- which(gene.id == "?")
exp.data.process.pro <- exp.data.process[-index,]
exp.matrix.count.pro <- exp.matrix.count[-index,]
gene.id <- gene.id[-index]
#去除重复--- SLC35E2
gene.id <- gene.id[-16273]
exp.data.process.pro <- exp.data.process.pro[-16273,]
exp.matrix.count.pro <- exp.matrix.count.pro[-16273,]
rownames(exp.data.process.pro) <- gene.id
rownames(exp.matrix.count.pro) <- gene.id
saveRDS(exp.data.process.pro, file = "Log2.expression.pro.rds") #20501*606
saveRDS(exp.matrix.count.pro, file = "exp.matrix.count.pro.rds") #20501*606

#### process the survival data
clin.data <- read.table("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/KIRC_clin.txt", header = T, stringsAsFactors = F, sep = "\t")
clin.data$Sample <- gsub("-", ".", clin.data$Sample)
index <- grep("01$", clin.data$Sample)
clin.data <- clin.data[index, ]
saveRDS(clin.data, file = "clin.data.rds")

##############################1.Differential expression analysis#########################################
# Distinguish between normal and cancer samples
sample.id <- colnames(exp.data.process.pro)
sample.id <- substr(sample.id, 1, 16)
sample.type <- substr(sample.id, 14, 15)

##normal
normal.index <- which(sample.type == "11") #72
cancer.index <- which(sample.type == "01") #533

####DEseq2 analysis
library(DESeq2)
exp.matrix.count.pro2 <- as.matrix(exp.matrix.count.pro)
exp.matrix.count.pro2 <- apply(exp.matrix.count.pro2, 2, as.numeric)
exp.matrix.count.pro2 <- round(exp.matrix.count.pro2)
index <- which(rowSums(exp.matrix.count.pro2)>0)
exp.matrix.count.pro2 <- exp.matrix.count.pro2[index,] #去除不表达的基因，
rownames(exp.matrix.count.pro2) <- rownames(exp.matrix.count.pro)[index] #20221*606
exp.matrix.count.pro2 <- exp.matrix.count.pro2[,c(normal.index, cancer.index)] #605个样本
condition <- rep("Tumor", ncol(exp.matrix.count.pro2))
condition[1:length(normal.index)] <- "Normal"
coldata <- data.frame(row.names = colnames(exp.matrix.count.pro2), condition = as.factor(condition))
dds <- DESeqDataSetFromMatrix(countData = exp.matrix.count.pro2, colData = coldata, design=~condition)
dds <- DESeq(dds)
DESeq2.result <- results(dds)
saveRDS(DESeq2.result, file = "DESeq2.result.rds")
#提取标化后的counts
DESeq2.normalized_counts <- counts(dds, normalized=TRUE)
saveRDS(DESeq2.normalized_counts, file = "DESeq2.normalized_counts.rds")