##Ref: Interplay of somatic alterations and immune infiltration modulates response to PD-1 blockade in advanced clear cell renal cell carcinoma. 2020 Nature Medicine

setwd("/data/active_data/lzl/RenalTumor-20200713/Data/ICB.therapy")

####Observe the value of program signature in ICB treatment data
normalized_expression <- read.table("/data/active_data/lzl/RenalTumor-20200713/Data/ICB.therapy/Checkmate.009.010.025.TPM.matrix.txt", header = T, check.names = F, stringsAsFactors = F, sep = "\t",)
# 43893*312
#Because some gene names have been converted to numbers by excel, they are directly removed
index <- which(normalized_expression[,1] %in% c("1-Mar", "2-Mar", "3-Mar", "4-Mar", "5-Mar", "6-Mar" , "7-Mar", "8-Mar", "9-Mar", "10-Mar"))
normalized_expression <- normalized_expression[-index,]
rownames(normalized_expression) <- normalized_expression[,1] #43881*312
normalized_expression <- normalized_expression[,-1]
#0.3188041 68.3377700 --- No log2ization
saveRDS(normalized_expression, file = "/data/active_data/lzl/RenalTumor-20200713/Data/ICB.therapy/normalized_expression.rds")

patient.info <- read.table("/data/active_data/lzl/RenalTumor-20200713/Data/ICB.therapy/Checkmate.009.010.025.patient.info.txt", header = T, stringsAsFactors = F, sep = "\t")
# 1006*121
##Corresponding sample clinical information
index <- match(colnames(normalized_expression), patient.info$RNA_ID)
patient.info.RNA <- patient.info[index,] # 311*121
write.table(patient.info.RNA, file = "/data/active_data/lzl/RenalTumor-20200713/Data/ICB.therapy/Checkmate.009.010.025.patient.RNA.info.txt", row.names = F , col.names = T, sep = "\t", quote = F)
saveRDS(patient.info.RNA, file = "/data/active_data/lzl/RenalTumor-20200713/Data/ICB.therapy/patient.info.RNA.rds")

#PCAAnalysis
library(factoextra)
res.pca <- prcomp(t(normalized_expression), scale = TRUE)
pdf("/data/active_data/lzl/RenalTumor-20200713/Data/ICB.therapy/PCA.pdf")
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
             )
dev.off()