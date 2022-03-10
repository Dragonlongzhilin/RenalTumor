#' @description: Speculate the copy number change of cancer cells

####Method 1: inferCNV
#2021.4.14
###Create GRCh38 gene location file --- Stored on 231 server
#Reference: https://github.com/broadinstitute/inferCNV/wiki/instructions-create-genome-position-file
cd /data/activate_data/longzhilin/RenalTumor-20200713/inferCNV_result
#Homo_sapiens.GRCh38.100.gtf
python /home/longzhilin/softwares/infercnv-master/scripts/gtf_to_position_file.py --attribute_name gene_name /home/longzhilin/softwares/ReferenceFasta/Annotation/Homo_sapiens.GRCh38.100.gtf GRCh38_gen_pos.txt
#'@note Note that you need to modify the GRCh38_gen_pos.txt file
#'1. Add "chr" to the second column to mark the chromosome
#'2. Refer to the gencode_v19_gene_pos.txt provided by it to remove genes other than chr1-22, X, Y, MT
#'3. It should be noted that the order of chromosomes may be reversed. In gencode_v19_gene_pos.txt it is chr1-22,MT,X,Y, but it seems that the order of chromosomes is not mandatory, so no additional processing is currently performed.
#R code
genecode.modify <- read.table("/data/activate_data/longzhilin/RenalTumor-20200713/inferCNV_result/GRCh38_gen_pos_modify.txt", header = F, stringsAsFactors = F, sep = "\t") #60591
chrom.order <- c(paste0("chr",1:22), "chrM", "chrX", "chrY")
index <- sapply(chrom.order, function(x){
    index <- which(genecode.modify[,2] == x)
    return(index)
})
index <- unlist(index)
genecode.modify.order <- genecode.modify[index,]
write.table(genecode.modify.order, "/data/activate_data/longzhilin/RenalTumor-20200713/inferCNV_result/GRCh38_gen_pos_modify_order.txt", row.name = F, col.names = F, quote = F, sep = "\t ")

####Build the file to be imported on the 91 server
##Predict the CNV of potential cancer cells based on the annotated normal cells
library(Seurat)
data.merge <- readRDS(file = "/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210721/scRNA/data.merge.harmony.stardard.PC50.feature3000.rds")

##Extract NK/NKT cell and CD8+ T cell as reference cells
index <- which(data.merge@meta.data$cellType %in% c("CD8+ T cell", "NK/NKT cell"))
Normal.cell <- rownames(data.merge@meta.data)[index]
Normal.cellType <- as.character(data.merge@meta.data$cellType[index])

###The first model: cancer cells
##Predict the CNV of potential cancer cells based on the annotated normal cells
library(Seurat)
data.merge <- readRDS(file = "/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/data.merge.pro.rds")

##Extract NK/NKT cell and CD8+ T cell as reference cells
index <- which(data.merge@meta.data$cellType %in% c("CD8+ T cell", "NK/NKT cell"))
Normal.cell <- rownames(data.merge@meta.data)[index] # 9949 cells
Normal.cellType <- as.character(data.merge@meta.data$cellType[index])

##Need to guess the cell population of the copy number
setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/3.CNV")
possible.cancer.cell <- rownames(data.merge@meta.data)[which(data.merge@meta.data$seurat_clusters %in% c("4"))] # 3564 cells
DefaultAssay(data.merge) <- "RNA"
cancer.count <- GetAssayData(data.merge, slot = "counts")
counts_matrix <- cancer.count[, c(Normal.cell, possible.cancer.cell)] #26118 13513
cell.id <- gsub("-", "\\.", colnames(counts_matrix))
colnames(counts_matrix) <- cell.id
write.table(counts_matrix, file = "countmatrix.txt", sep = "\t", row.names = T, col.names = T, quote = F)

##Building cell type annotation file
index <- match(possible.cancer.cell, rownames(data.merge@meta.data))
cellLabel1 <- as.character(data.merge$orig.ident[index])
cellLabel1 <- paste0("Sample_", cellLabel1)
annoFile <- data.frame(cell.id = cell.id, celltype = c(Normal.cellType, cellLabel1))
write.table(annoFile, file = "annoFile.txt", col.names = F, row.names = F, quote = F, sep = "\t")

###231 Server
setwd("/data/activate_data/longzhilin/RenalTumor-20200713/inferCNV_result/RNA_20210803")
library(infercnv) #Version:1.6
infercnv_obj = CreateInfercnvObject(raw_counts_matrix= "/data/activate_data/longzhilin/RenalTumor-20200713/inferCNV_result/RNA_20210803/countmatrix.txt",
                                    annotations_file= "/data/activate_data/longzhilin/RenalTumor-20200713/inferCNV_result/RNA_20210803/annoFile.txt",
                                    delim="\t",
                                    gene_order_file="/data/activate_data/longzhilin/RenalTumor-20200713/inferCNV_result/GRCh38_gen_pos_modify_order.txt",
                                    ref_group_names=c("CD8+ T cell", "NK/NKT cell"))
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,# cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir= "/data/activate_data/longzhilin/RenalTumor-20200713/inferCNV_result/RNA_20210803",
                             cluster_by_groups=T,
                             denoise=TRUE,
                             HMM=TRUE,
                             num_threads = 48)