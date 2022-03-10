#' @description: identify the target genes of TF

library(openxlsx)
library(circlize)
library(dplyr)
library(tibble)
library(data.table)
library(Signac)
library(Seurat)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVARmotifs)
data("human_pwms_v2")
library(ggpubr)

setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC")
scATAC.data <- readRDS("scATAC.data.pro.rds")
DefaultAssay(scATAC.data) <- "Peaks"
scRNA.data <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/data.merge.pro.rds")
DefaultAssay(scRNA.data) <- "RNA"

celltype <- "Tumor"
cCREs <- readRDS(paste0("6.Co-Accessible/", celltype,"_peak_gene_correlation.rds")) # 4368 cCREs
#gene-cCRE link --- enhancer
cCREs <- cCREs %>% subset(distance_bp/1000>=10 & peak2_type != 'Promoter')
cCREs$target_gene <- cCREs$peak1_nearestGene
cCREs$peak_gene <- paste0(cCREs$Peak2, '_', cCREs$peak1_nearestGene)
cCREs$Peak1 <- gsub("_", "-", cCREs$Peak1)
cCREs$Peak2 <- gsub("_", "-", cCREs$Peak2)

# motif information
motif.info <- data.frame(originName = names(scATAC.data@assays$Peaks@motifs@motif.names), TF = unlist(scATAC.data@assays$Peaks@motifs@motif.names))
rownames(motif.info) <- NULL
#                           originName     TF
#   ENSG00000008196_LINE2_TFAP2B_D_N1 TFAP2B
#      ENSG00000008197_LINE6_TFAP2D_D TFAP2D
#   ENSG00000087510_LINE7_TFAP2C_D_N3 TFAP2C
#     ENSG00000116819_LINE19_TFAP2E_I TFAP2E
#  ENSG00000137203_LINE29_TFAP2A_D_N3 TFAP2A
#  ENSG00000116017_LINE44_ARID3A_D_N1 ARID3A

# subset by celltype
scATAC.sub <- subset(scATAC.data, cellType_low == celltype)
scRNA.sub <- subset(scRNA.data, cellType_low == celltype)

scATAC.sub <- FindTopFeatures(scATAC.sub, min.cutoff=ncol(scATAC.sub)*0.05) # present in what % of cells, 要求95%的细胞都包含peak
scRNA.sub <- FindVariableFeatures(scRNA.sub, nfeatures = 3000)

# use chipseeker annotate each peak in scATAC-seq
peak.info <- readRDS("4.Peak/peak.annotation.simple.ChIPseeker.rds")

# TF target genes
# for each motif, find genes:
use_variable_genes <- F
motif_list <- list()
edge_df.list <- data.frame()
vertex_df.list <- data.frame()
for(i in 1:nrow(motif.info)){
  motif_name <- motif.info[i,2]
  motif_ID <- motif.info[i,1]
  cat("Analysis ", motif_ID, "\n")

  # get list of promoter and enhancer targets of these TFs
  DefaultAssay(scATAC.data) <- "Peaks"
  motif_accessible <- Motifs(scATAC.data)@data[,motif_ID]
  motif_accessible <- names(motif_accessible)[motif_accessible > 0]
  motif_accessible <- motif_accessible[motif_accessible %in% VariableFeatures(scATAC.sub)]

  # promoter locate in 1kb
  motif_accessible_promoters <- motif_accessible[motif_accessible %in% peak.info$peaks[which(peak.info$originType == "Promoter (<=1kb)")]]
  motif_target_genes <- peak.info$SYMBOL[match(motif_accessible_promoters, peak.info$peaks)]
  motif_target_genes <- as.character(motif_target_genes)

  # variable genes only?
  if(use_variable_genes){
    motif_target_genes <- motif_target_genes[motif_target_genes %in% VariableFeatures(scRNA.sub)] %>% as.character %>% unique
  } else{
    motif_target_genes<- motif_target_genes %>% as.character %>% unique
  }

  # enhancer target genes
  motif_accessible_enhancers <- motif_accessible[motif_accessible %in% cCREs$Peak2]
  motif_cCREs <- subset(cCREs, Peak2 %in% motif_accessible_enhancers)
  motif_enhancers_target_genes <- motif_cCREs$target_gene

  # variable genes only?
  if(use_variable_genes){
    motif_enhancers_target_genes <- motif_enhancers_target_genes[motif_enhancers_target_genes %in% VariableFeatures(scRNA.sub)] %>% as.character %>% unique
  }else{
    motif_enhancers_target_genes<- motif_enhancers_target_genes %>% as.character %>% unique
  }

  vertex_df <- data.frame(name = c(motif_name, as.character(unique(c(motif_enhancers_target_genes, motif_target_genes)))))

  # check if there are promoter targets:
  if(length(motif_target_genes) > 0){
    promoter_edge_df <- data.frame(
      from = motif_name,
      to = as.character(motif_target_genes),
      type = 'promoter'
    )
  } else{promoter_edge_df <- data.frame()}

  # check if there are enhancer targets:
  if(length(motif_enhancers_target_genes) > 0){
    enhancer_edge_df <- data.frame(
      from = motif_name,
      to = as.character(motif_enhancers_target_genes),
      type = 'enhancer'
    )
  } else{enhancer_edge_df <- data.frame()}

  #cur_edge_df <- rbind(cur_promoter_edge_df, cur_enhancer_edge_df, cur_repressors_edge_df)
  edge_df <- rbind(promoter_edge_df, enhancer_edge_df)

  edge_df.list <- rbind(edge_df.list, edge_df)
  vertex_df.list <- rbind(vertex_df.list, vertex_df)
}

vertex_df.list <- data.frame(name=na.omit(as.character(unique(vertex_df.list$name))))
vertex_df.list$name <- as.character(vertex_df.list$name)
edge_df.list <- na.omit(edge_df.list)
saveRDS(edge_df.list, file = paste0("7.TF.analysis/targets/", celltype,".All.TF.targets.rds"))