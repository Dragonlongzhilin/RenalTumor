#' @description: identify the cis-coassessibility peaks

library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(cicero)
library(monocle3)
set.seed(101)
library(ggpubr)
library(openxlsx)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(data.table)
library(future)
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 50000 * 1024^2) #50G

setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC")

scATAC.data <- readRDS("scATAC.data.pro.rds")
Idents(scATAC.data) <- scATAC.data$AnnotatedcellType
DefaultAssay(scATAC.data) <- "Peaks"

# input Seurat object
CreateCiceroCDS <- function(seurat_atac_obj, assay = "Peaks", reduction_method = c("UMAP", "tSNE")) {
    DefaultAssay(seurat_atac_obj) <- assay
    count_data <- GetAssayData(seurat_atac_obj, slots = "counts")
    summ <- summary(count_data)
    rownames(count_data) <- gsub("-", "_", rownames(count_data))
    summ_frame <- data.frame(peak = rownames(count_data)[summ$i],
                            cell.id = colnames(count_data)[summ$j],
                            count = summ$x)

    # create cell data set object with cicero constructor
    input_cds <- make_atac_cds(summ_frame, binarize = TRUE)
    input_cds <- detect_genes(input_cds)
    input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 
    input_cds <- estimate_size_factors(input_cds)
    input_cds <- preprocess_cds(input_cds, method="LSI")
    input_cds <- reduce_dimension(input_cds, reduction_method=reduction_method, preprocess_method="LSI")

    if(reduction_method == "UMAP"){
        umap_coords <- reducedDims(input_cds)$UMAP
        cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates=umap_coords)
    }

    if(reduction_method == "tSNE"){
        umap_coords <- reducedDims(input_cds)$tSNE
        cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates=umap_coords)
    }
    return(cicero_cds)
}

FindCiceroConnection <- function(seurat_atac_obj, cds, chrom = NULL) {
    # require default assay is peak assay
    
    # get the chromosome sizes from the Seurat object
    genome <- seqlengths(seurat_atac_obj)

    # run on the whole genome
    levels <- paste0("chr",c(seq(1,22),"X","Y"))
    genome <- genome[levels]

    # convert chromosome sizes to a dataframe
    genome.df <- data.frame("chr" = names(genome), "length" = genome)
    conns <- run_cicero(cds, genomic_coords = genome.df) 
    return(conns)
}

Get_Ccans <- function(seurat_atac_obj, clusterID = NULL, reduction_method = "UMAP") {
  if(!is.null(clusterID)) {
    print(paste0("Subsetting seurat object for: ",clusterID))
    seurat_atac_obj <- subset(seurat_atac_obj, ident = clusterID) # create a subset
  } 
  # convert seurat objects into cicero cell datasets in preparation for detecting cicero connections
  print("Preparing Cicero CDS")
  ciceroCds <- CreateCiceroCDS(seurat_atac_obj, reduction_method = reduction_method)
  
  # generate disease-specific CCANS for all chromsomes of a particular celltype
  print("Finding Cicero connections")
  conns <- FindCiceroConnection(seurat_atac_obj = seurat_atac_obj, cds = ciceroCds)
  
  # Finding cis-Co-accessibility Networks (CCANS)
  CCAN_assigns <- generate_ccans(conns)

  # create a column that identifies which connections belong to a CCAN
  ccan1 <- left_join(conns, CCAN_assigns, by=c("Peak1" = "Peak"), all.x=TRUE)
  colnames(ccan1)[4] <- "CCAN1"
  ccan2 <- left_join(conns, CCAN_assigns, by=c("Peak2" = "Peak"), all.x=TRUE)
  colnames(ccan2)[4] <- "CCAN2"
  df <- cbind(ccan1, CCAN2=ccan2$CCAN2) %>%
    dplyr::mutate(CCAN = ifelse(CCAN1 == CCAN2, CCAN1, 0)) %>%
    dplyr::select(-CCAN1, -CCAN2)
  res <- list(ciceroCds = ciceroCds, conns = conns, CCAN_assigns = CCAN_assigns, df = df)
  return(res)
}

####################################################################
# subset the snATACseq object by disease state within celltype of interest and find ccans
idents <- levels(scATAC.data)
# remove mast cell due to the cell number less than 100
idents <- idents[-12]
list.ccan <- lapply(idents, function(ident) {Get_Ccans(ident, seurat_atac_obj = scATAC.data)})

# write to file
names(list.ccan) <- idents
lapply(names(list.ccan), function(x) {
  fwrite(list.ccan[[x]]$df, file = paste0("6.Co-Accessible/ccans/ciceroConns.",x,".csv"), row.names = F)
})
saveRDS(list.ccan, file = "6.Co-Accessible/ccans/cellType.ccans.rds")

# calculate a global CCAN for all celltypes
ccan <- Get_Ccans(seurat_atac_obj = scATAC.data)
fwrite(ccan$df, file = "6.Co-Accessible/ccans/ciceroConns.allcells.csv", row.names = F)
saveRDS(ccan, file = "6.Co-Accessible/ccans/All.ccans.rds")

#### extract the tumor-ccan
tumor.ciceroCds <- list.ccan$Tumor$ciceroCds
tumor.conns <- list.ccan$Tumor$conns
tumor.ccans <- list.ccan$Tumor$CCAN_assigns
saveRDS(tumor.ciceroCds, file = "6.Co-Accessible/ccans/tumor.ciceroCds.rds")
saveRDS(tumor.conns, file = "6.Co-Accessible/ccans/tumor.conns.rds")