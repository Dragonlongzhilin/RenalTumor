#' @description: Motif enrichment with chromVAR

library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
set.seed(101)
library(ggpubr)
library(openxlsx)
library(chromVARmotifs)
data("human_pwms_v2")
library(ComplexHeatmap)
library(circlize)
library(future)
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 50000 * 1024^2) #50G

## We will explore two complementary options for performing motif analysis: 
#one by finding overrepresented motifs in a set of differentially accessible peaks, 
#one method performing differential motif activity analysis between groups of cells.

source(file = "/home/longzhilin/Analysis_Code/Plot_colorPaletters.R")
setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC")

scATAC.data <- readRDS("scATAC.data.pro.rds")
Idents(scATAC.data) <- scATAC.data$AnnotatedcellType
DefaultAssay(scATAC.data) <- "Peaks"

# Get a list of motif position frequency matrices from the JASPAR database
# pfm <- getMatrixSet(
#   x = JASPAR2020,
#   opts = list(species = 9606, all_versions = FALSE)
# )

# add motif information
scATAC.data <- AddMotifs(
  object = scATAC.data,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = human_pwms_v2
)

# Finding overrepresented motifs
GetMotifs <- function(cellType, scATAC.object) {
  print(paste0("Finding motifs for: ",cellType))
  overrepresented.motifs <- FindMarkers(scATAC.object,
                                        ident.1 = cellType,
                                        test.use = 'LR',
                                        min.pct = 0.05,
                                        latent.vars = "nCount_Peaks")
  enriched.motifs <- FindMotifs(object = scATAC.object, features = rownames(overrepresented.motifs[overrepresented.motifs$p_val_adj < 0.05, ]))
  return(enriched.motifs)
}
# FindMarkers and write to an xlsx file with default parameters
idents <- as.character(levels(Idents(scATAC.data)))
cellType.motifs <- lapply(idents, function(x) GetMotifs(x, scATAC.object = scATAC.data))
names(cellType.motifs) <- idents
write.xlsx(cellType.motifs, file = "5.Motif/motifs.celltype.human_pwms_v2.xlsx", sheetName = idents, rowNames = T)
saveRDS(cellType.motifs, file = "5.Motif/cellType.motifs.human_pwms_v2.rds")

library(ggplot2)
library(ggseqlogo)
PWMatrixToProbMatrix <- function(x){
        if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
        m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
        m <- t(t(m)/colSums(m))
        m
}

pdf("5.Motif/cellType.specific.TFs.logo.pdf")
cellType.TFs.PPM <- sapply(idents, function(x){
    cellType.top.TFs <- head(rownames(cellType.motifs[[x]]))
    cellType.top.TFs <- gsub("-", "_", cellType.top.TFs)
    index <- match(cellType.top.TFs, names(human_pwms_v2))
    PPM.list <- lapply(index, function(y){
        PPM <- PWMatrixToProbMatrix(human_pwms_v2[[y]])
    })
    names(PPM.list) <- as.character(unlist(scATAC.data@assays$Peaks@motifs@motif.names[index]))
    p <- ggseqlogo(PPM.list) + ggtitle(x) + theme(plot.title = element_text(hjust = 0.5))
    print(p)
    return(PPM.list)
})
dev.off()

##########################################ChromVAR
#We can also compute a per-cell motif activity score by running chromVAR. 
#This allows us to visualize motif activities per cell, and also provides an alternative method of identifying differentially-active motifs between cell types.

# Computing motif activities
scATAC.data <- RunChromVAR(
  object = scATAC.data,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
DefaultAssay(scATAC.data) <- 'chromvar'

#cellType.chromvar.activity <- FindAllMarkers(scATAC.data, group.by = "AnnotatedcellType", test.use = 'LR', latent.vars = "nCount_Peaks")

library(tidyverse)
GetChromvarActivities <- function(cellType, scATAC.object, motif.info) {
  print(paste0("Finding chromVAR activities for: ",cellType))
  differential.activity <- FindMarkers(scATAC.object, 
                     ident.1 = cellType,
                     test.use = 'LR',
                     logfc.threshold = 0, 
                     latent.vars = "nCount_Peaks") 
  motifs <- gsub("-", "_", rownames(differential.activity))
  motifNames <- sapply(motifs, function(x) motif.info@motif.names[[x]])
  return(cbind(differential.activity, gene = motifNames))
}
idents <- as.character(levels(Idents(scATAC.data)))
cellType.motifs.chromVAR <- lapply(idents, function(x) GetChromvarActivities(x, scATAC.object = scATAC.data, motif.info = scATAC.data@assays$Peaks@motifs))
names(cellType.motifs.chromVAR) <- idents
cellType.motifs.chromVAR <- lapply(cellType.motifs.chromVAR, function(x){
  up.x <- x %>% filter(avg_log2FC>0) %>% arrange(desc(avg_log2FC))
  down.x <- x %>% filter(avg_log2FC<=0) %>% arrange(avg_log2FC)
  x <- rbind(up.x, down.x)
  return(x)
})
names(cellType.motifs.chromVAR) <- idents
write.xlsx(cellType.motifs.chromVAR, file = "5.Motif/cellType.motifs.chromVAR.human_pwms_v2.xlsx", sheetName = idents, rowNames = T)
saveRDS(cellType.motifs.chromVAR, file = "5.Motif/cellType.motifs.chromVAR.human_pwms_v2.rds")

pdf("5.Motif/cellType.specific.TFs.chromVAR.logo.pdf")
cellType.TFs.PPM <- sapply(idents, function(x){
    cellType.top.TFs <- head(rownames(cellType.motifs.chromVAR[[x]]))
    cellType.top.TFs <- gsub("-", "_", cellType.top.TFs)
    index <- match(cellType.top.TFs, names(human_pwms_v2))
    PPM.list <- lapply(index, function(y){
        PPM <- PWMatrixToProbMatrix(human_pwms_v2[[y]])
    })
    names(PPM.list) <- as.character(unlist(scATAC.data@assays$Peaks@motifs@motif.names[index]))
    p <- ggseqlogo(PPM.list) + ggtitle(x) + theme(plot.title = element_text(hjust = 0.5))
    print(p)
    return(PPM.list)
})
dev.off()
saveRDS(scATAC.data, "scATAC.data.pro.rds")
