#' @description: Endothelium cells

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
library(dplyr)
library(future)
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 50000 * 1024^2) #50G

source(file = "/home/longzhilin/Analysis_Code/Visualization/colorPalettes.R")
setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC")
scATAC.data <- readRDS("scATAC.data.pro.rds")

###########################1.differential peak#################################
DefaultAssay(scATAC.data) <- "Peaks"
Idents(scATAC.data) <- scATAC.data$AnnotatedcellType
idents <- as.character(levels(scATAC.data))
Endothelium.DARs <- FindMarkers(scATAC.data, 
                             test.use = 'LR',
                             ident.1 = "Endothelium (VCAM1+)",
                             ident.2 = "Endothelium (VCAM1-)",
                             logfc.threshold = 0,
                             min.pct = 0.05, # often necessary to lower the min.pct threshold
                             latent.vars = "peak_region_fragments")
cf <- ClosestFeature(scATAC.data, regions = rownames(Endothelium.DARs)) # Find the closest feature to a given set of genomic regions
Endothelium.DARs <- cbind(Endothelium.DARs, gene=cf$gene_name, gene_biotype = cf$gene_biotype, type = cf$type, distance=cf$distance)
Endothelium.DARs <- arrange(Endothelium.DARs, desc(avg_log2FC))
write.xlsx(Endothelium.DARs, file = "4.Peak/Endothelium.all.DARs.xlsx", rowNames = F)
saveRDS(Endothelium.DARs, file = "4.Peak/Endothelium.all.DARs.rds")

Endothelium.sig.DARs <- Endothelium.DARs %>% filter(abs(avg_log2FC) >= 0.25 & p_val_adj < 0.05)
Endothelium.sig.DARs$genomicRegion <- rownames(Endothelium.sig.DARs)
#plot--- differentially accessible chromatin regions heatmap
sig.region <- Endothelium.sig.DARs %>% select(genomicRegion) %>% distinct() 
scATAC.Endothelium.data <- subset(scATAC.data, subset = AnnotatedcellType %in% c("Endothelium (VCAM1+)", "Endothelium (VCAM1-)"))
#average fragment of each peak in each cell type
sig.region.mean <- AverageExpression(scATAC.Endothelium.data, features = sig.region$genomicRegion, assays = "Peaks")
sig.region.mean.scale <- scale(t(sig.region.mean$Peaks))
pdf("4.Peak/Endothelium.sig.pos.DAR.pdf", height = 6, width = 6)
Heatmap(sig.region.mean$Peaks, name = "z-score", show_column_dend = F, show_row_dend = F, show_row_names = F, row_names_gp = gpar(fontsize = 10), width = unit(10, "cm"), height = unit(8, "cm"))
Heatmap(sig.region.mean.scale, name = "z-score", show_column_dend = F, show_row_dend = F, show_column_names = F, row_names_gp = gpar(fontsize = 10), width = unit(10, "cm"), height = unit(8, "cm"))
dev.off()

### differentially accessible chromatin and genes
## load DEGs
Endothelium.DEGs <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/2.Cluster/AnnotateCellType/Endothelium.all.DEGs.rds")
Endothelium.DEGs <- Endothelium.DEGs %>% filter(abs(avg_log2FC) > 0.25 & p_val_adj < 0.05)

# overlap gene
common.gene <- intersect(Endothelium.DEGs$gene, Endothelium.sig.DARs$gene)
Endothelium.sig.DARs.common <- Endothelium.sig.DARs[which(Endothelium.sig.DARs$gene %in% common.gene),]
Endothelium.DEGs.common <- Endothelium.DEGs[which(Endothelium.DEGs$gene %in% common.gene),]

#### volcano map
source("/home/longzhilin/Analysis_Code/Visualization/Plot.EnhancedVolcano.R")
pdf("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/2.Cluster/AnnotateCellType/Endothelium.EnhancedVolcano.pdf")
Plot.EnhancedVolcano(Endothelium.DEGs, x = "avg_log2FC", y = "p_val_adj", selected.showType = "Both", select.num = 0, interestGenes = Endothelium.DEGs.common$gene[which(abs(Endothelium.DEGs.common$avg_log2FC)>=1)], title = "Endothelium (VCAM1+) VS Endothelium (VCAM1-)")
dev.off()

###########################2.differential motif#################################
DefaultAssay(scATAC.data) <- 'chromvar'
motif.info <- scATAC.data@assays$Peaks@motifs
Endothelium.chromVAR <- FindMarkers(
    object = scATAC.data,
    ident.1 = 'Endothelium (VCAM1+)',
    ident.2 = 'Endothelium (VCAM1-)',
    test.use = 'LR',
    logfc.threshold = 0, 
    latent.vars = "nCount_Peaks"
)
Endothelium.chromVAR <- arrange(Endothelium.chromVAR, desc(avg_log2FC))
motifs <- gsub("-", "_", rownames(Endothelium.chromVAR))
motifNames <- sapply(motifs, function(x) motif.info@motif.names[[x]])
Endothelium.chromVAR$motifName <- motifNames
write.xlsx(Endothelium.chromVAR, file = "4.Peak/Endothelium.chromVAR.xlsx", rowNames = F)
saveRDS(Endothelium.chromVAR, file = "4.Peak/Endothelium.chromVAR.rds")