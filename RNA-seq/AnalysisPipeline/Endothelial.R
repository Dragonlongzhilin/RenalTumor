#' @description: Endothelium cells

library(Seurat)
library(harmony)
library(clustree)
library(ggpubr)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(vegan)
set.seed(101)
library(openxlsx)
library(future)
plan("multiprocess", workers = 10) 
options(future.globals.maxSize = 50000 * 1024^2) # set 50G RAM
setwd(dir = "/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA")
source(file = "/home/longzhilin/Analysis_Code/Visualization/colorPalettes.R")
source(file = "/home/longzhilin/Analysis_Code/code/ratio.plot.R")

scRNA.data <- readRDS("data.merge.pro.rds")

###########################1.differential gene#################################
Idents(scRNA.data) <- scRNA.data$cellType_low
idents <- as.character(levels(scRNA.data))
Endothelium.markers <- FindMarkers(scRNA.data, 
                                    ident.1 = "Endothelium (VCAM1+)",
                                    ident.2 = "Endothelium (VCAM1-)",
                                    logfc.threshold = 0, 
                                    min.pct = 0.1, 
                                    test.use = "MAST", 
                                    latent.vars = "orig.ident")
Endothelium.markers <- arrange(Endothelium.markers, desc(avg_log2FC))
Endothelium.markers$gene <- rownames(Endothelium.markers)
Endothelium.DEGs <- Endothelium.markers %>% filter(abs(avg_log2FC) >= 0.25 & p_val_adj < 0.05)
write.xlsx(Endothelium.markers, file = "2.Cluster/AnnotateCellType/Endothelium.all.DEGs.xlsx", rowNames = F, overwrite = T)
saveRDS(Endothelium.markers, file = "2.Cluster/AnnotateCellType/Endothelium.all.DEGs.rds")

#### volcano map
source("/home/longzhilin/Analysis_Code/Visualization/Plot.EnhancedVolcano.R")
pdf("2.Cluster/AnnotateCellType/Endothelium.EnhancedVolcano.pdf")
Plot.EnhancedVolcano(Endothelium.DEGs, x = "avg_log2FC", y = "p_val_adj", selected.showType = "Both", select.num = 15, title = "Endothelium (VCAM1+) VS Endothelium (VCAM1-)")
dev.off()

#### pathway enrichment
source("/home/longzhilin/Analysis_Code/PathwayEnrichment/clusterProfiler.enricher.R")
up.DEGs <- filter(Endothelium.DEGs, avg_log2FC > 1 & p_val_adj < 0.05)
down.DEGs <- filter(Endothelium.DEGs, avg_log2FC< -1 & p_val_adj < 0.05)
DEGs <- list(up = up.DEGs$gene, down = down.DEGs$gene)
pdf("2.Cluster/AnnotateCellType/Endothelium.program.pathway.enrichment.pdf")
res <- lapply(names(DEGs), function(x){
    y <- DEGs[[x]]
    res <- cluterProfiler.enricher(gene = y, geneType = "SYMBOL", db.type = "MsigDB", saveDir = paste0(getwd(),"/2.Cluster/AnnotateCellType"),
                            title = x, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
    # gene ratio
    res <- res$em.res.genesymbol@result %>% filter(p.adjust<0.05) #fdr adjust
    pathway.geneNum <- unlist(strsplit(res$BgRatio, "/"))[seq(1, nrow(res),by=2)]
    gene.ratio <- as.numeric(res$Count)/as.numeric(pathway.geneNum)
    res$gene.ratio <- gene.ratio
    return(res)
})
dev.off()

####cluster profiler enrichment result display
library(ggthemes)
#top10
enrich1 <- read.csv("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/2.Cluster/AnnotateCellType/up_enricherResult.csv", header = T)
enrich2 <- read.csv("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/2.Cluster/AnnotateCellType/down_enricherResult.csv", header = T)
enrich.list <- list(up = enrich1, down = enrich2)
enrich <- list(up = enrich1[c(1,3,8,9,16,17,18,21,22,23,25,26,27,28,42,49,50,72), c("ID")],
               down = enrich2[c(2:5,8:10,11,12,19,21,22), c("ID")])
enrich.path <- Reduce(function(x,y) c(x,y), enrich)
enrich.path <- unique(enrich.path)
enrich.path <- enrich.path[-c(3,6,26,28,29,4,13,11,5,20,23,27,30)]
enrich.list.select <- lapply(enrich.list, function(x){
  idx <- which(x$ID %in% enrich.path)
  return(x[idx, c("ID", "Count", "pvalue", "p.adjust")])
})
enrich.res <- Reduce(function(x,y) rbind(x,y), enrich.list.select)
a <- sapply(enrich.list.select, function(x){nrow(x)})
enrich.res$Type <- c(rep("Endothelium (VCAM1+) ", as.numeric(a[1])), rep("Endothelium (VCAM1-) ", as.numeric(a[2])))
enrich.res$p.adjust <- -log10(enrich.res$p.adjust)

# process gene name
enrich.res$ID <- gsub("_", " ", enrich.res$ID)
# factor ID
idx1 <- na.omit(match(enrich$up, enrich.path))
idx2 <- na.omit(match(enrich$down, enrich.path))
enrich.path <- enrich.path[c(idx1, idx2)]
enrich.res$ID <- factor(enrich.res$ID, levels = gsub("_", " ", enrich.path))
#source(file = "/home/longzhilin/Analysis_Code/Visualization/ggplot.bubble.plot.R")
pdf("2.Cluster/AnnotateCellType/Endothelium.enrichment.bubble.pdf", width = 5, height = 5)
p1 <- ggplot(enrich.res, 
            aes(x = Type, 
                y = ID, 
                size = Count, 
                fill = p.adjust))
p2 <- p1 + guides(color=FALSE) + geom_point(shape = 21, alpha = 0.7) + theme_few() + scale_size_continuous(breaks = c(3,5,7,9,11), labels = c(3,5,7,9,11)) + scale_fill_gradientn(colours = c("#fed71b", "#fcc41c", "#f27620", "#e92925"), breaks = c(1.69897, 4, 6, 8), labels = c(0.02, 0.0001, 0.000001, 0.00000001))
p2 <- p2 + theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1, size = 8), axis.text.y = element_text(size = 8)) + xlab("") + ylab("")
print(p2)
dev.off()
