#' @description: peak annotation and enrichment analysis

##reference：
#https://hbctraining.github.io/Intro-to-ChIPseq/lessons/12_functional_analysis.html
#http://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html
library(Signac)
library(Seurat)
library(ggpubr)
library(GenomicRanges)
library(tibble)
library(dplyr)
set.seed(101) 

setwd(dir = "/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC")
scATAC.data <- readRDS("scATAC.data.pro.rds")
DefaultAssay(scATAC.data) <- "Peaks"

## loading packages
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)
library(org.Hs.eg.db)

#####################################All Peak Annotation
peakSet <- GetAssayData(scATAC.data, slot = "data") # 212326 x 21271
All.gr <- StringToGRanges(rownames(peakSet))
#TSS (transcription start site) region, by default TSS is defined from -3kb to +3kb
All.peakAnno <- annotatePeak(All.gr, tssRegion=c(-3000, 3000), TxDb = txdb, annoDb="org.Hs.eg.db")
saveRDS(All.peakAnno, file = "4.Peak/peak.annotation.ChIPseeker.rds")
peakAnno.info <- All.peakAnno@anno
peakAnno.info$annotation <- gsub(" .*", "", peakAnno.info$annotation)
peakAnno.info$annotation <- gsub("'", "_UTR", peakAnno.info$annotation)
peak.info <- data.frame(peaks = rownames(Motifs(scATAC.data)@data),
                        originType = All.peakAnno@anno$annotation,
                        peakType = peakAnno.info$annotation, 
                        distanceToTSS = peakAnno.info$distanceToTSS,
                        SYMBOL = peakAnno.info$SYMBOL)
saveRDS(peak.info, file = "4.Peak/peak.annotation.simple.ChIPseeker.rds")
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(All.gr, windows=promoter)

pdf("4.Peak/Annotation/All.peakAnnotation.pdf")
plotAnnoPie(All.peakAnno)
plotAnnoBar(All.peakAnno)
vennpie(All.peakAnno) #some annotation overlap
upsetplot(All.peakAnno, vennpie=TRUE)
#Visualize distribution of TF-binding loci relative to TSS
plotDistToTSS(All.peakAnno, title="Distribution of transcription factor-binding loci\nrelative to TSS")
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="#4DBBD5")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)
dev.off()

#######################################annotated the cell-type specific peaks
cellType.sig.pos.DARs <- readRDS("4.Peak/cellType.sig.pos.DARs.rds")
cellType.sig.pos.gr <- StringToGRanges(unique(cellType.sig.pos.DARs$genomicRegion)) # 8797 DAR peaks
DAR.Anno <- annotatePeak(cellType.sig.pos.gr, tssRegion=c(-3000, 3000), TxDb = txdb, annoDb="org.Hs.eg.db")

idents <- unique(cellType.sig.pos.DARs$cellType)
cellType.gr.list <- lapply(idents, function(x){
  idx <- which(cellType.sig.pos.DARs$cellType == x)
  gr <- StringToGRanges(cellType.sig.pos.DARs$genomicRegion[idx])
})
names(cellType.gr.list) <- idents

#annotated DAR for each cell type
DAR.anno.list <- lapply(cellType.gr.list, function(x){
    DAR.cell <- annotatePeak(x, tssRegion=c(-3000, 3000), TxDb = txdb, annoDb="org.Hs.eg.db")
    return(DAR.cell)
})
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(cellType.sig.pos.gr, windows=promoter)

pdf("4.Peak/Annotation/cellType.DAR.annotation.pdf")
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="#4DBBD5")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)
plotAnnoPie(DAR.Anno)
plotAnnoBar(DAR.Anno)
vennpie(DAR.Anno) #some annotation overlap
upsetplot(DAR.Anno, vennpie=TRUE)

#Visualize distribution of TF-binding loci relative to TSS
plotDistToTSS(DAR.Anno, title="Distribution of transcription factor-binding loci\nrelative to TSS")

plotAnnoBar(DAR.anno.list)
plotDistToTSS(DAR.anno.list)
dev.off()

############################# write DARs to bed files
# write all peaks to file to serve as background set:
DefaultAssay(scATAC.data) <- "Peaks"
peak_ranges <- StringToGRanges(rownames(scATAC.data)) %>% sort
write.table(as.data.frame(peak_ranges)[,1:3], file='4.Peak/bed/allPeaks.bed', row.names=F, col.names=F, sep='\t', quote=F)

All.peakAnno <- readRDS("4.Peak/peak.annotation.ChIPseeker.rds")
peakAnno.info <- All.peakAnno@anno
peakAnno.info$annotation <- gsub(" .*", "", peakAnno.info$annotation)
peakAnno.info$annotation <- gsub("'", "_UTR", peakAnno.info$annotation)
peakAnno.info.frame <- as.data.frame(peakAnno.info)
peakAnno.info.frame$id <- paste0(peakAnno.info.frame[,1], "-", peakAnno.info.frame[,2], "-", peakAnno.info.frame[,3])

# Distal type
distal.peak <- peakAnno.info.frame[which(peakAnno.info.frame$annotation == "Distal"),] # 52102 peaks
write.table(distal.peak[,1:3], file='4.Peak/bed/distalPeaks.bed', row.names=F, col.names=F, sep='\t', quote=F)

# proximal type
proximal.peak <- peakAnno.info.frame[which(peakAnno.info.frame$annotation != "Distal"),] # 160224 peaks
write.table(proximal.peak[,1:3], file='4.Peak/bed/proximalPeaks.bed', row.names=F, col.names=F, sep='\t', quote=F)

########### each cell type
All.DARs <- readRDS("4.Peak/cellType.all.DARs.rds") # 224807 region
idents <- levels(All.DARs$cellType)
res < sapply(idents, function(x){
  cat("writing ", x, "\n")
  DARs.up <- All.DARs %>% dplyr::filter(cellType == x & avg_log2FC >= 0.25 & p_val_adj < 0.05)
  DARs.down <- All.DARs %>% dplyr::filter(cellType == x & avg_log2FC < -0.25 & p_val_adj < 0.05)

  if(nrow(DARs.up) >= 25){
    DARs.up.info <- peakAnno.info.frame[match(DARs.up$genomicRegion, peakAnno.info.frame$id),]
    # distal
    distal.peak.up <- DARs.up.info[which(DARs.up.info$annotation == "Distal"),]
    write.table(distal.peak.up[,1:3], file = paste0('4.Peak/bed/', x,'_up_distalPeaks.bed'), row.names=F, col.names=F, sep='\t', quote=F)
    cat("Distal Up: ", nrow(distal.peak.up), "\n")

    # proximal
    proximal.peak.up <- DARs.up.info[which(DARs.up.info$annotation != "Distal"),]
    write.table(proximal.peak.up[,1:3], file = paste0('4.Peak/bed/', x,'_up_proximalPeaks.bed'), row.names=F, col.names=F, sep='\t', quote=F)
    cat("proximal Up: ", nrow(proximal.peak.up), "\n")
  }

  if(nrow(DARs.down) >= 25){
    DARs.down.info <- peakAnno.info.frame[match(DARs.down$genomicRegion, peakAnno.info.frame$id),]
    # distal
    distal.peak.down <- DARs.down.info[which(DARs.down.info$annotation == "Distal"),]
    write.table(distal.peak.down[,1:3], file = paste0('4.Peak/bed/', x,'_down_distalPeaks.bed'), row.names=F, col.names=F, sep='\t', quote=F)
    cat("Distal Down: ", nrow(distal.peak.down), "\n")

    # proximal
    proximal.peak.down <- DARs.down.info[which(DARs.down.info$annotation != "Distal"),]
    write.table(proximal.peak.down[,1:3], file = paste0('4.Peak/bed/', x,'_down_proximalPeaks.bed'), row.names=F, col.names=F, sep='\t', quote=F)
    cat("proximal Down: ", nrow(proximal.peak.down), "\n")
  }
  return("Finish!")
})

##rGREAT enrichment analysis
source(file = "/home/longzhilin/Analysis_Code/PathwayEnrichment/Great.enrich.R")

idents <- levels(scATAC.data)
bg <- read.table("4.Peak/bed/allPeaks.bed", sep = "\t", stringsAsFactors = F)

pdf("4.Peak/bed/GREAT.enrich.pdf")
res <- lapply(idents, function(x){

  ## up regulate
  up.distal.bed <- tryCatch({
    read.table(paste0("4.Peak/bed/", x, "_up_distalPeaks.bed"), sep = "\t", stringsAsFactors = F)}, error = function(e){return(NULL)})
  if(!is.null(up.distal.bed)){
    cat("Analysis: ", x, " up distal region...\n")
    up_distalPeaks_enrich <- Great.enrich(gr = up.distal.bed, bg = bg, title = paste0(x, "_up_distalPeaks"), Type = "GO Biological Process")
  }else{
    up_distalPeaks_enrich <- NULL
  }
  up.proximal.bed <- tryCatch({
    read.table(paste0("4.Peak/bed/", x, "_up_proximalPeaks.bed"), sep = "\t", stringsAsFactors = F)}, error = function(e){return(NULL)})
  if(!is.null(up.proximal.bed)){
    cat("Analysis: ", x, " up proximal region...\n")  
    up_proximalPeaks_enrich <- Great.enrich(gr = up.proximal.bed, bg = bg, title = paste0(x, "_up_proximalPeaks"), Type = "GO Biological Process")
  }else{
    up_proximalPeaks_enrich <- NULL
  }

  ## down regulate
  down.distal.bed <- tryCatch({
    read.table(paste0("4.Peak/bed/", x, "_down_distalPeaks.bed"), sep = "\t", stringsAsFactors = F)}, error = function(e){return(NULL)})
  if(!is.null(down.distal.bed)){
    cat("Analysis: ", x, " down distal region...\n")
    down_distalPeaks_enrich <- Great.enrich(gr = down.distal.bed, bg = bg, title = paste0(x, "_down_distalPeaks"), Type = "GO Biological Process")
  }else{
    down_distalPeaks_enrich <- NULL
  }
  down.proximal.bed <- tryCatch({
    read.table(paste0("4.Peak/bed/", x, "_down_proximalPeaks.bed"), sep = "\t", stringsAsFactors = F)}, error = function(e){return(NULL)})
  if(!is.null(down.proximal.bed)){
    cat("Analysis: ", x, " down proximal region...\n")  
    down_proximalPeaks_enrich <- Great.enrich(gr = down.proximal.bed, bg = bg, title = paste0(x, "_down_proximalPeaks"), Type = "GO Biological Process")
  }else{
    down_proximalPeaks_enrich <- NULL
  }
  return(list(up_distalPeaks_enrich = up_distalPeaks_enrich, up_proximalPeaks_enrich = up_proximalPeaks_enrich, 
              down_distalPeaks_enrich = down_distalPeaks_enrich, down_proximalPeaks_enrich = down_proximalPeaks_enrich))
})
dev.off()
names(res) <- idents
saveRDS(res, file = "4.Peak/bed/GREAT.enrichment.rds")
library(openxlsx)
write.xlsx(res$Tumor, file = "4.Peak/bed/GREAT.Tumor.xlsx", sheetName = names(res$Tumor), rowNames = F)

####plot --- tumor.up
tumor.up.proximal <- res$Tumor$up_proximalPeaks_enrich
tumor.up.distal <- res$Tumor$up_distalPeaks_enrich
pdf("4.Peak/bed/GREAT.Tumor.selection.pdf", height = unit(5, "inches"))
p <- ggbarplot(tumor.up.proximal, x = "wrap", y = "Hyper_Fold_Enrichment", fill = "grey", color = "grey", width = 0.4, ylab = "Fold Enrichment", xlab = "", orientation = "horiz", sort.val = c("asc"))+theme(legend.position="none") 
print(p)
p <- ggbarplot(tumor.up.distal, x = "wrap", y = "Hyper_Fold_Enrichment", fill = "grey", color = "grey", width = 0.4, ylab = "Fold Enrichment", xlab = "", orientation = "horiz", sort.val = c("asc"))+theme(legend.position="none") 
print(p)
dev.off()
