#' @description: this script will annotate cicero ccan with the ChIPSeeker database to determine what genomic
# regions are linked by the connections and create circos plots

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(openxlsx)
library(circlize)
library(dplyr)
library(tibble)
library(data.table)
library(Signac)
library(Seurat)
library(stringr)
library(RColorBrewer)
library(cicero)

source(file = "/home/longzhilin/Analysis_Code/Plot_colorPaletters.R")
source(file = "/home/longzhilin/Analysis_Code/SingleCell/user.CoveragePlot.R")
setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC")

# read in cell-type-specific CCAN and filter by coaccess threshold > 0.2
dar_files <- list.files("6.Co-Accessible/ccans", recursive = TRUE, pattern = "ciceroConns.*.csv$", full.names = TRUE)
idents <- gsub(".*ciceroConns.", "", dar_files)
idents <- gsub(".csv", "", idents)
list.dar <- lapply(dar_files, function(file_path) {
  fread(file_path) %>%
    dplyr::filter(coaccess > 0.2) %>%
    dplyr::select(Peak1, Peak2)
})
names(list.dar) <- idents
# count cis-regulatory number of coaccess > 0.2
peak.num <- sapply(list.dar, function(x){return(nrow(x))})
names(peak.num) <- idents

# convert the DAR to GRanges objects to annotate
list.dar.peak1.gr <- lapply(seq(list.dar), function(x) {
  df <- list.dar[[x]]
  gr <- StringToGRanges(df$Peak1, sep = c("_","_"))
  return(gr)
})
names(list.dar.peak1.gr) <- idents

list.dar.peak2.gr <- lapply(seq(list.dar), function(x) {
  df <- list.dar[[x]]
  gr <- StringToGRanges(df$Peak2, sep = c("_","_"))
  return(gr)
})
names(list.dar.peak2.gr) <- idents

# annotate the list of GRanges DAR for each cell type with the peak location
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
list.peak1.Anno <- lapply(list.dar.peak1.gr, annotatePeak, TxDb = txdb,
                        tssRegion = c(-3000, 3000), annoDb="org.Hs.eg.db", verbose = FALSE)
saveRDS(list.peak1.Anno, file = "6.Co-Accessible/ccans/list.peak1.Anno.rds")                        
list.peak1.loc <- lapply(seq(list.peak1.Anno), function(x) {
  loc <- list.peak1.Anno[[x]]@anno$annotation 
  loc <- str_split(loc, pattern = " ", simplify = TRUE)[,1]
  loc <- str_replace(loc, pattern="3'", replacement="3p_UTR")
  loc <- str_replace(loc, pattern="5'", replacement="5p_UTR")
  loc <- str_replace(loc, pattern="Distal", replacement="Intergenic")
  return(loc)
})

list.peak2.Anno <- lapply(list.dar.peak2.gr, annotatePeak, TxDb = txdb,
                          tssRegion = c(-3000, 3000), annoDb="org.Hs.eg.db", verbose = FALSE)
saveRDS(list.peak2.Anno, file = "6.Co-Accessible/ccans/list.peak2.Anno.rds")
list.peak2.loc <- lapply(seq(list.peak2.Anno), function(x) {
  loc <- list.peak2.Anno[[x]]@anno$annotation 
  loc <- str_split(loc, pattern = " ", simplify = TRUE)[,1]
  loc <- str_replace(loc, pattern="3'", replacement="3p_UTR")
  loc <- str_replace(loc, pattern="5'", replacement="5p_UTR")
  loc <- str_replace(loc, pattern="Distal", replacement="Intergenic")
  return(loc)
})

# create df with peak1 and peak2 locations
list.peaks.loc.df <- lapply(seq(idents), function(x) {
  df <- cbind(list.peak1.loc[[x]], list.peak2.loc[[x]]) %>% as.data.frame()
  colnames(df) <- c("Peak1","Peak2")
  counts <- dplyr::count(df, Peak1, Peak2) %>% as.data.frame()
})
names(list.peaks.loc.df) <- idents

# generate circos plots of predicted chromatin-chromatin interactions and save to plots directory
Plot_Cicero_Anno <- function(ident, list.peaks.loc.df) {
  clusterID <- ident
  toplot <- as.data.frame(list.peaks.loc.df[[ident]])
  colnames(toplot) <- c("Peak1","Peak2","n")
  print(clusterID)
  # convert to an adjacency list with a value indicating the number of connections between
  # each of the unique genomic location pairs
  unique_combos <- !duplicated(t(apply(toplot, 1, sort)))
  toplot <- toplot[unique_combos, ]
  toplot <- toplot[order(toplot$n), ]
  # toplot$n <- log2(toplot$n)
  
  # grid.col = c('3p_UTR'="red", '5p_UTR'="black", Distal_Intergenic="blue",
  # Downstream="grey", Exon="purple", Intron="orange", Promoter="green")
  grid.col = brewer.pal(7, "Paired")

  circos.clear()
  pdf(paste0("6.Co-Accessible/ccans/plots/ccans.annotation.",clusterID,".pdf"), width=10, height=5)
  par(mar=c(0.5,0.5,0.5,0.5))
  circos.par(gap.after = 5)
  chordDiagram(toplot,
             grid.col = grid.col,
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(toplot))))))

  # we go back to the first track and customize sector labels
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter,
                CELL_META$ylim[1],
                CELL_META$sector.index,
                facing = "clockwise",
                niceFacing = TRUE,
                adj = c(0, 0.5))
   }, bg.border = NA) # here set bg.border to NA is important
  dev.off()
  }

# save circlize plots to file
lapply(idents, function(ident) {
  Plot_Cicero_Anno(ident, list.peaks.loc.df=list.peaks.loc.df)
  })

##############################identify specific cis-regulatory elements
# focus on the tumor or Macrophage cell 
cell.type <- "Tumor"
#Computing average expression & accessibility
scRNA <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/data.merge.pro.rds")
DefaultAssay(scRNA) <- "RNA"
Idents(scRNA) <- scRNA$cellType_low
scATAC <- readRDS("scATAC.data.pro.rds")
DefaultAssay(scATAC) <- "Peaks"
Idents(scATAC) <- scATAC$AnnotatedcellType
average_exp_all <- AverageExpression(scRNA, assays = "RNA")
average_acc_all <- AverageExpression(scATAC, assays = "Peaks")
# extract the same cell type between scRNA and scATAC
average_exp_all$RNA <- average_exp_all$RNA[, colnames(average_acc_all$Peaks)]

####require one of the peaks overlap a gene's promoter (1kb)
peak1.anno <- list.peak1.Anno[[cell.type]]@anno
peak2.anno <- list.peak2.Anno[[cell.type]]@anno
conns <- fread(paste0("6.Co-Accessible/ccans/ciceroConns.", cell.type,".csv")) # 250302
conns <- conns %>% dplyr::filter(coaccess > 0.2)
conns <- conns %>% mutate(peak1_type = peak1.anno$annotation, peak1_distanceToTSS = peak1.anno$distanceToTSS, peak1_nearestGene = peak1.anno$SYMBOL,
                                      peak2_type = peak2.anno$annotation, peak2_distanceToTSS = peak2.anno$distanceToTSS, peak2_nearestGene = peak2.anno$SYMBOL)
bg <- na.omit(unique(c(conns$peak1_nearestGene, conns$peak2_nearestGene))) # 作为background gene：17567
saveRDS(bg, file = paste0("6.Co-Accessible/bg.", cell.type,".rds")) # coaccess > 0.2

#select the peaks overlap the gene's promoter in conns
cCREs <- conns %>% dplyr::filter(str_detect(peak1_type, pattern = "Promoter \\(<=1kb\\)"))
cCREs <- cCREs[!is.na(cCREs$peak1_nearestGene),]
conns_list <- group_split(cCREs, peak1_nearestGene) #基于peak1_nearestGene分组，通常按字母顺序排序
names(conns_list) <- sort(unique(cCREs$peak1_nearestGene)) 
cTargetGenes <- names(conns_list)[names(conns_list) %in% rownames(average_exp_all$RNA)] # 11453 genes

# loop over all genes to compute correlation between accessibility & expression
df <- data.frame()
corr_list <- lapply(cTargetGenes, function(x){
  cat("calculated gene:", x, "\n")
  # subset by cur_gene
  cur_conns <- conns_list[[x]]
  cur_conns <- cur_conns[!(cur_conns$Peak2 %in% cur_conns$Peak1),]
  cur_conns <- subset(cur_conns, coaccess >= 0.2)

  # skip this gene if there are no co-accessible connections
  if(nrow(cur_conns) == 0){return(data.frame())}

  # get average exp and acc for this gene and peaks that are co-accessible
  cur_conns$Peak2 <- gsub("_", "-", cur_conns$Peak2)
  average_acc <- average_acc_all$Peaks[as.character(cur_conns$Peak2),, drop = F]
  
  average_exp <- average_exp_all$RNA[x,]
  # correlation between expression and accessibility:
  cor_mat <- apply(average_acc, 1, function(x){
    correlation <- cor.test(as.numeric(average_exp), as.numeric(x), method='pearson')
    data.frame("pval"=correlation$p.value, "pcc"=correlation$estimate)
  })
  # collapse individual correlation dfs, add info
  cor_df <- Reduce(rbind, cor_mat)
  cur_conns$pcc <- cor_df$pcc
  cur_conns$pval <- cor_df$pval
  return(cur_conns)
})

# combine into one df and remove incomplete entries:
df <- Reduce(rbind, corr_list)
df <- df[!is.na(df$pcc),]

# compute FDR:
df$FDR <- p.adjust(df$pval, method='fdr')
df <- df %>% dplyr::filter(FDR < 0.05) %>% arrange(desc(pcc))
df$peak2_type <- gsub(" .*", "", df$peak2_type)
df$peak2_type <- gsub("'", "_UTR", df$peak2_type)

# distance between peak and target gene
peak1_ranges <- StringToGRanges(gsub("_", "-", df$Peak1))
peak2_ranges <- StringToGRanges(gsub("_", "-", df$Peak2))
df$distance_bp <- abs(start(peak1_ranges) - start(peak2_ranges))

# write df to file:
write.table(df, file = paste0("6.Co-Accessible/", cell.type,"_peak_gene_correlation.csv"), sep = ",", quote=FALSE, row.names=F)
saveRDS(df, file = paste0("6.Co-Accessible/", cell.type,"_peak_gene_correlation.rds"))