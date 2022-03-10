#' @description: regulated network in CD8+ T cell

library(Signac)
library(Seurat)
library(harmony)
library(SeuratWrappers)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86) #---GRCh38 (hg38)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(clustree)
library(vegan)
library(reshape2)
library(openxlsx)
library(tidyverse)
library(data.table)
set.seed(101)
library(future)
plan("multiprocess", workers = 5)
options(future.globals.maxSize = 50000 * 1024^2) #50G
library(GenomicRanges)
library(ggpubr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVARmotifs)
data("human_pwms_v2")
source(file = "/home/longzhilin/Analysis_Code/Plot_colorPaletters.R")
source("/home/longzhilin/Analysis_Code/SingleCell/scATAC.Integrate.multipleSample.R")
source(file = "/home/longzhilin/Analysis_Code/SingleCell/Integrate.scRNA.scATAC.R")
setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC")
source(file = "/home/longzhilin/Analysis_Code/code/ratio.plot.R")

sub.scATAC.harmony <- readRDS("8.Immune/CD8T/sub.scATAC.harmony.pro.rds")
sub.scATAC.harmony$cellType <- sub.scATAC.harmony$cellType3
Idents(sub.scATAC.harmony) <- sub.scATAC.harmony$cellType
cellType.colors <- c("#F8766D", "#00BA38", "#619CFF")
DefaultAssay(sub.scATAC.harmony) <- "Peaks"

sub.scRNA.harmony <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/5.Immune/CD8T/sub.scRNA.harmony.pro.rds")
sub.scRNA.harmony$cellType <- sub.scRNA.harmony$cellType3
Idents(sub.scRNA.harmony) <- sub.scRNA.harmony$cellType
DefaultAssay(sub.scRNA.harmony) <- "RNA"

####--------------------------------------------------------------------------- 5. TF regulated network
library(ggseqlogo)
interest.genes <- "CD27"
celltype <- "Exhaustion"

# motif information
motif.info <- data.frame(originName = names(sub.scATAC.harmony@assays$Peaks@motifs@motif.names), TF = unlist(sub.scATAC.harmony@assays$Peaks@motifs@motif.names))
rownames(motif.info) <- NULL
# subset by celltype
scATAC.sub <- subset(sub.scATAC.harmony, cellType == celltype) # 376
scRNA.sub <- subset(sub.scRNA.harmony, cellType == celltype)  # 1803

scATAC.sub <- FindTopFeatures(scATAC.sub, min.cutoff=ncol(scATAC.sub)*0.05) # present in what % of cells
scRNA.sub <- FindVariableFeatures(scRNA.sub, nfeatures = 3000)

# use chipseeker annotate each peak in scATAC-seq
peak.info <- readRDS("4.Peak/peak.annotation.simple.ChIPseeker.rds")

PWMatrixToProbMatrix <- function(x){
        if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
        m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
        m <- t(t(m)/colSums(m))
        m
}

################################################################################
# Find promoters & genes with accessible TF binding sites
################################################################################
# settings for featureplot
order_values <- TRUE
reduction <- 'umap'

res <- sapply(interest.genes, function(x){
    motif_name <- x
    motif_ID <- motif.info$originName[match(x, motif.info$TF)] # e.g., ENSG00000008196_LINE2_TFAP2B_D_N1

    # get all regions with motif binding site:
    DefaultAssay(scATAC.sub) <- "Peaks"
    motif_accessible <- Motifs(scATAC.sub)@data[, motif_ID]
    motif_accessible <- names(motif_accessible)[motif_accessible > 0]

    # subset this list by top features
    motif_accessible <- motif_accessible[motif_accessible %in% VariableFeatures(scATAC.sub)]

    # which of these peaks are at promoters?
    motif_accessible_promoters <- motif_accessible[motif_accessible %in% peak.info$peaks[which(peak.info$originType == "Promoter (<=1kb)")]]
    # which genes are associated with these promoters?
    motif_target_genes <- peak.info$SYMBOL[match(motif_accessible_promoters, peak.info$peaks)]

    # optional:
    # which of these genes are highly expressed in snRNA-seq?
    motif_target_genes.RNA <- motif_target_genes[motif_target_genes %in% VariableFeatures(scRNA.sub)] %>% as.character %>% unique

    # remove genes that are not in the seurat obj
    motif_target_genes <- as.character(motif_target_genes)
    motif_target_genes <- motif_target_genes[motif_target_genes %in% rownames(scRNA.sub)]

    # output 
    target.genelists <- list(data.frame(motif_target_genes.RNA), data.frame(motif_target_genes))
    write.xlsx(target.genelists, file = paste0("8.Immune/CD8T/TF.analysis/", x, "_targetGenes.xlsx"), sheetName = c("overlap with hvgs in scRNA-seq", "all target genes"))
    
    ################################################################################
    # Compute module score for these target genes
    ################################################################################
    gene_list <- list(motif_target_genes)
    names(gene_list) <- paste0(motif_name, '_targets')

    # add in whole data
    sub.scRNA.harmony <- AddModuleScore(
        sub.scRNA.harmony,
        features=gene_list,
        pool = rownames(sub.scRNA.harmony),
        name=paste0(motif_name, '_targets')
    )

    ################################################################################
    # plot module score feature plot:
    ################################################################################

    # plot promoter target gene module score for this TF:
    p1 <- FeaturePlot(sub.scRNA.harmony, features = paste0(motif_name, '_targets1'), order = order_values, reduction = reduction) +
    scale_color_gradient2(low=scales::muted('blue'), mid='white', high=scales::muted('red')) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "in"))  +
    ggtitle(paste0(motif_name, ' target score'))
    pdf(paste0('8.Immune/CD8T/TF.analysis/', celltype, '_', motif_name, '_targets.pdf'), width=5, height=4, useDingbats=FALSE)
    print(p1)

    # plot chromVAR deviation for this TF:
    DefaultAssay(sub.scATAC.harmony) <- "chromvar"
    p2 <- FeaturePlot(sub.scATAC.harmony, features = gsub("_", "-", motif_ID), order = order_values, reduction = reduction) +
    scale_color_gradient2(
        low=rgb(32, 67, 37, maxColorValue=255), mid='white', high=rgb(58, 22, 72, maxColorValue=255)) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "in")) +
        ggtitle(paste0(motif_name, ' motif'))
    print(p2)

    ################################################################################
    # cluster violin plot for target expression modules
    ################################################################################
    my_comparisons <- as.list(as.data.frame(combn(levels(sub.scRNA.harmony$cellType),2)))
    a <- sub.scRNA.harmony@meta.data[, c(paste0(motif_name, '_targets1'), "cellType")]
    colnames(a) <- c("Signature score", "group")
    p3 <- ggboxplot(a, x = "group", y = "Signature score", title = x, color = "black", fill = "group", add = "jitter", add.params = list(color = "black", size = 0.1)) + 
         xlab("")  + ylab(paste0(motif_name, ' targets')) + rotate_x_text(angle = 45, vjust = 1) + scale_fill_manual(values = cellType.colors) + 
         stat_compare_means(comparisons = my_comparisons, label = "p.signif") + NoLegend()
    print(p3)

    DefaultAssay(sub.scATAC.harmony) <- "chromvar"
    a <- FetchData(sub.scATAC.harmony, vars = gsub("_", "-", motif_ID), slot = "data")
    a$group <- sub.scATAC.harmony$cellType
    my_comparisons <- as.list(as.data.frame(combn(levels(sub.scATAC.harmony$cellType),2)))
    colnames(a) <- c("Signature score", "group")
    p4 <- ggboxplot(a, x = "group", y = "Signature score", title = x, color = "black", fill = "group", add = "jitter", add.params = list(color = "black", size = 0.1)) + 
         xlab("")  + ylab(paste0(motif_name, ' deviation')) + rotate_x_text(angle = 45, vjust = 1) + scale_fill_manual(values = cellType.colors) + 
         stat_compare_means(comparisons = my_comparisons, label = "p.signif") + NoLegend()
    print(p4)

    # motif logo
    PPM <- PWMatrixToProbMatrix(human_pwms_v2[[motif_ID]])
    p5 <- ggseqlogo(PPM) + ggtitle(motif_name) + theme(plot.title = element_text(hjust = 0.5))
    print(p5)    
    dev.off()
})

################################################################################
# co-accessibility
################################################################################
library(cicero)
library(monocle3)
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
idents <- levels(sub.scATAC.harmony)
# 不分析mast cell，由于细胞数目小于100
list.ccan <- lapply(idents, function(ident) {Get_Ccans(ident, seurat_atac_obj = sub.scATAC.harmony)})

# write to file
names(list.ccan) <- idents
lapply(names(list.ccan), function(x) {
  fwrite(list.ccan[[x]]$df, file = paste0("8.Immune/CD8T/Co-Accessible/ciceroConns.",x,".csv"), row.names = F)
})
saveRDS(list.ccan, file = "8.Immune/CD8T/Co-Accessible/cellType.ccans.rds")

# calculate a global CCAN for all celltypes
ccan <- Get_Ccans(seurat_atac_obj = sub.scATAC.harmony)
fwrite(ccan$df, file = "8.Immune/CD8T/Co-Accessible/ciceroConns.allcells.csv", row.names = F)
saveRDS(ccan, file = "8.Immune/CD8T/Co-Accessible/ccans/All.ccans.rds")

#### 提取感兴趣的ccan
ciceroCds <- list.ccan[["celltype"]]$ciceroCds
conns <- list.ccan[["celltype"]]$conns
ccans <- list.ccan[["celltype"]]$CCAN_assigns
saveRDS(ciceroCds, file = paste0("8.Immune/CD8T/Co-Accessible/", celltype,".ciceroCds.rds"))
saveRDS(conns, file = paste0("8.Immune/CD8T/Co-Accessible/", celltype,".conns.rds"))
saveRDS(ccans, file = paste0("8.Immune/CD8T/Co-Accessible/", celltype,".ccans.rds"))
fwrite(list.ccan[[celltype]]$df, file = paste0("8.Immune/CD8T/Co-Accessible/ciceroConns.", celltype, ".csv"), row.names = F)

################################################################################
# annotated ccans
################################################################################
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# read in cell-type-specific CCAN and filter by coaccess threshold > 0.2
dar_files <- list.files("8.Immune/CD8T/Co-Accessible", recursive = TRUE, pattern = "ciceroConns.*.csv$", full.names = TRUE)
idents <- gsub(".*ciceroConns.", "", dar_files)
idents <- gsub(".csv", "", idents)
list.dar <- lapply(dar_files, function(file_path) {
  fread(file_path) %>%
    dplyr::filter(coaccess > 0.2) %>%
    dplyr::select(Peak1, Peak2)
})
names(list.dar) <- idents

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
save(list.peak1.Anno, list.peak2.Anno, file = "8.Immune/CD8T/Co-Accessible/ciceroConns.annotation.RData")

################################################################################
# identify cis-regulated elements
################################################################################
####one of the peaks overlap a gene's promoter (1kb)
peak1.anno <- list.peak1.Anno[[celltype]]@anno
peak2.anno <- list.peak2.Anno[[celltype]]@anno
conns <- fread(paste0("8.Immune/CD8T/Co-Accessible/ciceroConns.", celltype, ".csv"))
conns <- conns %>% dplyr::filter(coaccess > 0.2) # 249226 peak-pairs
conns <- conns %>% mutate(peak1_type = peak1.anno$annotation, peak1_distanceToTSS = peak1.anno$distanceToTSS, peak1_nearestGene = peak1.anno$SYMBOL,
                                      peak2_type = peak2.anno$annotation, peak2_distanceToTSS = peak2.anno$distanceToTSS, peak2_nearestGene = peak2.anno$SYMBOL)
bg <- na.omit(unique(c(conns$peak1_nearestGene, conns$peak2_nearestGene))) # 作为background gene：17567
saveRDS(bg, file = paste0("8.Immune/CD8T/Co-Accessible/bg.", celltype,".rds")) # coaccess > 0.2

#筛选出包含promoter的conns
cCREs <- conns %>% dplyr::filter(str_detect(peak1_type, pattern = "Promoter \\(<=1kb\\)"))
cCREs <- cCREs[!is.na(cCREs$peak1_nearestGene),]
conns_list <- group_split(cCREs, peak1_nearestGene) #基于peak1_nearestGene分组，通常按字母顺序排序
names(conns_list) <- sort(unique(cCREs$peak1_nearestGene))

sub.scRNA <- subset(sub.scRNA.harmony, subset = cellType == "Exhausted")
source(file = "/home/longzhilin/Analysis_Code/code/Filter.gene.R")
sub.scRNA <- Filter.gene(sub.scRNA)
cTargetGenes <- names(conns_list)[names(conns_list) %in% rownames(sub.scRNA)] # 11453 genes

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
  return(cur_conns)
})
# combine into one df and remove incomplete entries:
df <- Reduce(rbind, corr_list)
df$peak2_type <- gsub(" .*", "", df$peak2_type)
df$peak2_type <- gsub("'", "_UTR", df$peak2_type)

# distance between peak and target gene
peak1_ranges <- StringToGRanges(gsub("_", "-", df$Peak1))
peak2_ranges <- StringToGRanges(gsub("_", "-", df$Peak2))
df$distance_bp <- abs(start(peak1_ranges) - start(peak2_ranges))
# write df to file:
write.table(df, file = paste0("8.Immune/CD8T/Co-Accessible/", celltype,"_peak_gene_correlation.csv"), sep = ",", quote=FALSE, row.names=F)
saveRDS(df, file = paste0("8.Immune/CD8T/Co-Accessible/", celltype,"_peak_gene_correlation.rds"))

#### defined the enhancer
cCREs <- df
#gene-cCRE links
cCREs <- cCREs %>% subset(distance_bp/1000>=10 & peak2_type != 'Promoter')
cCREs$target_gene <- cCREs$peak1_nearestGene
cCREs$peak_gene <- paste0(cCREs$Peak2, '_', cCREs$peak1_nearestGene)

################################################################################
# Construct TF nets
################################################################################
library(igraph)
library(RColorBrewer)

TFs.info <- motif.info[match(interest.genes, motif.info$TF),]

# for each motif, find genes:
use_variable_genes <- T
motif_list <- list()
edge_df.list <- data.frame()
vertex_df.list <- data.frame()
for(i in 1:nrow(TFs.info)){
  motif_name <- TFs.info[i,2]
  motif_ID <- TFs.info[i,1]

  # get list of promoter and enhancer targets of these TFs
  DefaultAssay(scATAC.sub) <- "Peaks"
  motif_accessible <- Motifs(scATAC.sub)@data[,motif_ID]
  motif_accessible <- names(motif_accessible)[motif_accessible > 0]
  motif_accessible <- motif_accessible[motif_accessible %in% VariableFeatures(scATAC.sub)]  
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

################################################################################
# visual settings for network
################################################################################
# add DEGs
sub.scRNA.DEGs <- readRDS("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/5.Immune/CD8T/cellType.DE.pro.rds")
sub.scRNA.DEGs <- filter(sub.scRNA.DEGs, abs(avg_log2FC) > 0.25 & p_val_adj < 0.05)
cellType.DEGs <- sub.scRNA.DEGs[which(sub.scRNA.DEGs$cluster == celltype),]

# color vertices based on Diagnosis DEGs:
up_genes <- cellType.DEGs %>% filter(p_val_adj < 0.05 & avg_log2FC >= 0.5) %>% .$gene
down_genes <- cellType.DEGs %>% filter(p_val_adj < 0.05 & avg_log2FC < -0.5) %>% .$gene

# remove labels if gene is not DE, or not a TF:
de_targets <- as.character(vertex_df.list$name[vertex_df.list$name %in% unique(c(up_genes, down_genes))])
vertex_df.list$label <- ifelse(vertex_df.list$name %in% de_targets, vertex_df.list$name, '')
vertex_df.list$label <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), vertex_df.list$name, vertex_df.list$label)
#vertex_df.list$label <- ifelse(vertex_df.list$name %in% gwas_genes, vertex_df.list$name, vertex_df.list$label)

# set node color based on control vs AD DEGs:
vertex_df.list$color <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), '#1E90FF', "#FFFFFF")
vertex_df.list$color <- ifelse(vertex_df.list$name %in% up_genes, "#E87D72",  "#FFFFFF")
vertex_df.list$color <- ifelse(vertex_df.list$name %in% down_genes, '#55BCC2', vertex_df.list$color)
vertex_df.list$color <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), '#1E90FF',vertex_df.list$color)
#vertex_df.list$color <- ifelse(vertex_df.list$name %in% gwas_genes, 'orange', vertex_df.list$color)

# italics font for genes:
vertex_df.list$font <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), 2, 1)

# set size to larger if the gene is a TF:
vertex_df.list$size <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), 10, 2)

other_tfs <- setdiff(motif.info$TF, interest.genes)

#vertex_df.list$size <- ifelse((vertex_df.list$name %in% de_targets | vertex_df.list$name %in% gwas_genes | vertex_df.list$name %in% other_tfs), 5, 2)
vertex_df.list$size <- ifelse((vertex_df.list$name %in% de_targets | vertex_df.list$name %in% other_tfs), 5, 2)
vertex_df.list$size <- ifelse(vertex_df.list$name %in% interest.genes, 10, vertex_df.list$size)

################################################################################
# graph all nodes
################################################################################
enhancer_color <- '#FFC125'
promoter_color <- '#00CED1'

# repressor_color <- 'gray'

g <- igraph::graph_from_data_frame(edge_df.list, directed=TRUE, vertices = vertex_df.list)
l <- layout_with_fr(g)

edge_colors <- ifelse(E(g)$type == 'promoter', promoter_color, enhancer_color)
# edge_colors <- ifelse(E(g)$type == 'repressor', repressor_color, edge_colors)

pdf(paste0('8.Immune/CD8T/TF.analysis/', celltype, '_TF_interaction_graph.pdf'), width=10, height=10, useDingbats=FALSE)
plot(
  g, layout=l,
  vertex.size=vertex_df.list$size,
  edge.color=edge_colors,
  edge.alpha=0.5,
  vertex.color=vertex_df.list$color,
  vertex.label=vertex_df.list$label, vertex.label.family='Helvetica', vertex.label.font=vertex_df.list$font,
  vertex.label.color = 'black',
  vertex.frame.color = "grey",
  vertex.alpha = 0.8,
  edge.arrow.size=0.25
)
dev.off()
