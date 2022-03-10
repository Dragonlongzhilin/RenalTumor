#' @description Merge scATAC data of three patients

setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC")
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86) #---GRCh38 (hg38)
library(ggploT1)
library(patchwork)
set.seed(101)
library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) #50G
library(GenomicRanges)
library(ggpubr)
source(file = "/home/longzhilin/Analysis_Code/Plot_colorPaletters.R")

### Creating a common peak set by merging:
T1.peaks <- read.table(
  file = "/data/raw_data/lzl/RenalTumor-20200713/cellRanger_result/scATAC-V1.2/T1/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
T2.peaks <- read.table(
  file = "/data/raw_data/lzl/RenalTumor-20200713/cellRanger_result/scATAC-V1.2/T2/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
T3.peaks <- read.table(
  file = "/data/raw_data/lzl/RenalTumor-20200713/cellRanger_result/scATAC-V1.2/T3/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)

# convert to genomic ranges
gr.T1 <- makeGRangesFromDataFrame(T1.peaks)
gr.T2 <- makeGRangesFromDataFrame(T2.peaks)
gr.T3 <- makeGRangesFromDataFrame(T3.peaks)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.T1, gr.T2, gr.T3))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20] 

# load metadata
md.T1 <- read.table(
  file = "/data/raw_data/lzl/RenalTumor-20200713/cellRanger_result/scATAC-V1.2/T1/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
md.T2 <- read.table(
  file = "/data/raw_data/lzl/RenalTumor-20200713/cellRanger_result/scATAC-V1.2/T2/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
md.T3 <- read.table(
  file = "/data/raw_data/lzl/RenalTumor-20200713/cellRanger_result/scATAC-V1.2/T3/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

# perform an initial filtering of low count cells
md.T1 <- md.T1[md.T1$passed_filters > 500, ] # 59599*17
md.T2 <- md.T2[md.T2$passed_filters > 500, ] # 43387*17
md.T3 <- md.T3[md.T3$passed_filters > 500, ] # 43162*17

# create fragment objects
frags.T1 <- CreateFragmentObject(
  path = "/data/raw_data/lzl/RenalTumor-20200713/cellRanger_result/scATAC-V1.2/T1/outs/fragments.tsv.gz",
  cells = rownames(md.T1)
)
frags.T2 <- CreateFragmentObject(
  path = "/data/raw_data/lzl/RenalTumor-20200713/cellRanger_result/scATAC-V1.2/T2/outs/fragments.tsv.gz",
  cells = rownames(md.T2)
)
frags.T3 <- CreateFragmentObject(
  path = "/data/raw_data/lzl/RenalTumor-20200713/cellRanger_result/scATAC-V1.2/T3/outs/fragments.tsv.gz",
  cells = rownames(md.T3)
)

# Quantify peaks in each dataset
T1.counts <- FeatureMatrix(
  fragments = frags.T1,
  features = combined.peaks,
  cells = rownames(md.T1)
)
T2.counts <- FeatureMatrix(
  fragments = frags.T2,
  features = combined.peaks,
  cells = rownames(md.T2)
)
T3.counts <- FeatureMatrix(
  fragments = frags.T3,
  features = combined.peaks,
  cells = rownames(md.T3)
)

# Create the objects
T1_assay <- CreateChromatinAssay(T1.counts, fragments = frags.T1)
T1 <- CreateSeuratObject(T1_assay, assay = "ATAC", meta.data = md.T1)
T2_assay <- CreateChromatinAssay(T2.counts, fragments = frags.T2)
T2 <- CreateSeuratObject(T2_assay, assay = "ATAC", meta.data = md.T2)
T3_assay <- CreateChromatinAssay(T3.counts, fragments = frags.T3)
T3 <- CreateSeuratObject(T3_assay, assay = "ATAC", meta.data = md.T3)

# Merge object
# add information to identify dataset of origin
T1$dataset <- 'T1'
T2$dataset <- 'T2'
T3$dataset <- 'T3'

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = T1,
  y = list(T2, T3),
  add.cell.ids = c("T1", "T2", "T3")
)
combined[["ATAC"]]

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
# add the gene information to the object
Annotation(combined) <- annotations

# Quality Control
# compute nucleosome signal score per cell
combined <- NucleosomeSignal(object = combined)
# compute TSS enrichment score per cell
combined <- TSSEnrichment(object = combined, fast = FALSE, assay = 'ATAC', verbose = T)
# add blacklist ratio and fraction of reads in peaks
combined$pct_reads_in_peaks <- combined$peak_region_fragments / combined$passed_filters * 100
combined$blacklist_ratio <- combined$blacklist_region_fragments / combined$peak_region_fragments
combined$high.tss <- ifelse(combined$TSS.enrichment > 2, 'High', 'Low')
combined$nucleosome_group <- ifelse(combined$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
pdf("1.QualityControl/merge.QC.pdf")
TSSPlot(combined, group.by = 'high.tss') + NoLegend()
##Fragment length periodicity:
FragmentHistogram(object = combined, group.by = 'nucleosome_group')
VlnPlot(
  object = combined,
  features = c('pct_reads_in_peaks', 'peak_region_fragments','TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0,
  ncol = 5
)
VlnPlot(
  object = combined,
  features = c('pct_reads_in_peaks'),
  pt.size = 0,
  y.max = c(20)
)
VlnPlot(
  object = combined,
  features = c('peak_region_fragments'),
  pt.size = 0,
  y.max = c(2500)
)
VlnPlot(
  object = combined,
  features = c('TSS.enrichment'),
  pt.size = 0,
  y.max = c(10)
)
VlnPlot(
  object = combined,
  features = c('blacklist_ratio'),
  pt.size = 0,
  y.max = c(0.05)
)
VlnPlot(
  object = combined,
  features = c('nucleosome_signal'),
  pt.size = 0,
  y.max = c(15)
)
dev.off()
saveRDS(combined, file = "scATAC.merge.rds") # 146148

## remove outlier cells based on QC metrics:
combined.pro <- subset(
  x = combined,
  subset = peak_region_fragments > 1000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3
)
# 88392*24173

# process orig.ident information
combined.pro@meta.data$orig.ident <- gsub("_.*", "", rownames(combined.pro@meta.data))
saveRDS(combined.pro, file = "scATAC.merge.pro.rds")

##---- plot sample distribution
cell.number <- as.data.frame(table(combined.pro@meta.data$orig.ident))
pdf("1.QualityControl/highQuality.pdf")
ggbarplot(cell.number, x="Var1", y="Freq", fill = "Var1", color = "Var1", palette = "npg",
          sort.by.groups=FALSE,
          label = T, xlab = "", ylab = "Cell Number") + theme(legend.position="none") 
dev.off()