#' @description: Projecting scRNAseq data onto a reference map of Tumour-Infiltrating Lymphocytes
# ref:Interpretation of T cell states from single-cell transcriptomics data using reference atlases


library(ProjecTILs)
ref <- load.reference.map("/data/ExtraDisk/sdd/longzhilin/Data/ref_LCMV_Atlas_mouse_v1.rds")
seurat.object <- sub.scRNA.harmony

#### Let’s explore the reference atlas
refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000", "#87f6a5", "#e812dd")
DimPlot(ref,label = T, cols = refCols)

markers <- c("Cd4","Cd8a","Ccr7","Tcf7","Pdcd1","Havcr2","Tox","Izumo1r","Cxcr6","Xcl1","Gzmb","Gzmk","Ifng","Foxp3")
VlnPlot(ref,features=markers,stack = T,flip = T,assay = "RNA")

#### Run Projection algorithm
query.projected <- make.projection(querydata, ref=ref)
plot.projection(ref, query.projected)

# each batch
data.seurat.list <- SplitObject(seurat.object, split.by = "orig.ident")
query.projected.list <- make.projection(data.seurat.list, ref=ref, ncores=10, future.maxSize=10000, filter.cells = F)
pdf("5.Immune/CD8T/ProjecTILs.LCMVCD8.pdf")
p1 <- plot.projection(ref, query.projected.list$T2) + ggtitle("T2")
p2 <- plot.projection(ref, query.projected.list$T3) + ggtitle("T3")
p3 <- plot.projection(ref, query.projected.list$T4) + ggtitle("T4")
p4 <- plot.projection(ref, query.projected.list$T5) + ggtitle("T5")
ggarrange(p1,p2,p3,p4, ncol=2, nrow =2)
dev.off()

# Predict cell states
query.projected <- cellstate.predict(ref=ref, query=query.projected)
table(query.projected$functional.cluster)
plot.statepred.composition(ref, query.projected,metric = "Percent")
plot.states.radar(ref, query=query.projected, min.cells=30)

# merge the functional cluster
T2 <- cellstate.predict(ref=ref, query=query.projected.list$T2)
T3 <- cellstate.predict(ref=ref, query=query.projected.list$T3)
T4 <- cellstate.predict(ref=ref, query=query.projected.list$T4)
T5 <- cellstate.predict(ref=ref, query=query.projected.list$T5)

T2.info <- data.frame(cell.id = colnames(T2), functional.cluster = T2$functional.cluster)
T3.info <- data.frame(cell.id = colnames(T3), functional.cluster = T3$functional.cluster)
T4.info <- data.frame(cell.id = colnames(T4), functional.cluster = T4$functional.cluster)
T5.info <- data.frame(cell.id = colnames(T5), functional.cluster = T5$functional.cluster)

anno.info <- rbind(T2.info,rbind(T3.info, T4.info))
idx <- match(colnames(seurat.object), anno.info$cell.id)
seurat.object$functional.cluster <- anno.info$functional.cluster[idx]
DimPlot(seurat.object, group.by = "functional.cluster")