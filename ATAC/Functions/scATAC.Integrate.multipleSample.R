#' @param scATAC.object: scATAC object after RunSVD function

# The first PC is usually highly correlated with the sequencing depth, so it generally needs to be removed
ElbowPlot.user <- function(object, ndims = 20, reduction = "pca"){
    data.use <- Stdev(object = object, reduction = reduction)
    if (length(x = data.use) == 0) {
        stop(paste("No standard deviation info stored for", reduction))
    }
    if (ndims > length(x = data.use)) {
        warning("The object only has information for ", length(x = data.use),
            " reductions")
        ndims <- length(x = data.use)
    }
    stdev <- "Standard Deviation"
    plot <- ggplot(data = data.frame(dims = 2:ndims, stdev = data.use[2:ndims])) +
        geom_point(mapping = aes_string(x = "dims", y = "stdev")) +
        labs(x = gsub(pattern = "_$", replacement = "", x = Key(object = object[[reduction]])),
            y = stdev)
    return(plot)
}

####Seurat Integration
Seurat.integration.reduceDimension <- function(scATAC.object, set.resolutions, groups = "dataset", PC = 20, assay = "ATAC", npcs = 30){
    DefaultAssay(scATAC.object) <- assay
    scATAC.object <- RunSVD(scATAC.object, n = npcs)
    p <- DepthCor(scATAC.object)
    print(p)
    p <- ElbowPlot.user(object = scATAC.object, ndims = npcs, reduction = "lsi")
    print(p)    

    scATAC.lists <- SplitObject(object = scATAC.object, split.by = groups)
    integration.anchors <- FindIntegrationAnchors(
      object.list = scATAC.lists,
      anchor.features = rownames(scATAC.lists[[1]]), #common feature
      reduction = "rlsi",
      dims = 2:PC
    )
    # integrate LSI embeddings
    seurat.integrated.scATAC <- IntegrateEmbeddings(
        anchorset = integration.anchors,
        reductions = scATAC.object[["lsi"]],
        new.reduction.name = "integrated_lsi",
        dims.to.integrate = 1:PC
    )
    # create a new UMAP using the integrated embeddings
    seurat.integrated.scATAC <- RunUMAP(seurat.integrated.scATAC, reduction = "integrated_lsi", dims = 2:PC)
    seurat.integrated.scATAC <- RunTSNE(seurat.integrated.scATAC, reduction = "integrated_lsi", dims = 2:PC)
    p <- DimPlot(seurat.integrated.scATAC, group.by = groups, pt.size = 0.1, reduction = 'umap')
    print(p)
    p <- DimPlot(seurat.integrated.scATAC, group.by = groups, pt.size = 0.1, reduction = 'tsne')
    print(p)

    ##聚类分析
    seurat.integrated.scATAC <- FindNeighbors(object = seurat.integrated.scATAC, reduction = 'integrated_lsi', dims = 2:PC)
    seurat.integrated.scATAC <- FindClusters(object = seurat.integrated.scATAC, resolution = set.resolutions, verbose = FALSE, algorithm = 3)
    p <- clustree(seurat.integrated.scATAC)
    print(p)
    merge.res <- sapply(set.resolutions, function(x){
        p <- DimPlot(object = seurat.integrated.scATAC, reduction = 'umap',label = TRUE, group.by = paste0(assay, "_snn_res.", x)) + NoLegend()
        print(p)
        p <- DimPlot(object = seurat.integrated.scATAC, reduction = 'tsne',label = TRUE, group.by = paste0(assay, "_snn_res.", x)) + NoLegend()
        print(p)
    })
    return(seurat.integrated.scATAC)
}

####Harmony
Harmony.integration.reduceDimension <- function(scATAC.object, set.resolutions, groups = "dataset", PC = 20, assay = "ATAC", npcs = 30){
    DefaultAssay(scATAC.object) <- assay
    scATAC.object <- RunSVD(scATAC.object, n = npcs)
    p <- DepthCor(scATAC.object)
    print(p)
    p <- ElbowPlot.user(object = scATAC.object, ndims = npcs, reduction = "lsi")
    print(p)

    Harmony.scATAC <- RunHarmony(object = scATAC.object, group.by.vars = groups, assay.use = assay, reduction = "lsi", project.dim = FALSE, verbose = FALSE)
    Harmony.scATAC <- RunUMAP(Harmony.scATAC, reduction = "harmony", dims = 2:PC, verbose = FALSE)
    Harmony.scATAC <- RunTSNE(Harmony.scATAC, reduction = "harmony", dims = 2:PC, verbose = FALSE)

    p <- DimPlot(Harmony.scATAC, group.by = groups, pt.size = 0.1, reduction = 'umap')
    print(p)
    p <- DimPlot(Harmony.scATAC, group.by = groups, pt.size = 0.1, reduction = 'tsne')
    print(p)

    ##聚类分析
    Harmony.scATAC <- FindNeighbors(object = Harmony.scATAC, reduction = 'harmony', dims = 2:PC)
    Harmony.scATAC <- FindClusters(object = Harmony.scATAC, resolution = set.resolutions, verbose = FALSE, algorithm = 3)

    p <- clustree(Harmony.scATAC)
    print(p)
    merge.res <- sapply(set.resolutions, function(x){
        p <- DimPlot(object = Harmony.scATAC, reduction = 'umap',label = TRUE, group.by = paste0(assay, "_snn_res.", x)) + NoLegend()
        print(p)
        p <- DimPlot(object = Harmony.scATAC, reduction = 'tsne',label = TRUE, group.by = paste0(assay, "_snn_res.", x)) + NoLegend()
        print(p)
    })
    return(Harmony.scATAC)
}