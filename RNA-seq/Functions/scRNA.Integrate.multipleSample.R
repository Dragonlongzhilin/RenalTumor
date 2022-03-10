#' @author longzhilin
#' @decscription using harmony method to integrate mulitiple sample and adjust patient bias
#' @param seurat.object: the merged seurat object after quality control (NormalizeData, FindVariableFeatures, ScaleData)
#' @param set.resolutions: the resolutions for FindClusters function
#' @param assay: RNA or SCT
#' @param PC: PCA numbers to reduce dimension and cluster
#' @param npcs: RunPCA parameters and ElbowPlot show the all PCs

####SCT or RNA
####Harmony correction batch
Harmony.integration.reduceDimension <- function(seurat.object, set.resolutions, assay = "RNA", nfeatures = 3000, PC = 50, npcs = 100){
     require(Seurat)
     require(harmony)
     require(clustree)
     require(dplyr)
     #OK nfeatures
     DefaultAssay(seurat.object) <- assay

     seurat.object <- RunPCA(seurat.object, verbose = FALSE, npcs = npcs)
     p <- ElbowPlot(object = seurat.object, ndims = npcs)
     print(p)
     #If dims.use is not specified, all PCs will be used by default
     #assay.use defaults to RNA. If the SCTransform standardization method is used, you need to specify assay.use="SCT" in RunHarmony
     #Compare the differences between the two modes and find that the default is better and more complete
    seurat.object <- RunHarmony(object = seurat.object, group.by.vars = "orig.ident", assay.use=assay, verbose = FALSE)
    seurat.object <- RunUMAP(seurat.object, reduction = "harmony", dims = 1:PC, verbose = FALSE)
    seurat.object <- RunTSNE(seurat.object, reduction = "harmony", dims = 1:PC, verbose = FALSE)
    seurat.object <- FindNeighbors(seurat.object, dims = 1:PC, reduction = "harmony", verbose = FALSE) #使用harmony替代PCA
    seurat.object <- FindClusters(seurat.object, resolution = set.resolutions, verbose = FALSE)
    p <- clustree(seurat.object)
    print(p)
    p <- DimPlot(object = seurat.object, reduction = 'umap',label = TRUE, group.by = "orig.ident")
    print(p)
    merge.res <- sapply(set.resolutions, function(x){
        p <- DimPlot(object = seurat.object, reduction = 'umap',label = TRUE, group.by = paste0(assay, "_snn_res.", x)) + NoLegend()
        print(p)
    })
    p <- DimPlot(object = seurat.object, reduction = 'tsne',label = TRUE, group.by = "orig.ident")
    print(p)
    merge.res <- sapply(set.resolutions, function(x){
        p <- DimPlot(object = seurat.object, reduction = 'tsne',label = TRUE, group.by = paste0(assay, "_snn_res.", x)) + NoLegend()
        print(p)
    })
    return(seurat.object)
}

####Seurat's own integration method
Seurat.integration.reduceDimension <- function(seurat.lists, set.resolutions, assay = "RNA", PC = 50, nfeatures = 3000, npcs = 100){
     require(Seurat)
     require(clustree)
     require(dplyr)
     #Determine natures and assays
    seurat.features <- SelectIntegrationFeatures(object.list = seurat.lists, nfeatures = nfeatures, verbose = FALSE)

    if(assay == "SCT"){
        seurat.lists <- PrepSCTIntegration(object.list = seurat.lists, anchor.features = seurat.features, verbose = FALSE)
        seurat.anchors <- FindIntegrationAnchors(object.list = seurat.lists, normalization.method = "SCT", anchor.features = seurat.features, verbose = FALSE)
        data.merge.seurat <- IntegrateData(anchorset = seurat.anchors, normalization.method = "SCT", verbose = FALSE)           
    }else{
        seurat.anchors <- FindIntegrationAnchors(object.list = seurat.lists, anchor.features = seurat.features, verbose = FALSE)
        data.merge.seurat <- IntegrateData(anchorset = seurat.anchors, verbose = FALSE)        
    }
    DefaultAssay(data.merge.seurat) <- "integrated"
    data.merge.seurat <- ScaleData(data.merge.seurat, verbose = FALSE)
    data.merge.seurat <- RunPCA(data.merge.seurat, verbose = FALSE, npcs = npcs)
    p <- ElbowPlot(object = data.merge.seurat, ndims = npcs)
    print(p)    
    data.merge.seurat <- RunUMAP(data.merge.seurat, reduction = "pca", dims = 1:PC, verbose = FALSE)
    data.merge.seurat <- RunTSNE(data.merge.seurat, reduction = "pca", dims = 1:PC, verbose = FALSE)
    data.merge.seurat <- FindNeighbors(data.merge.seurat, reduction = "pca", dims = 1:PC, verbose = FALSE)
    data.merge.seurat <- FindClusters(data.merge.seurat, resolution = set.resolutions, verbose = FALSE)
    p <- clustree(data.merge.seurat)
    print(p)
    p <- DimPlot(object = data.merge.seurat, reduction = 'umap',label = TRUE, group.by = "orig.ident")
    print(p)
    merge.res <- sapply(set.resolutions, function(x){
        p <- DimPlot(object = data.merge.seurat, reduction = 'umap',label = TRUE, group.by = paste0("integrated_snn_res.", x)) + NoLegend()
        print(p)
    })
    p <- DimPlot(object = data.merge.seurat, reduction = 'tsne',label = TRUE, group.by = "orig.ident")
    print(p)
    merge.res <- sapply(set.resolutions, function(x){
        p <- DimPlot(object = data.merge.seurat, reduction = 'tsne',label = TRUE, group.by = paste0("integrated_snn_res.", x)) + NoLegend()
        print(p)
    })
    return(data.merge.seurat)
}