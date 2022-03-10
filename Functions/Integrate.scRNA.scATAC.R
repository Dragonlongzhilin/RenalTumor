#ref: Single cell transcriptional and chromatin accessibility profiling redefine cellular heterogeneity in the adult human kidney

#'@param ref.npcs = 30, dims = 30, are the parameters in the FindTransferAnchors function, used to control the number of dimensions used by reference data and anchor search space respectively
#'@param scATAC.assay = "ACTIVITY", scATAC.normalize = T, PC = 30, is used to process scATAC data.
# scATAC.assay indicates the name of the assay where the gene activity is stored in the scATAC number
# scATAC.normalize Whether the gene activity has been logNormalized
# PC indicates the dimension of scATAC dimensionality reduction, which comes from the lsi of scATAC.object itself. Note that the default dimensionality reduction assay of scATAC.object is lsi

Integrate.scRNA.scATAC <- function(scATAC.object, scRNA.object, class.label = "cellType", scRNA.assay = "RNA", ref.npcs = 30, dims = 30, scATAC.assay = "ACTIVITY", scATAC.normalize = T, PC = 30, nfeatures = 3000, SCT = F){

    ##1. load data
    #scRNA-seq data
    DefaultAssay(scRNA.object) <- scRNA.assay
    scRNA.object <- NormalizeData(scRNA.object)
    if(SCT){
        #Need to use sctransform to normalize data in advance
        var.genes <- scRNA.object@assays$SCT@var.features[1:nfeatures]
    }else{
        var.genes <- VariableFeatures(scRNA.object)[1:nfeatures]
    }
    #scATAC-seq data
    DefaultAssay(scATAC.object) <- scATAC.assay
    if(!scATAC.normalize){
        scATAC.object <- NormalizeData(scATAC.object)
    }
    scATAC.object <- ScaleData(scATAC.object, features = rownames(scATAC.object))

    ##2. Identify anchors
    #npcs = 30; dim = 1:30
    transfer.anchors <- FindTransferAnchors(reference = scRNA.object, 
                                            query = scATAC.object, 
                                            features = var.genes, #SCT model：scRNA.merge@assays$SCT@var.features
                                            reference.assay = scRNA.assay,
                                            query.assay = scATAC.assay, 
                                            reduction = "cca",
                                            npcs = ref.npcs, # Number of PCs to compute on reference if reference
                                            dims = 1:dims) # Which dimensions to use from the reduction to specify the neighbor search space
    #saveRDS(transfer.anchors, file = file.path(saveDir, "transfer.anchors.rds"))

    #The default dims = NULL, it is consistent with FindTransferAnchors
    celltype.predictions <- TransferData(anchorset = transfer.anchors, 
                                refdata = scRNA.object@meta.data[,class.label], 
                                weight.reduction = scATAC.object[["lsi"]], 
                                dims = 2:PC)    


    scATAC.object <- AddMetaData(scATAC.object, metadata = celltype.predictions)
    scATAC.object <- AddMetaData(scATAC.object, metadata = celltype.predictions$predicted.id, class.label)

    p <- gghistogram(scATAC.object@meta.data, x = "prediction.score.max", y = "..density..", main = "Prediction Score for scATAC", add_density = TRUE)
    print(p)    
    p <- DimPlot(scATAC.object, group.by = class.label, label = TRUE, repel = TRUE, reduction = "umap") +
        ggtitle("scATAC-seq Predicted Celltypes Before Harmony") + NoLegend()
    print(p)
    p <- DimPlot(scATAC.object, group.by = class.label, label = TRUE, repel = TRUE, reduction = "umap") +
        ggtitle("scATAC-seq Predicted Celltypes Before Harmony")
    print(p)    
    p <- DimPlot(scATAC.object, group.by = class.label, label = TRUE, repel = TRUE, reduction = "tsne") +
        ggtitle("scATAC-seq Predicted Celltypes Before Harmony") + NoLegend()
    print(p)
    p <- DimPlot(scATAC.object, group.by = class.label, label = TRUE, repel = TRUE, reduction = "umap") +
        ggtitle("scATAC-seq Predicted Celltypes Before Harmony")
    print(p)
    return(list(transfer.anchors = transfer.anchors, scATAC.object = scATAC.object))
}

#' @param dims Indicates the dimensionality reduction used when coembedding
Coembedding.scATAC.scRNA <- function(transfer.anchors, scATAC.object, scRNA.assay = "RNA", scRNA.object, class.label = "cellType", PC = 30, dims = 30){
    
    ####Co-embedding scRNA-seq and scATAC-seq datasets
    # note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
    # full transcriptome if we wanted to
    genes.use <- transfer.anchors@anchor.features
    DefaultAssay(scRNA.object) <- scRNA.assay
    scRNA.object <- NormalizeData(scRNA.object)
    refdata <- GetAssayData(scRNA.object, assay = scRNA.assay, slot = "data")[genes.use, ]
    # refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
    # (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
    imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = scATAC.object[["lsi"]], dims = 2:PC)
    scATAC.object[["RNA"]] <- imputation
    DefaultAssay(scATAC.object) <- "RNA"     

    scRNA.object <- AddMetaData(scRNA.object, metadata = rep("scRNA", nrow(scRNA.object@meta.data)), col.name = "type")
    scATAC.object <- AddMetaData(scATAC.object, metadata = rep("scATAC", nrow(scATAC.object@meta.data)), col.name = "type")

    coembed <- merge(x = scRNA.object, y = scATAC.object)
    DefaultAssay(coembed) <- "RNA"
    # Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
    # datasets
    coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
    coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
    coembed <- RunHarmony(object = coembed, group.by.vars = "orig.ident", verbose = FALSE)
    coembed <- RunUMAP(coembed, dims = 1:dims, reduction = 'harmony', verbose = FALSE)
    coembed <- RunTSNE(coembed, dims = 1:dims, reduction = 'harmony', check_duplicates = FALSE, verbose = FALSE)
    #saveRDS(coembed, file = file.path(saveDir, "coembed.rds"))

    p <- ElbowPlot(object = coembed, ndims = dims)
    print(p)
    p <- DimPlot(coembed, group.by = c("orig.ident"), reduction = "umap", label = T) + NoLegend()
    print(p)
    p <- DimPlot(coembed, group.by = c("orig.ident"), reduction = "tsne", label = T) + NoLegend()
    print(p)
    p <- DimPlot(coembed, group.by = c("type"), reduction = "umap", label = T) + NoLegend()
    print(p)
    p <- DimPlot(coembed, group.by = class.label, reduction = "umap", label = T) + NoLegend()
    print(p)
    p <- DimPlot(coembed, group.by = c("type"), reduction = "tsne", label = T) + NoLegend()
    print(p)
    p <- DimPlot(coembed, group.by = class.label, reduction = "tsne", label = T) + NoLegend()
    print(p)

    p <- DimPlot(coembed, group.by = c("orig.ident"), reduction = "umap", label = T)
    print(p)
    p <- DimPlot(coembed, group.by = c("orig.ident"), reduction = "tsne", label = T)
    print(p)
    p <- DimPlot(coembed, group.by = c("type"), reduction = "umap", label = T)
    print(p)
    p <- DimPlot(coembed, group.by = class.label, reduction = "umap", label = T)
    print(p)
    p <- DimPlot(coembed, group.by = c("type"), reduction = "tsne", label = T)
    print(p)
    p <- DimPlot(coembed, group.by = class.label, reduction = "tsne", label = T)
    print(p)
    return(list(coembed = coembed, scATAC.object = scATAC.object))
}
