#' @description get variable genes and evalute the cell cycle effect
#' @param seurat.lists: multi seurat object list
#' @parma method: seurat process method
#' @param: return.only.var.gene, vars.to.regress: parameter for SCTranform
#' @param: nfeatures: the number of variable features

variableFeatureSelection <- function(seurat.lists, method = c("Standard", "SCT"), return.only.var.genes = T, nfeatures = 3000, vars.to.regress = c("nCount_RNA", "percent.mt")){

    for (i in 1:length(seurat.lists)) {
        DefaultAssay(seurat.lists[[i]]) <- "RNA"
        if(method == "SCT"){
            seurat.lists[[i]] <- SCTransform(seurat.lists[[i]], variable.features.n = nfeatures, vars.to.regress = vars.to.regress, verbose = FALSE, return.only.var.genes = return.only.var.genes)
        }else{
            seurat.lists[[i]] <- NormalizeData(seurat.lists[[i]], verbose = FALSE)
            seurat.lists[[i]] <- FindVariableFeatures(seurat.lists[[i]], selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)            
        }
    }
    return(seurat.lists)
}