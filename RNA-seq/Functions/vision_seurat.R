# https://yoseflab.github.io/VISION/articles/web_only/Seurat.html
#'@decription: VISION working with seurat
#'reference: Functional interpretation of single cell similarity maps. 2019 Nature commmun

#For the sigData vector, the names represent gene names and the values (1 or -1) represent the ‘sign’ of the gene.
#For an unsigned signature, just use 1 for the value of every gene.

#Based on MsigDB: Hallmark, KEGG, and GO Biological Process gene sets

#seurat.object, try to request that it has been processed by logNormailzed, PCA and UMAP; The original (normalized) counts will be used as the expression input
#customize.signature
#gene type
vision_seurat <- function(seurat.object, customize.signature = NULL, positive = 1, min_signature_genes = 5, signature.dir = "/data/ExtraDisk/sdd/longzhilin/Data/pathway/MsigDB", signature.names = c("h.all.v7.4.symbols.gmt", "c5.go.bp.v7.4.symbols.gmt", "c2.cp.reactome.v7.4.symbols.gmt", "c2.cp.kegg.v7.4.symbols.gmt", "c2.cp.biocarta.v7.4.symbols.gmt"), mc.cores = 32){

    require(VISION)
    options(mc.cores = mc.cores)
    DefaultAssay(seurat.object) <- "RNA"

    #construct signature set
    if(is.null(customize.signature)){
        signature.genesets <- paste0(signature.dir, "/", signature.names)
    }else{
        signature.genesets <- c()
        for(i in names(table(customize.signature$Type))){
            genes <- customize.signature$Gene[which(customize.signature$Type==i)]
            sigData <- rep(positive, length(genes))
            names(sigData) <- genes
            sig <- createGeneSignature(name = i, sigData = sigData)
            signature.genesets <- c(signature.genesets, sig)
        }
    }

    vision.obj <- Vision(seurat.object, signatures = signature.genesets, min_signature_genes = min_signature_genes)
    vision.obj <- analyze(vision.obj)
    return(vision.obj)
}