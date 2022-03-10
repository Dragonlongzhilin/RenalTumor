# Filter genes that express less than 3 cells

Filter.gene <- function(Seurat.object, slot = "counts", min.cell = 3, assay = "RNA"){
    # Remove genes less than 3 cells expressed
    expression_matrix <- GetAssayData(Seurat.object, slot = slot, assay = assay)
    # The gene is required to be expressed in at least 3 cells
    exp.num <- apply(expression_matrix, 1, function(x){
        return(length(which(x>0)))
    })
    expression_matrix  <- expression_matrix[which(exp.num>min.cell),]
    Seurat.object.new <- subset(x = Seurat.object, features = rownames(expression_matrix))
    return(Seurat.object.new)
}