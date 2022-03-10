#### convert gene ID
#author: Zhilin long
#time: 2021.1.19

#'@author Zhilin Long
#'@param gene Character vector, ID vector to be converted
#'@param fromType Character vector, initial ID type
#'@param toType Character vector, ID type to be converted
#'
#'@return data.frame has removed duplication and NA, and can be used directly


IDConvert <- function(genes, 
                       method = c("clusterProfiler", "org.Hs.eg.db"),
                       fromType = c("ENTREZID", "SYMBOL", "ENSEMBL"), 
                       toType = c("ENTREZID", "SYMBOL", "ENSEMBL")
                       ){
    require(org.Hs.eg.db)
    require(clusterProfiler)

    method <- match.arg(method)
    fromType <- match.arg(fromType)
    toType <- match.arg(toType)

    if(method == "clusterProfiler"){
        gene.df <- bitr(genes, fromType = fromType, toType = toType, OrgDb = "org.Hs.eg.db")
        #By default, IDs that are not mapped will be filtered, that is, the entries will be reduced
    }else{
        gene.df <- mapIds(org.Hs.eg.db, 
                          keys = genes, 
                          column = toType, 
                          keytype = fromType,
                          multiVals="first")
        gene.df <- data.frame(ENSEMBL = names(gene.df), SYMBOL = gene.df)
        colnames(gene.df) <- c(fromType, toType)
        # If there is no map, it is NA, that is, the entry remains unchanged
        gene.df <- na.omit(gene.df) 
    }

    #remove duplicates
    gene.df <- gene.df[!duplicated(gene.df[,2]),]
    gene.df <- gene.df[!duplicated(gene.df[,1]),]
    return(gene.df)
}