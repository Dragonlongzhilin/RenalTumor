#Draw the expression of features in each cluster

#'@param data Generally requires the expression matrix after log, behavioral genes, listed as samples
#'@param features character vector,
#'@param clusters The sample classification situation is required to be a vector, that is, names (clusters) correspond to the column names of the expression matrix, and you need to arrange the display order yourself

plot.violin <- function(data, features, clusters, clusters.levels = NULL, xlab = NULL, ylab = NULL, title = NULL, normalizedMethod = "min-max"){
    
    # Corresponding gene expression
    index <- match(features, rownames(data))
    data <- data[na.omit(index),]
    if(normalizedMethod == "min-max"){
        data <- decostand(data, "range", 1)
    }
    # Corresponding cell classification
    index <- match(names(clusters), colnames(data))
    data <- data[, index]
    data <- as.data.frame(data)
    data$gene <- rownames(data)

    require(tidyr)
    data.reshape <- gather(data, sample, expression, colnames(data)[-length(colnames(data))])
    data.reshape$cluster <- rep(clusters, each = nrow(data))
    if(!is.null(clusters.levels)){
        data.reshape$cluster <- factor(data.reshape$cluster, levels = clusters.levels)
    }

    require(ggpubr)
    p <- ggviolin(data.reshape, x = "cluster", y = "expression", fill = "cluster", color = "cluster", add = "mean", add.params = list(color = "black", size = 0.01), xlab = xlab, ylab = ylab, title = title)
    p <- p + facet_grid(gene ~ ., scales="free_y") + theme(legend.position="none")
    return(p)
}