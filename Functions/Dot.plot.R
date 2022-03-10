#' @param exp.matrix Expression matrix, row is genes, column is samples
#' @param group.info Vector, grouping information
Dot.plot <- function(exp.matrix, group.info, features, group.by, id.levels = NULL, x = "feature", y = "groups",
                    cols = c("lightgrey","blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
                    scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA){
    require(ggplot2)
    require(ggpubr)
    require(tidyr)
    require(cowplot)
    require(Seurat)

    scale.func <- switch(EXPR = scale.by, size = scale_size,
        radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
    # ATAC.matrix <- getMatrixFromProject(projATAC, useMatrix = useMatrix)
    # #expression matrix
    # exp.matrix <- assay(ATAC.matrix)
    # #rowname
    # row.info <- rowData(ATAC.matrix)
    # rownames(exp.matrix) <- row.info$name

    # Extract the expression of a specific feature
    index <- match(features, rownames(exp.matrix))
    index.na <- which(is.na(index))
    if(length(index.na)>0){
        cat(features[index.na], "not exist!\n")
    }
    features.plot <- features[which(!is.na(index))]
    feature.matrix <- exp.matrix[na.omit(index),]

    # Extract group information
    # index <- match(colnames(feature.matrix), rownames(projATAC@cellColData))
    # meta.info <- projATAC@cellColData[index,]
    # group.info <- as.character(meta.info[, group.by])

    # Calculate the average expression and expression ratio of each cluter marker
    # ref: https://rdrr.io/github/satijalab/seurat/src/R/visualization.R --- seurat
    data.plot <- apply(feature.matrix, 1, function(x){
        avg.exp <- tapply(x, group.info, function(y){
            return(mean(expm1(y)))
        })
        pct.exp <- tapply(x, group.info, function(y){
            r <- length(which(y>0))/length(y)
            return(r)
        })
        return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    })
    avg.exp <- sapply(data.plot, function(x){
        return(x$avg.exp)
    })

    pct.exp <- sapply(data.plot, function(x){
        return(x$pct.exp)
    })
    data.plot <- list()
    data.plot$avg.exp <- avg.exp
    data.plot$pct.exp <- pct.exp   
    data.plot$features.plot <- features.plot

    # Standardized average value as a color scale
    data.use <- data.plot$avg.exp
    if(scale){
        data.use <- scale(x = data.use)
        data.use <- MinMax(data = data.use, min = col.min, max = col.max)
    }else{
        data.use <- log(x = data.use)
    }
    data.plot$avg.exp.scaled <- data.use

    data.plot$features.plot <- factor(x = data.plot$features.plot, levels = features.plot)
    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100

    if (!is.na(x = scale.min)) {
        data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }

    avg.plot <- data.frame(groups = rownames(data.plot$avg.exp), data.plot$avg.exp)
    avg.plot <- gather(avg.plot, "feature", "avg.exp", -1)
    pct.plot <- data.frame(groups = rownames(data.plot$pct.exp), data.plot$pct.exp)
    pct.plot <- gather(pct.plot, "feature", "pct.exp", -1)
    avg.scale.plot <- data.frame(groups = rownames(data.plot$avg.exp.scaled), data.plot$avg.exp.scaled)
    avg.scale.plot <- gather(avg.scale.plot, "feature", "avg.scale.exp", -1)   

    data.plot.pro <- merge(avg.plot, pct.plot, by = c("groups", "feature")) #合并数据
    data.plot.pro <- merge(data.plot.pro, avg.scale.plot, by = c("groups", "feature"))
    data.plot.pro$feature <- factor(data.plot.pro$feature, levels = features.plot)
    
    if(is.null(id.levels)){
        data.plot.pro$groups <- factor(data.plot.pro$groups)
    }else{
        data.plot.pro$groups <- factor(data.plot.pro$groups, levels = id.levels)
    }

    plot <- ggplot(data = data.plot.pro, mapping = aes_string(x = x, y = y)) + geom_point(mapping = aes_string(size = "pct.exp", color = "avg.scale.exp")) 
    plot <- plot + scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + labs(x = "", y = "")+ theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) + theme_cowplot()

    if(length(cols) == 2){
        plot <- plot + scale_color_gradient(low = cols[1], high = cols[2]) + guides(color = guide_colorbar(title = "Average Expression"))
    }else{
        plot <- plot + scale_color_gradient2(low = cols[1], mid = cols[2], high = cols[3]) + guides(color = guide_colorbar(title = "Average Expression"))
    }
    return(plot)
}