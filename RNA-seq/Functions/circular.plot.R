circular.plot <- function(gene.a, gene.b, cellphone.interaction.sig, data.merge, factors, size = 2){
    source(file = "/home/longzhilin/Analysis_Code/Plot.CircularChart.R")
    pair1 <- paste(gene.a, gene.b, sep = "_")
    pair2 <- paste(gene.b, gene.a, sep = "_")    
    index <- which(cellphone.interaction.sig$interacting_pair %in% c(pair1, pair2))
    res <- cellphone.interaction.sig[index,]

    # factors <- c("Mononuclear Phagocyte", "Mononuclear Phagocyte like", "Dendritic cell", "Mast cell",
    #          "Proliferating T cell", "CD8+ T cell", "NK cell", "Th cell", "Regulatory T cell", "B cell", "Plasma cell",
    #          "Cancer cell", "Endothelium", "Mesangial cell")
    interactions.cellType <- apply(res, 1, function(x){
        a <- unlist(strsplit(x[12], "\\|"))
        a <- c(a, as.numeric(x[15])*size)
        return(a)
    })
    interactions.cellType <- t(interactions.cellType)
    track.matrix <- matrix(c(rep(0, nrow(total.interaction)), total.interaction$Number), ncol = 2)
    rownames(track.matrix) <- rownames(total.interaction)
    interaction.matrix <- interactions.cellType
    Plot.CircularChart(track.matrix = track.matrix, interaction.matrix = interaction.matrix, bg.col = NULL, cex = 0.8,
                    link.col = NULL, arr.lwd = 0.3)
    p1 <- FeaturePlot(data.merge, features = gene.a, cols = c("lightgrey", "red"), label = F)
    p2 <- FeaturePlot(data.merge, features = gene.b, cols = c("lightgrey", "red"), label = F)
    print(p1)
    print(p2)
    return(res)
}