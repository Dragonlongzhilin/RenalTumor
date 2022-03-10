vision.plot <- function(vision.obj, groupName, signature.names = NULL, title = "", subgroup = T, fontsize = 10, dot.size = 1, wid = 2, hei = 4, jitter.size = 0.3){
    require(ggsci)
    if(is.null(signature.names)){
        sigScore <- as.data.frame(vision.obj@SigScores)
    }else{
        sigScore <- as.data.frame(vision.obj@SigScores[, signature.names])
    }
    sigScore <- scale(sigScore)

    sigScore$group <- vision.obj@metaData[, groupName]
    require(reshape2)
    vision.group.score <- melt(sigScore, id.vars = "group", variable.name = "Type", value.name = "Score")
    
    #box plot
    p1 <- ggboxplot(vision.group.score, x = "group", y = "Score", color = "Type", xlab = "", palette = "npg", title = title, ylab = "Score", add="jitter", add.params = list(size = jitter.size)) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + stat_compare_means(aes(group = Type))
    print(p1)
    umap <- as.data.frame(vision.obj@Projections$Seurat_umap)
    umap.data <- cbind(umap, sigScore)
    if(subgroup){
        subgroup.label <- colnames(sigScore)[-ncol(sigScore)]
        for(i in subgroup.label){
            index <- which(vision.group.score$Type == i)
            p2 <- ggboxplot(vision.group.score[index,], x = "group", y = "Score", palette = "npg", color = "group", xlab = "", title = i, ylab = "Score", add="jitter", add.params = list(size = jitter.size)) + theme(legend.position="none", axis.text.x = element_text(angle = 30, hjust = 1)) + stat_compare_means(label.y = max(vision.group.score[index,"Score"]))
            print(p2)

            p3 <- ggscatter(umap.data, x = "UMAP_1", y = "UMAP_2", color = i, fill = i, alpha = 0.6, size = dot.size) + gradient_color(c("lightgrey", "red"))
            print(p3)
        }
    }

    #heatmap plot
    vision.meanScore <- apply(sigScore[,-ncol(sigScore)], 2, function(x) {
        a <- tapply(x, vision.obj@metaData[, groupName], mean)
        return(a)
    })
    require(ComplexHeatmap)
    require(circlize)
    p4 <- Heatmap(vision.meanScore, width = unit(wid, "cm"), height = unit(hei, "cm"), show_column_dend = F, show_row_dend = F, row_names_gp = gpar(fontsize = fontsize), column_names_gp = gpar(fontsize = fontsize))
    print(p4)
}
