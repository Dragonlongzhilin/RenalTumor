ratio.plot <- function(seurat.object, id.vars1 = "orig.ident", id.vars2 = "seurat_clusters", angle = 45, color.len = c("#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#E5B350")){

    color.len <- color.len[1:length(unique(seurat.object@meta.data[,id.vars1]))]
    a <- sort(table(seurat.object@meta.data[,id.vars2]), decreasing=T)
    cluster.number <- data.frame(Cluster = names(a), number = as.numeric(a))
    cluster.number$Cluster <- factor(cluster.number$Cluster, cluster.number$Cluster)

    patient.Cluster.ratio <- table(seurat.object@meta.data[,c(id.vars1, id.vars2)])
    patient.Cluster.ratio <- as.data.frame(patient.Cluster.ratio)
    cell.num <- tapply(patient.Cluster.ratio$Freq, patient.Cluster.ratio[, id.vars2], sum)
    cell.num <- rep(cell.num, each = length(unique(seurat.object@meta.data[,id.vars1]))) #each指定每个重复多少次
    patient.Cluster.ratio$ratio <- patient.Cluster.ratio$Freq/cell.num
    patient.Cluster.ratio[, id.vars2] <- factor(patient.Cluster.ratio[, id.vars2], levels = cluster.number$Cluster)
    xlabs <- paste0(cluster.number$Cluster, " (n=", cluster.number$number, " cells)")
    library(ggplot2)
    colnames(patient.Cluster.ratio)[2] <- "Type2"
    colnames(patient.Cluster.ratio)[1] <- "Type1"
    p <- ggplot(data = patient.Cluster.ratio, aes(x = Type2, y = ratio, fill = Type1)) + 
    theme_bw()+
    geom_bar(stat ='identity', position ='fill',width = 0.5)+ #Stacked graph, position = fill means stacked graph
     labs(x ='',y ='Ratio',fill = NULL)+ #Define the axis and legend title
     scale_fill_manual(values = color.len) +#Custom colors, you can use `library(RColorBrewer);display.brewer.all()` to display all the options
     scale_y_continuous(labels = c("0%","25%","50%", "75%", "100%")) + ## Percentage axis (need to load scales package)
     scale_x_discrete(labels = xlabs) +
     theme(axis.text.x = element_text(angle = angle, hjust=1), #x-axis label deflection 45°, and drop 0.5
             panel.grid = element_blank(),
             legend.position ='right',
             legend.key.height = unit(0.6,'cm')) # Aspect ratio, expressed in y / x
     print(p)
}