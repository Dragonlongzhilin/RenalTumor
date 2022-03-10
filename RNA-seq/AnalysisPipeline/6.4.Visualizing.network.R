#' @description: Visualized network view
#Combine R package---RCy3 (V2.10.2) and Cytoscape (V 3.8.2)
#Since cytoscape cannot be started on the server, it is currently being visualized locally

############################################## ########All interactions############################
####Step 1: Get network data
#91Server
#Since network.txt cannot distinguish the directional interaction, it is based on your own processing
library(dplyr)
library(tidyverse)
library(ggpubr)
setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/6.CrossTalk/CellPhoneDB")
count_network <- readRDS("All.interaction.number.rds") # mean>=1 & pvalue < 0.05
node.info <- count_network[,1:2]

mypvals <- readRDS("mypvals.usr.rds")
mymeans <- readRDS("mymeans.usr.rds")
mymeans <- mymeans[,-c(1,3:10,12)]
mypvals <- mypvals[,-c(1,3:10,12)]
mymeans %>% reshape2::melt(id.vars = c("interacting_pair", "ligand")) -> meansdf
colnames(meansdf)<- c("interacting_pair","ligand","CC","means")
mypvals %>% reshape2::melt(id.vars = c("interacting_pair", "ligand"))-> pvalsdf
colnames(pvalsdf)<- c("interacting_pair","ligand","CC","pvals")
pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair, "_", pvalsdf$CC)
meansdf$joinlab<- paste0(meansdf$interacting_pair, "_", meansdf$CC)
pldf <- merge(pvalsdf,meansdf,by = "joinlab")
pldf <- pldf[,-c(6:8)]
pldf.pro <- pldf%>% filter(means>=1 & pvals<0.05) # 13421head
pldf.pro[,4] <- as.character(pldf.pro[,4])

#Adjust the way of interaction --- ligand-receptor
pldf.pro <- apply(pldf.pro, 1, function(x){
    if(x[3] == "Receptor"){

        pairs <- unlist(strsplit(x = x[2], split = "_"))
        x[2] <- paste0(pairs[2], "_", pairs[1])

        interactions <- unlist(strsplit(x = x[4], split = "\\|"))
        x[4] <- paste0(interactions[2], "|", interactions[1])

        x[1] <- paste0(x[2], "|", x[4])

        x[3] <- "Ligand"
    }
    return(x)
})
pldf.pro <- as.data.frame(t(pldf.pro))

edge.num <- as.data.frame(table(as.character(pldf.pro$CC.x)))
edge.info <- apply(edge.num, 1, function(x){
    inter <- unlist(strsplit(x[1], "\\|"))
    return(c(inter, x[2]))
})
edge.info <- as.data.frame(t(edge.info)) #189
colnames(edge.info) <- c("Source", "Target", "Number")
#Replace CD8+ T-IEG back to CD8+ T-Exhausted IEG
edge.info$Source <- gsub("CD8\\+ T-IEG", "CD8\\+ T-Exhausted IEG", edge.info$Source)
edge.info$Target <- gsub("CD8\\+ T-IEG", "CD8\\+ T-Exhausted IEG", edge.info$Target)
node.info$cellType <- gsub("CD8\\+ T-IEG", "CD8\\+ T-Exhausted IEG", node.info$cellType)

##Building network information
class.type <- rep("Lymphoid", nrow(node.info))
class.type[c(1, 8, 13:16, 19)] <- "Myeloid"
class.type[10] <- "Tumor"
class.type[c(4, 6, 12)] <- "Other"
nodes <- data.frame(id = node.info[,1],
                    group = class.type, # categorical strings
                    score = as.integer(node.info[,2]), # integers
                    stringsAsFactors=FALSE)
edges <- data.frame(source = edge.info$Source,
                    target = edge.info$Target,
                    weight = as.numeric(edge.info$Number), # numeric
                    stringsAsFactors=FALSE)
write.table(nodes, file = "nodes.txt", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(edges, file = "edges.txt", sep = "\t", row.names = F, col.names = T, quote = F)


library(RCy3)
cytoscapePing()
cytoscapeVersionInfo()
setwd("D:/work/Renal Tumor/Result/DataAnalysis-2021.8.3/scRNA/6.CrossTalk/CellPhoneDB")
nodes <- read.table("nodes.txt", header = T, stringsAsFactors = F, sep = "\t")
edges <- read.table("edges.txt", header = T, stringsAsFactors = F, sep = "\t")

tumor.colors <- rep(0, nrow(edges))
#tumor- target
idx <- which(edges$target=="Tumor")
tumor.colors[idx] <- 2
#tumor-source
idx <- which(edges$source=="Tumor")
tumor.colors[idx] <- 1
edges$tumor <- tumor.colors

#Standardize the points to a size of 1-10
library(vegan)
nodes$weight <- round(nodes[,3]/100)
createNetworkFromDataFrames(nodes = nodes,edges = edges, title="cross-talk", collection="M1")

style.name = "myStyle"
defaults <- list(NODE_SHAPE="ELLIPSE",
                 EDGE_TRANSPARENCY=120)
nodeLabels <- mapVisualProperty('node label','id','p')
nodesizes <- mapVisualProperty("node size","weight","p")
#node.border.colors <- mapVisualProperty("node border paint","group","d",c("Lymphoid","Myeloid", "CancerCell"), c("#34A047","#00B3F1", "#EF7F48"))
#nodeFills <- mapVisualProperty("node fill color","group","d",c("Lymphoid","Myeloid", "CancerCell"), c("#34A047","#00B3F1", "#EF7F48"))
edgeWidth <- mapVisualProperty("edge width","weight","p")
createVisualStyle(style.name, defaults, list(nodeLabels, nodesizes, edgeWidth))

#pal_npg("nrc")(4)
# [1] "#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF"
library(ggplot2)
library(ggsci)

node.colors <- nodes$group
node.colors <- gsub("Lymphoid", "#00A087", node.colors)
node.colors <- gsub("Myeloid", "#4DBBD5", node.colors)
node.colors <- gsub("Tumor", "#E64B35", node.colors)
node.colors <- gsub("Other", "#3C5488", node.colors)
setNodeColorBypass(node.names = nodes$id, new.colors = node.colors)
setNodeBorderColorMapping(table.column = 'group', table.column.values = c("Lymphoid","Myeloid", "Tumor", "Other"), colors = c("#00A087","#4DBBD5", "#E64B35", "#3C5488"), mapping.type = 'd', style.name = style.name)
setEdgeTargetArrowShapeDefault(new.shape = "ARROW", style.name = style.name)
setEdgeColorMapping(table.column = c('tumor'), table.column.values = c(0,1,2), colors = c("#CCCCCC", "#F15F30", "#1e90ff"), mapping.type = 'd', style.name = style.name)
setVisualStyle(style.name)


##########################################################Tumor interactions#############################
library(dplyr)
library(tidyverse)
library(ggpubr)
setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/6.CrossTalk/CellPhoneDB")
count_network <- readRDS("Tumor/tumor.interaction.number.rds") # mean>=1 & pvalue < 0.05
node.info <- count_network[,1:2]

mypvals <- readRDS("mypvals.usr.rds")
mymeans <- readRDS("mymeans.usr.rds")
mymeans <- mymeans[,-c(1,3:10,12)]
mypvals <- mypvals[,-c(1,3:10,12)]
mymeans %>% dplyr::select("interacting_pair", "ligand", starts_with("Tumor"), ends_with("Tumor")) %>% reshape2::melt() -> meansdf
colnames(meansdf)<- c("interacting_pair","ligand","CC","means")
mypvals %>% dplyr::select("interacting_pair", "ligand", starts_with("Tumor"), ends_with("Tumor")) %>% reshape2::melt()-> pvalsdf
colnames(pvalsdf)<- c("interacting_pair", "ligand", "CC","pvals")
pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair, "_", pvalsdf$CC)
meansdf$joinlab<- paste0(meansdf$interacting_pair, "_", meansdf$CC)
pldf <- merge(pvalsdf,meansdf,by = "joinlab")
pldf <- pldf[,-c(6:8)]
pldf.pro <- pldf%>% filter(means>=1 & pvals<0.05)
pldf.pro[,4] <- as.character(pldf.pro[,4])

pldf.pro <- apply(pldf.pro, 1, function(x){
    if(x[3] == "Receptor"){
        pairs <- unlist(strsplit(x = x[2], split = "_"))
        x[2] <- paste0(pairs[2], "_", pairs[1])

        interactions <- unlist(strsplit(x = x[4], split = "\\|"))
        x[4] <- paste0(interactions[2], "|", interactions[1])

        #调整joinlab
        x[1] <- paste0(x[2], "|", x[4])

        x[3] <- "Ligand"
    }
    return(x)
})
pldf.pro <- as.data.frame(t(pldf.pro))
edge.num <- as.data.frame(table(as.character(pldf.pro$CC.x)))
edge.info <- apply(edge.num, 1, function(x){
    inter <- unlist(strsplit(x[1], "\\|"))
    return(c(inter, x[2]))
})
edge.info <- as.data.frame(t(edge.info)) #189
colnames(edge.info) <- c("Source", "Target", "Number")
edge.info$Source <- gsub("CD8\\+ T-IEG", "CD8\\+ T-Exhausted IEG", edge.info$Source)
edge.info$Target <- gsub("CD8\\+ T-IEG", "CD8\\+ T-Exhausted IEG", edge.info$Target)
node.info$cellType <- gsub("CD8\\+ T-IEG", "CD8\\+ T-Exhausted IEG", node.info$cellType)

class.type <- rep("Lymphoid", nrow(node.info))
class.type[c(7,10,12:13,16:18)] <- "Myeloid"
class.type[20] <- "Tumor"
class.type[c(8:9, 11)] <- "Other"
nodes <- data.frame(id = node.info[,1],
                    group = class.type, # categorical strings
                    score = as.integer(node.info[,2]), # integers
                    stringsAsFactors=FALSE)
edges <- data.frame(source = edge.info$Source,
                    target = edge.info$Target,
                    weight = as.numeric(edge.info$Number), # numeric
                    stringsAsFactors=FALSE)
write.table(nodes, file = "Tumor/nodes.txt", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(edges, file = "Tumor/edges.txt", sep = "\t", row.names = F, col.names = T, quote = F)

library(RCy3)
cytoscapePing()
cytoscapeVersionInfo()
setwd("D:/work/Renal Tumor/Result/DataAnalysis-2021.8.3/scRNA/6.CrossTalk/CellPhoneDB/Tumor")
nodes <- read.table("nodes.txt", header = T, stringsAsFactors = F, sep = "\t")
edges <- read.table("edges.txt", header = T, stringsAsFactors = F, sep = "\t")

tumor.colors <- rep(0, nrow(edges))
#tumor- target
idx <- which(edges$target=="Tumor")
tumor.colors[idx] <- 2
#tumor-source
idx <- which(edges$source=="Tumor")
tumor.colors[idx] <- 1
edges$tumor <- tumor.colors

library(vegan)
nodes$weight <- round(nodes[,3]/10)
idx <- which(nodes$id == "Tumor")
nodes$weight[idx] <- 15
createNetworkFromDataFrames(nodes = nodes,edges = edges, title="cross-talk", collection="M1")

style.name = "myStyle"
defaults <- list(NODE_SHAPE="ELLIPSE",
                 EDGE_TRANSPARENCY=120)
nodeLabels <- mapVisualProperty('node label','id','p')
nodesizes <- mapVisualProperty("node size","weight","p")
#node.border.colors <- mapVisualProperty("node border paint","group","d",c("Lymphoid","Myeloid", "CancerCell"), c("#34A047","#00B3F1", "#EF7F48"))
#nodeFills <- mapVisualProperty("node fill color","group","d",c("Lymphoid","Myeloid", "CancerCell"), c("#34A047","#00B3F1", "#EF7F48"))
edgeWidth <- mapVisualProperty("edge width","weight","p")
createVisualStyle(style.name, defaults, list(nodeLabels, nodesizes, edgeWidth))

library(ggplot2)
library(ggsci)

node.colors <- nodes$group
node.colors <- gsub("Lymphoid", "#00A087", node.colors)
node.colors <- gsub("Myeloid", "#4DBBD5", node.colors)
node.colors <- gsub("Tumor", "#E64B35", node.colors)
node.colors <- gsub("Other", "#3C5488", node.colors)
setNodeColorBypass(node.names = nodes$id, new.colors = node.colors)
setNodeBorderColorMapping(table.column = 'group', table.column.values = c("Lymphoid","Myeloid", "Tumor", "Other"), colors = c("#00A087","#4DBBD5", "#E64B35", "#3C5488"), mapping.type = 'd', style.name = style.name)
setEdgeTargetArrowShapeDefault(new.shape = "ARROW", style.name = style.name)
setEdgeColorMapping(table.column = c('tumor'), table.column.values = c(0,1,2), colors = c("#CCCCCC", "#F15F30", "#1e90ff"), mapping.type = 'd', style.name = style.name)
setVisualStyle(style.name)