<<<<<<< HEAD
#'@param data: data.frame, store three columns
#ALX3 CYB5A 1
#HESX1 CYB5A 1
#MYF6 CYB5A 1
#VSX1 CYB5A 1
#'@param LinkGroup: character, indicates the variable used for link classification, the default is target

Plot.Sankey <- function(data, saveDir, filename = "Sankey", LinkGroup = "target"){
    require(tidyverse)
    require(viridis)
    require(patchwork)
    require(networkD3)
    require(RColorBrewer)

    #Read the data frame just now and name the column
    colnames(data) <- c("source","target","value")

    #All nodes in the graph
    data.nodes <- data.frame(name=c(as.character(data$source), as.character(data$target)) %>% unique())

    #Renumber the source and target
    data$IDsource <- match(data$source, data.nodes$name)-1
    data$IDtarget <- match(data$target, data.nodes$name)-1

    #This step is to customize the color for flow, and define a column separately
    data$group="flow"

    #Configure the color part, the form is relatively fixed, no need to entangle specific grammatical rules
    #The color corresponding to the code can be viewed with the RColorBrewer package combined with the scales package
    my_color <-'d3.scaleOrdinal()
                .domain(["A", "B", "C", "D", "E", "F", "G", "H", "I", "J","flow"])
                .range(["#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#BC80BD", " #CCEBC5","#BDBDBD"])'

    #Drawing
    Sankey.p <- sankeyNetwork(Links = data, Nodes = data.nodes,
                    Source = "IDsource", Target = "IDtarget",
                    Value = "value", NodeID = "name", LinkGroup = LinkGroup,colourScale = my_color,
                    sinksRight=FALSE, nodeWidth=25, nodePadding=10, fontSize=13,width=900)
    # The main parameters:
    # Links refers to a data frame, including source, target, and value columns. Where source and target use the code replacement or directly correspond to the name.
    # Nodes refers to the names of all points, you can get the names in the links or correspond to the codes in the links yourself.
    # Source, target, value corresponds to the value in the links.
    # NodeID corresponds to the name in Nodes. Here, if the corresponding ID is required, the nodes in links need to be numbered starting from 0.
    # NodeGroup, LinkGroup refers to the changes in the colors of the corresponding nodes and connecting lines. If they are grouped, the colors between different groups will be marked differently.
    # Nodewidth refers to the width of the node.
    #Draw the graphics, you can also directly adjust the upper and lower positions of the square in the R window, but it is still recommended to save it as pdata and use AI to adjust
    #nodeWidth The width of each node block
    #nodePadding The vertical interval between each column of node squares
    #width is the horizontal width of the graphic
    library(htmlwidgets)
    saveWidget(Sankey.p, file = file.path(saveDir, paste0(filename, ".html")))
    library(webshot)
    webshot(file.path(saveDir, paste0(filename, ".html")), file.path(saveDir, paste0(filename, ".pdf")))
    return(NULL)
}