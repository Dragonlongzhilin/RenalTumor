
####Drawing an interactive circular plot, currently focusing on the interaction of the ligand-receptor in each cell
#refer to:
#https://www.r-graph-gallery.com/229-several-circular-plots-in-a-figure.html
#https://jokergoo.github.io/circlize_book/book/graphics.html#links

#'@author Zhilin Long
#'@param track.matrix matrix, specify the length of each track, and its row name is circular drawing order. Two columns are required, starting point coordinates (usually 0), end point coordinates, the difference between the two is the track length
#'[,1] [,2]
#Activated CD4+ T cells-C1 0 1101
#Activated CD4+ T cells-C2 0 1121
#Activated CD4+ T cells-C3 0 1307
#'@param interaction.matrix matrix, specific interaction pair. Two columns are required, the first column is the starting point, and the second column is the ending point; and the name should be consistent with the row name in the track; if lwd.set is T, you need to add a column to set the interaction strength
#'CC1 CC2
# "Cancer cells" "Activated CD4+ T cells-C1" "2.506"
# "Cancer cells" "Activated CD4+ T cells-C2" "2.172"
# "Cancer cells" "Activated CD4+ T cells-C3" "2.25"
# "Cancer cells" "cDC2" "2.532"
#'@param factors character vector, specify the drawing order of the track in circular
#'@param axis.scale Numerical vector, specify the scale of the link movement, easy to distinguish links from different sources. Usually interact around 100
#'@param bg.col character vector, specify the color of the track
#'@param ylim defaults to c(0,1), currently I don’t know how to use it
#'@param cex value, specify the font size of the track label
#'@param link.col Character vector, which specifies the color of the link. rand_color(1, transparency = 0.5) --- https://github.com/davidmerfield/randomColor
#'@param directional value, indicating the directionality of the link. 0 means undirected, 1 means point1-point2, 2 means bidirectional, -1 means point2-point1
#'@param h.ratio value, indicating the height of the link, used when there are more links. The smaller the value, the more concentrated the link is in the center
#'@param arr.type, arr.lty, arr.lwd, arr.col is used to set the attributes of the arrow

Plot.CircularChart <- function(track.matrix, interaction.matrix, factors = NULL, axis.scale = NULL,
                               bg.col = NULL, ylim = c(0, 1), cex = 1, bg.border = NA, transparency = 0.5,
                               link.col = NULL, directional = 1, h.ratio = 0.5, lwd.set = T, lty = par("lty"),
                               arr.type = "triangle", arr.lty = lty, arr.lwd = lwd, arr.col = col, seeds = 100){
    require(circlize)
    require(tidyr)
    set.seed(seeds)

    tracks <- rownames(track.matrix)
    circos.initialize(tracks, xlim = track.matrix)
    interactions.all <- track.matrix[,2] - track.matrix[,1]

    if(is.null(factors)){
        factors <- tracks
    }else{
        factors <- factors
    }

    if(is.null(axis.scale)){
        axis.scale <- mean(track.matrix[,2])/30
    }else{
        axis.scale <- axis.scale
    }

    if(lwd.set){
        interaction.strength <- as.data.frame(interaction.matrix)
        colnames(interaction.strength) <- c("x", "y", "z")
        lwd.info <- spread(as.data.frame(interaction.strength), y, z)
        rownames(lwd.info) <- lwd.info[,1]
        lwd.info <- lwd.info[,-1, drop = F] ##row is source, column is target
    }else{
        lwd.info <- par("lwd")
    }

    if(is.null(bg.col)){
        color.info <- rand_color(nrow(track.matrix), transparency = transparency) 
        circos.trackPlotRegion(ylim = ylim, 
                                track.height = 0.05, 
                                bg.col = color.info, 
                                bg.border = bg.border, 
                                panel.fun=function(x,y){
                                    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = cex)
                                    })                                  
    }else{
        circos.trackPlotRegion(ylim = ylim, 
                        track.height = 0.05, 
                        bg.col = bg.col, 
                        bg.border = bg.border, 
                        panel.fun=function(x,y){
                            circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = cex)
                            }) 
        color.info <- bg.col                          
    }
    names(color.info) <- tracks 

## add links
    #Need to design a list to record the position of the current link, so as to facilitate the corresponding moving position of the subsequent link
    position.list <- list()
    position.lists <- sapply(tracks, function(x){
        return(c(position.list, list(0)))
    })

    for(x in tracks){
        index <- which(interaction.matrix[, 1] == x) #The starting point is the number of interactions of the specified type

        #If there is no corresponding interaction, skip
        if(length(index)==0){
            next
        }else{
            specific.interactions <- interaction.matrix[index,] #Extract interaction entries

            #When there is only one entry, the data structure is character
            if(class(specific.interactions) == "character"){
                specific.interactions <- matrix(specific.interactions, ncol=3)
            }
            #Multiple entries, then matrix
            for(i in 1:nrow(specific.interactions)){
                index1 <- which(factors == specific.interactions[i, 1]) #Get the position of the starting point in the circular plot, corresponding to factors
                index2 <- which(factors == specific.interactions[i, 2]) #Get the position of the termination point in the circular plot, corresponding to factors

                if(lwd.set){
                    # When interacting with oneself, the starting point and ending point are both self, and the position can only reach nrow(specific.interactions)-1
                    if(index1 == index2){
                        circos.link(specific.interactions[i, 1], c(position.lists[[specific.interactions[i, 1]]], position.lists[[specific.interactions[i, 1]]]+axis.scale ), #Starting point position and width
                                    specific.interactions[i, 2], c(interactions.all[index1]-1, interactions.all[index1]), #If there is self-interaction, the end is at the end
                                    col = color.info[specific.interactions[i, 1]], #color and transparency
                                    directional = directional, lwd = as.numeric(lwd.info[specific.interactions[i, 1], specific.interactions[i, 2]]), lty = lty,
                                    h.ratio = h.ratio)

                        position.lists[[specific.interactions[i, 1]]] <- position.lists[[specific.interactions[i, 1]]]+axis.scale #Coordinate movement, specify the position for subsequent links
                    }else{
                    #Other Interactions
                        circos.link(specific.interactions[i, 1], c(position.lists[[specific.interactions[i, 1]]], position.lists[[specific.interactions[i, 1]]]+axis.scale ),
                                    specific.interactions[i, 2], c(position.lists[[specific.interactions[i, 2]]], position.lists[[specific.interactions[i, 2]]]+axis.scale),
                                    col = color.info[specific.interactions[i, 1]],
                                    directional = directional, lwd = as.numeric(lwd.info[specific.interactions[i, 1], specific.interactions[i, 2]]), lty = lty,
                                    h.ratio = h.ratio)
                        position.lists[[specific.interactions[i, 1]]] <- position.lists[[specific.interactions[i, 1]]]+axis.scale
                        position.lists[[specific.interactions[i, 2]]] <- position.lists[[specific.interactions[i, 2]]]+axis.scale #Coordinate movement, specify the position for subsequent links
                    }
                }else{
                    #When interacting with oneself, the starting point and ending point are both self, and the position can only reach nrow(specific.interactions)-1
                    if(index1 == index2){
                        circos.link(specific.interactions[i, 1], c(position.lists[[specific.interactions[i, 1]]], position.lists[[specific.interactions[i, 1]]]+axis.scale ), #Starting point position and width
                                    specific.interactions[i, 2], c(interactions.all[index1]-1, interactions.all[index1]), #If there is self-interaction, the end is at the end
                                    col = color.info[specific.interactions[i, 1]], #color and transparency
                                    directional = directional, lwd = as.numeric(lwd.info), lty = lty,
                                    h.ratio = h.ratio)

                        position.lists[[specific.interactions[i, 1]]] <- position.lists[[specific.interactions[i, 1]]]+axis.scale #Coordinate movement, specify the position for subsequent links
                    }else{
                    #Other Interactions
                        circos.link(specific.interactions[i, 1], c(position.lists[[specific.interactions[i, 1]]], position.lists[[specific.interactions[i, 1]]]+axis.scale ),
                                    specific.interactions[i, 2], c(position.lists[[specific.interactions[i, 2]]], position.lists[[specific.interactions[i, 2]]]+axis.scale),
                                    col = color.info[specific.interactions[i, 1]],
                                    directional = directional, lwd = as.numeric(lwd.info), lty = lty,
                                    h.ratio = h.ratio)
                        position.lists[[specific.interactions[i, 1]]] <- position.lists[[specific.interactions[i, 1]]]+axis.scale
                        position.lists[[specific.interactions[i, 2]]] <- position.lists[[specific.interactions[i, 2]]]+axis.scale #Coordinate movement, specify the position for subsequent links
                    }
                }

            }
        }
    }
    circos.clear()
}

