
user.CoveragePlot <- function(scATAC.object, interest.genes, interested.conns, ccans, gene_model, idents = NULL, extend.upstream = 1000, extend.downstream = 1000, coaccess_cutoff = 0.2, signac = T, cicero = T, heights = NULL){
    require(GenomicRanges)
    res <- sapply(interest.genes, function(x){
        cat("Plot coverage: ", x, "\n")
        idx1 <- which(interested.conns$peak1_nearestGene == x)
        idx2 <- which(interested.conns$peak2_nearestGene == x)
        idx <- unique(idx1, idx2)
        conns <- data.frame(Peak1 = gsub("_", "-", interested.conns$Peak1[idx]),
                            Peak2 = gsub("_", "-", interested.conns$Peak2[idx]),
                            coaccess = as.numeric(interested.conns$coaccess[idx]))
        region.highlights <- c(conns$Peak1, conns$Peak2)
        region.highlights <- StringToGRanges(region.highlights)
        start.loc <- min(start(region.highlights))
        end.loc <- max(end(region.highlights))
        region.highlights$color <- rep(c("#ffa500", "#c71585"), each = nrow(conns)) #第一个颜色标记promoter

        links <- ConnectionsToLinks(conns = conns, ccans = ccans)
        Links(scATAC.object) <- links

        plot.regions <- paste0(as.character(seqnames(region.highlights))[1], "-", start.loc, "-", end.loc)
        if(signac){
            if(is.null(idents)){
                p <- CoveragePlot(scATAC.object, 
                    region = plot.regions,
                    region.highlight = region.highlights, 
                    extend.upstream = extend.upstream, extend.downstream = extend.downstream,
                    heights = heights)
            }else{
                p <- CoveragePlot(scATAC.object, 
                    region = plot.regions,
                    idents = idents,
                    region.highlight = region.highlights, 
                    extend.upstream = extend.upstream, extend.downstream = extend.downstream,
                    heights = heights)                
            }
            print(p)
        }

        #cicero plot
        if(cicero){
            require(cicero)
            require(monocle3)
            p <- plot_connections(conns,
                as.character(seqnames(region.highlights))[1], start.loc, end.loc,
                gene_model = gene_model,
                viewpoint = interested.conns$Peak1[idx[1]],
                coaccess_cutoff = coaccess_cutoff, 
                connection_width = .5, include_axis_track = F,
                collapseTranscripts = "longest")
            print(p)
        }
        return(conns)
    })
    return(res)
}