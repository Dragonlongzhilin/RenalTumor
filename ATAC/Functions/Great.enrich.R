
wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}
Great.enrich <- function(gr, bg, title, species = "hg38", fdr.threshold = 0.05, Type = "GO Biological Process", n_terms = 20){
    require(rGREAT)
    require(tidyverse)
    require(cowplot)
    theme_set(theme_cowplot())
    set.seed(101)
    job <- submitGreatJob(gr = gr, bg = bg, species = species)
    tb <- getEnrichmentTables(job)
    region.graph <- plotRegionGeneAssociationGraphs(job)

    # "GO Molecular Function" "GO Biological Process" "GO Cellular Component"
    sig.tb <- lapply(names(tb), function(x){
        enr <- tb[[x]]
        enr$Term <- rep(x, nrow(enr))
        sig <- which(enr$Hyper_Adjp_BH < fdr.threshold)
        if(length(sig>0)){
            return(enr[sig,])
        }else{
            return(data.frame())
        }
    })
    rGREAT.res <- Reduce(rbind, sig.tb)
    if(nrow(rGREAT.res) > 0){
        if(!is.null(Type)){
            rGREAT.res <- rGREAT.res[which(rGREAT.res$Term == Type),]
        }

        if(nrow(rGREAT.res) > 0){
            if(nrow(rGREAT.res) >= n_terms){
            rGREAT.res <- rGREAT.res[1:n_terms,]
            }
            rGREAT.res$wrap <- wrapText(rGREAT.res$name, 45)
            p <- ggplot(rGREAT.res, aes(y = log(Hyper_Fold_Enrichment), x = reorder(wrap, log(Hyper_Fold_Enrichment)))) +
                geom_bar(stat='identity', aes(color = name, fill = name)) +
                xlab('') + theme(legend.position="none") + ggtitle(title) + coord_flip()
            print(p)
        }
    }
    return(rGREAT.res)
}
