#### plot ligand-receptor interaction
#ref：https://github.com/Teichlab/cellphonedb/blob/master/cellphonedb/src/plotters/R/plot_dot_by_column_name.R

#' @param selected_rows: specific ligand-receport pairs. note, Cannot be a factor variable
#' @param selected_columns: specific cell types
#' @param size: font size
#' @param max.dot.size: set max dot size
#' @param breaks: set the position of ticks, require the limit parameter
#' @param limit: set the range of legend bar

dot_plot = function(selected_rows = NULL,
                    selected_columns = NULL,
                    size = 8, max.dot.size = 2, breaks = c(-2, -1, 0, 1, 2), limit = c(-2, 2),
                    means_path = './means.txt',
                    pvalues_path = './pvalues.txt',
                    means_separator = '\t',
                    pvalues_separator = '\t',
                    filter.columns = c(1:12),
                    order.column = NULL,
                    legend.key.size = unit(0.2, "inches")
){
  require(ggplot2)

  all_pval <- read.table(pvalues_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
  all_means <- read.table(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)

  intr_pairs <- all_pval$interacting_pair
  ligand <- data.frame(pair = intr_pairs, ligand = all_pval$ligand)
  all_pval <- all_pval[,-filter.columns]
  all_means <- all_means[,-filter.columns]

  if(is.null(selected_rows)){
    selected_rows <- intr_pairs
  }

  if(is.null(selected_columns)){
    selected_columns <- colnames(all_pval)
  }

  sel_pval <- all_pval[match(selected_rows, intr_pairs), selected_columns] #通过select_columns和select_rows筛选出显著的互作
  sel_means <- all_means[match(selected_rows, intr_pairs), selected_columns]

  df_names <- expand.grid(selected_rows, selected_columns)
  pval <- unlist(sel_pval)
  pval[pval==0] <- 0.0009
  pval <- -log10(pval)
  plot.data <- cbind(df_names, pval)
  pr <- unlist(as.data.frame(sel_means))
  pr[pr==0] <- 1 #Padded the mean of 0 to 1, to ensure that the log2 conversion will not go wrong
  plot.data <- cbind(plot.data,log2(pr))
  colnames(plot.data) <- c('pair', 'clusters', 'pvalue', 'mean') 

 #Based on the ligand information, adjust the interaction information, distinguish between ligand and receptor
  plot.data <- apply(plot.data, 1, function(x){
    idx <- which(ligand$pair == x[1])
    if(ligand[idx,2] == "Receptor"){
      a <- unlist(strsplit(x[1], "_"))
      x[1] <- paste0(a[2], "_", a[1])

      a <- unlist(strsplit(x[2], "\\|"))
      x[2] <- paste0(a[2], "|", a[1])
    }
    return(x)
  })
  plot.data <- as.data.frame(t(plot.data))
  if(is.null(order.column)){
    plot.data$clusters <- factor(plot.data$clusters, levels = selected_columns)
  }else{
    idx <- which(plot.data$clusters %in% order.column)
    plot.data <- plot.data[idx,]
    plot.data$clusters <- factor(plot.data$clusters, levels = order.column)
  }

  plot.data$pair <- gsub("_", "-", plot.data$pair)
  plot.data$pvalue <- as.numeric(plot.data$pvalue)
  plot.data$mean <- as.numeric(plot.data$mean)

  my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
#  my_palette <- c(colorRampPalette(colors = c("black", "blue"),  alpha=TRUE)(length(bk)/2), colorRampPalette(colors = c("yellow", "red"), alpha=TRUE)(length(bk)/2))

  p <- ggplot(plot.data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=pvalue,color=mean)) + 
  scale_color_gradientn('Log2 mean\n(Molecule 1, Molecule 2)', colors=my_palette, breaks = breaks, limit = limit) +
  scale_size_continuous(name = '-Log10(p value)', range = c(0, max.dot.size))+
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=size, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=size, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"),
        legend.text = element_text(size = size),
        legend.title = element_text(size = size*1.25),
        legend.key.size = legend.key.size)
  return(p)
#   if (output_extension == '.pdf') {
#       ggsave(filename, width = width, height = height, device = cairo_pdf, limitsize=F)
#   }
#   else {
#       ggsave(filename, width = width, height = height, limitsize=F)
#   }
}