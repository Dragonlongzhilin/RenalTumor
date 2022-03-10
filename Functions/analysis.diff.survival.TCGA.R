analysis.diff.survival.TCGA <- function(interest.gene, diff.gene.pro, exp.data.process, clin.data, group.by = "median", EnhancedVolcano.plot = T, Box.plot = T, main = NULL, meta.signature = F, single.signature = T){
  require(ggpubr)

  # require interest.gene to be data.frame
  if(EnhancedVolcano.plot){
    require(EnhancedVolcano)
    gene.overlap <- intersect(rownames(diff.gene.pro), rownames(interest.gene))
    diff.gene.pro.sig.gene <- diff.gene.pro[gene.overlap,]
    interest.gene.sig <- interest.gene[gene.overlap,]

    up.gene <- which(diff.gene.pro.sig.gene$log2FoldChange>0)
    if(length(up.gene)>0){
      down.gene <- setdiff(1:nrow(diff.gene.pro.sig.gene), up.gene)
      interest.gene.sig$label <- rep("down", nrow(diff.gene.pro.sig.gene))
      interest.gene.sig$label[up.gene] <- "up"
    }else{
      down.gene <- 1:nrow(diff.gene.pro.sig.gene)
      interest.gene.sig$label <- rep("down", nrow(diff.gene.pro.sig.gene))
    }
    
    interest.gene.sig$TCGA_log2FoldChange <- diff.gene.pro.sig.gene$log2FoldChange
    interest.gene.sig$TCGA_padj <- diff.gene.pro.sig.gene$padj
    
    p <- EnhancedVolcano(diff.gene.pro.sig.gene,
                  lab = rownames(diff.gene.pro.sig.gene),
                  selectLab = rownames(diff.gene.pro.sig.gene),
                  x = 'log2FoldChange',
                  y = 'padj', pCutoff = 0.05, drawConnectors = TRUE, FCcutoff = 1,
                  widthConnectors = 0.2,
                  title = "Tumor VS Normal", xlab = "", ylab = "-Log10(padj)")
    print(p)

    if(length(gene.overlap) < nrow(interest.gene)){
      a <- setdiff(rownames(interest.gene), gene.overlap)
      cat("gene ", paste(a, sep = ";"), " not in TCGA!\n")
    }else{
      cat("All gene in TCGA data!\n")
    }
  }else{ 
    #Require interest.gene to be a character vector
    gene.overlap <- intersect(interest.gene, rownames(diff.gene.pro))
    diff.gene.pro.sig.gene <- diff.gene.pro[gene.overlap,]
    interest.gene.sig <- data.frame(Gene = gene.overlap)
    if(length(gene.overlap) < length(interest.gene)){
      a <- setdiff(interest.gene, gene.overlap)
      cat("gene ", paste(a, sep = ";"), " not in TCGA!\n")
    }else{
      cat("All gene in TCGA data!\n")
    }
  }

	##Draw the expression barplot of each gene
	exp.data.interest <- exp.data.process[gene.overlap,]
	sample.id <- colnames(exp.data.process)
	sample.id <- substr(sample.id, 1, 16) 
	sample.type <- substr(sample.id, 14, 15)

	normal.index <- which(sample.type == "11") #72
	cancer.index <- which(sample.type == "01") #533
	patient.type <- rep("Tumor", ncol(exp.data.process))
	patient.type[normal.index] <- "Normal"	

  if(Box.plot){
  cat("Normal sample: ", length(normal.index), "\n")
  cat("Tumor sample: ", length(cancer.index), "\n")
    box.res <- sapply(gene.overlap, function(x){
      interest.gene.matrix <- exp.data.process[x,]
      genes <- data.frame(type = patient.type, expressionValue = interest.gene.matrix)
      p <- ggboxplot(genes, x = "type", y = "expressionValue", ylab = "Normalized Count", color = "type", palette = "jco", add = "jitter", title = x) + stat_compare_means()
      print(p)
	  })
  }		

	##对应表达和生存数据
	sample.id <- substr(sample.id, 1, 15)
	patient.overlap <- intersect(clin.data$Sample, sample.id) # 533 tumor samples
	exp.data.process.pro <- exp.data.process[,match(patient.overlap, sample.id)]
	clin.info <- clin.data[match(patient.overlap, clin.data$Sample),]
	
  #all gene as the signature
  if(meta.signature){
    interest.matrix <- exp.data.process.pro[gene.overlap,]
    signature.score <- colSums(interest.matrix)/nrow(interest.matrix)
    median.score <- median(signature.score)
    high.group <- which(signature.score > median.score)
    sample.label <- rep("Low group", length(signature.score))
    sample.label[high.group] <- "High group"
    OS.data <- data.frame(Patient_ID = patient.overlap, event = clin.info$OS, time = clin.info$OS_time, sample.label = sample.label)
    DFS.data <- data.frame(Patient_ID = patient.overlap, event = clin.info$DFS, time = clin.info$DFS_time, sample.label = sample.label)
    p1 <- plot.surv(OS.data, risk.table = T, HR = T, ylab = "Overall Survival", main = main, surv.median.line = "hv", xlab = "Time (Month)")
    print(p1)
    p2 <- plot.surv(DFS.data, risk.table = T, HR = T, ylab = "Disease-Free Survival", main = main, surv.median.line = "hv", xlab = "Time (Month)")
    print(p2)
  }

	## single gene model
  if(single.signature){
    suv.res <- sapply(gene.overlap, function(x){
      gene.expression <- as.numeric(exp.data.process.pro[x,])	

      if(group.by == "top33"){
        quantile.value <- quantile(gene.expression, seq(0,1,0.33))
        top33 <- as.numeric(quantile.value[3])
        bottom33 <- as.numeric(quantile.value[2])
        index1 <- which(gene.expression >= top33)
        index2 <- which(gene.expression <= bottom33)
        sample.label <- c(rep("Low group", length(index1)), rep("High group", length(index2)))
        clin.info <- clin.info[c(index2,index1),]
        OS.data <- data.frame(Patient_ID = patient.overlap[c(index2,index1)], event = clin.info$OS, time = clin.info$OS_time, sample.label = sample.label)
        DFS.data <- data.frame(Patient_ID = patient.overlap[c(index2,index1)], event = clin.info$DFS, time = clin.info$DFS_time, sample.label = sample.label)
      }else{
        # median as threshold
        med.exp <- median(gene.expression)
        high.group <- which(gene.expression>med.exp)
        sample.label <- rep("Low group", length(gene.expression))
        sample.label[high.group] <- "High group"
        OS.data <- data.frame(Patient_ID = patient.overlap, event = clin.info$OS, time = clin.info$OS_time, sample.label = sample.label)
        DFS.data <- data.frame(Patient_ID = patient.overlap, event = clin.info$DFS, time = clin.info$DFS_time, sample.label = sample.label)
      }
      p1 <- plot.surv(OS.data, risk.table = T, HR = T, ylab = "Overall Survival", main = x, surv.median.line = "hv", xlab = "Time (Month)")
      print(p1)
      p2 <- plot.surv(DFS.data, risk.table = T, HR = T, ylab = "Disease-Free Survival", main = x, surv.median.line = "hv", xlab = "Time (Month)")
      print(p2)
      
      OS.obj <- coxph(Surv(time, event)~sample.label, data=OS.data)
      DFS.obj <- coxph(Surv(time, event)~sample.label, data=DFS.data)
      return(c(summary(OS.obj)$logtest[3], summary(DFS.obj)$logtest[3]))
    })
    interest.gene.sig$OS_logrank_p <- suv.res[1,]
    interest.gene.sig$DFS_logrank_p <- suv.res[2,] 
  }

  return(interest.gene.sig)
}

plot.surv <- function(clinical.data, upper.time = NULL, xscale = 1, xlab = "Time", median.time = TRUE, 
                      surv.median.line = "none", HR = FALSE, risk.table = TRUE, pval = TRUE, 
                      conf.int = FALSE, main = NULL, ylab = "Survival probability", colors = c("#D95319", "#3B6793","#EA4335","#4285F4","#34A853","#000000")) {
  
  #load R package
  require(survival)
  require(survminer)
  require(RColorBrewer)
  require(gridExtra)
  
  # Determine the type of event and the unit of time
  # survival.event <- survival.event[1];
  # unit.xlabel <- unit.xlabel[1];
  
  # If upper.time is set, the samples whose survival time exceeds upper.time will be removed
  if (!is.null(upper.time)) clinical.data <- clinical.data[clinical.data$time <= upper.time,]
  
  # #date format
  # xSL <- data.frame(xScale=c(1,7,30,365.25),xLab=c("Days","Weeks","Months","Years"), stringsAsFactors=FALSE)
  # switch(unit.xlabel, year={xScale <- 365.25;}, month={xScale <- 30;}, week={xScale <- 7;}, day={xScale <- 1})
  # xLab <- xSL[which(xSL[,1]==xScale),2];
  
  # color
  if (!is.factor(clinical.data$sample.label)) 
    clinical.data$sample.label <- as.factor(clinical.data$sample.label)
  
  t.name <- levels(clinical.data$sample.label)
  
  if (length(t.name) > 6) stop("Sample grouping >6, exceeding the function acceptance range!")
  t.col <- colors[1:length(t.name)]
  
  # 构造生存对象
  km.curves <- survfit(Surv(time, event)~sample.label, data=clinical.data)
  
  # HR and 95%CI
  if (length(t.name) == 2) {
    if (HR) {
      cox.obj <- coxph(Surv(time, event)~sample.label, data=clinical.data)
      tmp <- summary(cox.obj)
      HRs <- round(tmp$coefficients[ ,2], digits = 2)
      HR.confint.lower <- round(tmp$conf.int[,"lower .95"], 2)
      HR.confint.upper <- round(tmp$conf.int[,"upper .95"], 2)
      HRs <- paste0(HRs, " (", HR.confint.lower, "-", HR.confint.upper, ")")	  
    }
  }
  
  # Construct the legend display text in the survival image
  legend.content <- substr(names(km.curves$strata), start = 14, stop = 1000)
  
  # x-axis scale unit conversion
  if (is.numeric(xscale) | (xscale %in% c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"))) {
    xscale = xscale
  } else {
    stop('xscale should be numeric or one of c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m").')
  }
  
  # conversion of survival time units
  .format_xticklabels <- function(labels, xscale){
    # 1 year = 365.25 days
    # 1 month = 365.25/12 = 30.4375 days
    if (is.numeric(xscale)) xtrans <- 1/xscale
    else
      xtrans <- switch(xscale,
                       d_m = 12/365.25,
                       d_y = 1/365.25,
                       m_d = 365.25/12,
                       m_y = 1/12,
                       y_d = 365.25,
                       y_m = 12,
                       1
      )
    round(labels*xtrans,2)
  }
  
  # Add the median survival time and its 95% CI in the figure and place it in the subtitle position
  subtitle <- NULL
  if (median.time) {
    if (is.numeric(xscale)) {
      median.km.obj = km.curves
    } else if (xscale %in% c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m")) {
      clinical.data$time <- .format_xticklabels(labels = clinical.data$time, xscale = xscale)
      median.km.obj <- survfit(Surv(time, event)~sample.label, data=clinical.data)
    }
    survival.time.info <- NULL
    survival.time.info <- rbind(survival.time.info, summary(median.km.obj)$table)
    median.survival <- round(survival.time.info[!duplicated(survival.time.info[,7:9]),7:9], digits = 2) # 注意：这里取得的置信区间上界可能为NA
    if (length(levels(clinical.data$sample.label)) == 1) {
      tmp1 <- levels(clinical.data$sample.label)
    } else {
      tmp1 <- do.call(rbind,strsplit(rownames(summary(median.km.obj)$table), split = "="))[,2]
    }
    tmp2 <- paste(median.survival[,1], "(", median.survival[,2], "-", median.survival[,3], ")")
    subtitle <- paste(tmp1, tmp2, sep = ":", collapse = "\n")
  }
  
  # ggsurvplot
  ggsurv <- ggsurvplot(km.curves,               # survfit object with calculated statistics.
                       data = clinical.data,             # data used to fit survival curves.
                       palette = t.col,
                       
                       risk.table = risk.table,       # show risk table.
                       pval = pval,             # show p-value of log-rank test.
                       surv.median.line = surv.median.line,  # add the median survival pointer.
                       title = main,     #main title
                       subtitle = subtitle, #sub title
                       font.main = 15,                 
                       xlab = xlab,   # customize X axis label.
                       ylab = ylab,   # customize Y axis label
                       xscale = xscale,
                       
                       
                       #legend
                       legend.title = "", 
                       legend.labs = legend.content, 
                       legend = c(0.8,0.9), 
                       font.legend = 9,     
                       
                       #risk table
                       tables.theme = theme_cleantable(),
                       risk.table.title = "No. at risk:",
                       risk.table.y.text.col = T, 
                       risk.table.y.text = FALSE, 
                       tables.height = 0.15,      
                       risk.table.fontsize = 3  
  )
  if (length(t.name) == 2) {
    if (HR) 
      ggsurv$plot <- ggsurv$plot + ggplot2::annotate("text", x = max(km.curves$time)/12, 
                                                     y = 0.15, size = 5, label = paste("HR=", HRs))
  } 

  ggsurv$plot <- ggsurv$plot + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(size = 10), 
                                     plot.margin = unit(c(5.5, 5.5, 5.5, 50), "points"))

  ggsurv$table <- ggsurv$table + theme(plot.title = element_text(hjust = -0.04),
                                       plot.margin = unit(c(5.5, 5.5, 5.5, 50), "points"))
  

  if(length(t.name) > 2) {
    # pairwise: log rank P value
    res <- pairwise_survdiff(Surv(time, event)~sample.label, data=clinical.data);
    pairwise.pvalue <- round(res$p.value, digits = 4);
    pairwise.pvalue[which(pairwise.pvalue < 0.0001)] <- "<0.0001";
    pairwise.pvalue[is.na(pairwise.pvalue)] <- "-"
    
    # add table
    tt <- ttheme_minimal(core = list(fg_params = list(col = "black"),bg_params = list(fill = NA, col = "black")),
                         colhead = list(fg_params = list(col = NA),bg_params = list(fill = t.col, col = "black")),
                         rowhead = list(fg_params = list(col = NA, hjust = 1),bg_params = list(fill = c("white",t.col[-1]), col = "black"))
    )
    pairwise.table <- tableGrob(pairwise.pvalue, theme = tt)
    ggsurv <- ggarrange(ggarrange(ggsurv$plot, ggsurv$table, nrow=2, heights=c(2,0.5)),
                        pairwise.table, nrow=2, heights = c(2,0.5),
                        labels = c("","p from pairwise comparisons"),
                        hjust = 0, font.label = list(size = 15, face = "plain"))
  }
  
  ggsurv	
}
