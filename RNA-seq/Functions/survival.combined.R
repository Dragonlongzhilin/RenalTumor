source(file = "/home/longzhilin/Analysis_Code/code/analysis.diff.survival.TCGA.R")
DESeq2.normalized_counts <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/DESeq2.normalized_counts.rds")
DESeq2.normalized_counts <- log2(DESeq2.normalized_counts+1)
DESeq2.result <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/DESeq2.result.rds")
clin.data <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/clin.data.rds")

survival.combined <- function(geneA, geneB, clin.data, DESeq2.result, DESeq2.normalized_counts, risk.table = T){
    gene.overlap <- intersect(c(geneA, geneB), rownames(DESeq2.result))
    diff.gene.pro.sig.gene <- DESeq2.result[gene.overlap,]
    interest.gene.sig <- data.frame(Gene = gene.overlap)
    sample.id <- colnames(DESeq2.normalized_counts)
    sample.id <- substr(sample.id, 1, 16)
    sample.type <- substr(sample.id, 14, 15)
    sample.id <- substr(sample.id, 1, 15)
    patient.overlap <- intersect(clin.data$Sample, sample.id)
    exp.data.process.pro <- DESeq2.normalized_counts[,match(patient.overlap, sample.id)]
    clin.info <- clin.data[match(patient.overlap, clin.data$Sample),]

    exp.A <- exp.data.process.pro[geneA,]
    exp.B <- exp.data.process.pro[geneB,]

    med.exp <- median(exp.A)
    high.group <- which(exp.A>med.exp)
    exp.A.label <- rep(paste0("Low ", geneA), length(exp.A))
    exp.A.label[high.group] <- paste0("High ", geneA)

    med.exp <- median(exp.B)
    high.group <- which(exp.B>med.exp)
    exp.B.label <- rep(paste0("Low ", geneB), length(exp.B))
    exp.B.label[high.group] <- paste0("High ", geneB)

    groups <- paste0(exp.A.label, " + ", exp.B.label)
    OS.data <- data.frame(Patient_ID = patient.overlap, event = clin.info$OS, time = clin.info$OS_time, sample.label = groups)
    DFS.data <- data.frame(Patient_ID = patient.overlap, event = clin.info$DFS, time = clin.info$DFS_time, sample.label = groups)
    p1 <- plot.surv(OS.data, risk.table = risk.table, HR = T, ylab = "Overall Survival", main = paste0(geneA, " + ", geneB), surv.median.line = "hv", xlab = "Time (Month)", colors = c("#D95319", "#F39B7F", "#3B6793","#4285F4"))
    print(p1)
    p2 <- plot.surv(DFS.data, risk.table = risk.table, HR = T, ylab = "Disease-Free Survival", main = paste0(geneA, " + ", geneB), surv.median.line = "hv", xlab = "Time (Month)", colors = c("#D95319", "#F39B7F", "#3B6793","#4285F4"))
    print(p2)
    return(list(OS.data = OS.data, DFS.data = DFS.data))
}

plot.surv <- function(clinical.data, upper.time = NULL, xscale = 1, xlab = "Time", median.time = TRUE, 
                      surv.median.line = "none", HR = FALSE, risk.table = TRUE, pval = TRUE, 
                      conf.int = FALSE, main = NULL, ylab = "Survival probability", colors = c("#D95319", "#F39B7F", "#3B6793","#4285F4")) {
  
#Load related R packages
  require(survival)
  require(survminer)
  require(RColorBrewer)
  require(gridExtra)

  #Determine the unit of event type and time
  # survival.event <- survival.event[1];
  # unit.xlabel <- unit.xlabel[1];

  #If upper.time is set, the samples whose survival time exceeds upper.time will be removed
  if (!is.null(upper.time)) clinical.data <- clinical.data[clinical.data$time <= upper.time,]

  #set color
  if (!is.factor(clinical.data$sample.label)) 
    clinical.data$sample.label <- as.factor(clinical.data$sample.label)
  
  t.name <- levels(clinical.data$sample.label)
  
if (length(t.name)> 6) stop("Sample grouping>6, exceeding the function acceptance range")
  t.col <- colors[1:length(t.name)]
  
# Construct a living object
   km.curves <- survfit(Surv(time, event)~sample.label, data=clinical.data)
  
   # Calculate HR value and 95% CI
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
  
   # Implicit function: conversion of survival time unit
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
                       data = clinical.data,    # data used to fit survival curves.
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
