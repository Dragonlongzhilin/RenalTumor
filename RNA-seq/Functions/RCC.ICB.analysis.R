source(file = "/home/longzhilin/Analysis_Code/SurvivalAnalysis/Cox.function.R")

RCC.icb.analysis <- function(signature.list, expresionMatrix, clincal.info){
    library(ggpubr)    
    signature.list.score <- lapply(names(signature.list), function(y){

        genes <- signature.list[[y]]
        programName <- y
        index <- match(genes, rownames(expresionMatrix))
        if(length(na.omit(index)) < length(index)){
            cat(y, " lack ", genes[which(is.na(index))], " gene\n")
        }
        signature.score <- colMeans(expresionMatrix[na.omit(index),])

        ICB.treatment.idx <- which(clincal.info$Arm == "NIVOLUMAB")
        signature.score.icb <- signature.score[ICB.treatment.idx]
        clincal.info.icb <- clincal.info[ICB.treatment.idx,]
        plot.data.icb <- data.frame(score = signature.score.icb, Cohort = clincal.info.icb$Cohort, Benefit = clincal.info.icb$Benefit)
        p <- ggboxplot(plot.data.icb, x = "Cohort", y = "score", color = "Benefit", xlab = "", ylab = "Signature score", title = paste0("ICB threapy", "-", programName), add = "jitter") + stat_compare_means(aes(group = Benefit), label = "p.format")
        print(p)


        index <- which(plot.data.icb$Benefit=="ICB")
        plot.data.icb.CB.VS.ICB <- plot.data.icb[-index,]
        p <- ggboxplot(plot.data.icb.CB.VS.ICB, x = "Cohort", y = "score", color = "Benefit", xlab = "", ylab = "Signature score", title = paste0("ICB threapy", "-", programName), add = "jitter") + stat_compare_means(aes(group = Benefit), label = "p.format")
        print(p)    

        mTOR.treatment.idx <- which(clincal.info$Arm == "EVEROLIMUS")
        signature.score.mTOR <- signature.score[mTOR.treatment.idx]
        clincal.info.mTOR <- clincal.info[mTOR.treatment.idx,]
        plot.data.mTOR <- data.frame(score = signature.score.mTOR, Cohort = clincal.info.mTOR$Cohort, Benefit = clincal.info.mTOR$Benefit)
        p <- ggboxplot(plot.data.mTOR, x = "Cohort", y = "score", color = "Benefit", xlab = "", ylab = "Signature score", title = paste0("mTOR threapy", "-", programName), add = "jitter") + stat_compare_means(aes(group = Benefit), label = "p.format")
        print(p)

        index <- which(plot.data.mTOR$Benefit=="ICB")
        plot.data.mTOR.CB.VS.ICB <- plot.data.mTOR[-index,]
        p <- ggboxplot(plot.data.mTOR.CB.VS.ICB, x = "Cohort", y = "score", color = "Benefit", xlab = "", ylab = "Signature score", title = paste0("mTOR threapy", "-", programName), add = "jitter") + stat_compare_means(aes(group = Benefit), label = "p.format")
        print(p)        

        #survival plot
        #Based on median grouping
        survival.res <- sapply(c("CM-009", "CM-010", "CM-025"), function(x){

            index <- which(clincal.info.icb$Cohort == x)
            clincal.info.icb.cohort <- clincal.info.icb[index,]
            signature.score.icb.cohort <- signature.score.icb[index]

            mid.score.icb <- median(signature.score.icb.cohort)
            icb.label <- rep("High score", length(signature.score.icb.cohort))
            icb.label[which(signature.score.icb.cohort<=mid.score.icb)] <- "Low score"
            if(length(table(icb.label)) <= 1){
                cox.res <- "NA"
            }else{
                if(x == "CM-025"){
                    # different treatment

                    icb.clincal.OS <- data.frame(Patient_ID = clincal.info.icb.cohort$RNA_ID, event = clincal.info.icb.cohort$OS_CNSR, time = clincal.info.icb.cohort$OS, sample.label = icb.label)
                    p1 <- plot.surv(icb.clincal.OS, risk.table = T, HR = T, ylab = "Overall Survival", main = paste0("ICB ", programName, " ", x), surv.median.line = "hv", xlab = "Time (Month)")
                    print(p1)
                    icb.clincal.PFS <- data.frame(Patient_ID = clincal.info.icb.cohort$RNA_ID, event = clincal.info.icb.cohort$PFS_CNSR, time = clincal.info.icb.cohort$PFS, sample.label = icb.label)
                    p1 <- plot.surv(icb.clincal.PFS, risk.table = T, HR = T, ylab = "Progression-Free Survival", main = paste0("ICB ", programName, " ", x), surv.median.line = "hv", xlab = "Time (Month)")
                    print(p1)
                    
                    cox.clincal.data <- data.frame(Patient_ID = clincal.info.icb.cohort$RNA_ID, time = clincal.info.icb.cohort$OS, event = clincal.info.icb.cohort$OS_CNSR, sample.label = icb.label, Sex = clincal.info.icb.cohort$Sex, Age = clincal.info.icb.cohort$Age, Prior_Therapies = clincal.info.icb.cohort$Number_of_Prior_Therapies, Metastasis = clincal.info.icb.cohort$Tumor_Sample_Primary_or_Metastasis)
                    cox.clincal.data <- na.omit(cox.clincal.data)
                    icb.cox.OS <- Cox.function(time = cox.clincal.data$time, event = cox.clincal.data$event, clinical.data = cox.clincal.data)
                    cox.clincal.data <- data.frame(Patient_ID = clincal.info.icb.cohort$RNA_ID, time = clincal.info.icb.cohort$PFS, event = clincal.info.icb.cohort$PFS_CNSR, sample.label = icb.label, Sex = clincal.info.icb.cohort$Sex, Age = clincal.info.icb.cohort$Age, Prior_Therapies = clincal.info.icb.cohort$Number_of_Prior_Therapies, Metastasis = clincal.info.icb.cohort$Tumor_Sample_Primary_or_Metastasis)
                    cox.clincal.data <- na.omit(cox.clincal.data)
                    icb.cox.PFS <- Cox.function(time = cox.clincal.data$time, event = cox.clincal.data$event, clinical.data = cox.clincal.data)

                    ####mTOR 
                    mid.score.mTOR <- median(signature.score.mTOR)
                    mTOR.label <- rep("High score", length(signature.score.mTOR))
                    mTOR.label[which(signature.score.mTOR<=mid.score.mTOR)] <- "Low score"
                    mTOR.clincal.OS <- data.frame(Patient_ID = clincal.info.mTOR$RNA_ID, event = clincal.info.mTOR$OS_CNSR, time = clincal.info.mTOR$OS, sample.label = mTOR.label)

                    p1 <- plot.surv(mTOR.clincal.OS, risk.table = T, HR = T, ylab = "Overall Survival", main = paste0("mTOR ", programName, " ", x), surv.median.line = "hv", xlab = "Time (Month)")
                    print(p1)
                    mTOR.clincal.PFS <- data.frame(Patient_ID = clincal.info.mTOR$RNA_ID, event = clincal.info.mTOR$PFS_CNSR, time = clincal.info.mTOR$PFS, sample.label = mTOR.label)
                    p1 <- plot.surv(mTOR.clincal.PFS, risk.table = T, HR = T, ylab = "Progression-Free Survival", main = paste0("mTOR ", programName, " ", x), surv.median.line = "hv", xlab = "Time (Month)")
                    print(p1)

                    cox.clincal.data <- data.frame(Patient_ID = clincal.info.mTOR$RNA_ID, time = clincal.info.mTOR$OS, event = clincal.info.mTOR$OS_CNSR, sample.label = mTOR.label, Sex = clincal.info.mTOR$Sex, Age = clincal.info.mTOR$Age, Prior_Therapies = clincal.info.mTOR$Number_of_Prior_Therapies, Metastasis = clincal.info.mTOR$Tumor_Sample_Primary_or_Metastasis)
                    cox.clincal.data <- na.omit(cox.clincal.data)
                    mTOR.cox.OS <- Cox.function(time = cox.clincal.data$time, event = cox.clincal.data$event, clinical.data = cox.clincal.data)
                    cox.clincal.data <- data.frame(Patient_ID = clincal.info.mTOR$RNA_ID, time = clincal.info.mTOR$PFS, event = clincal.info.mTOR$PFS_CNSR, sample.label = mTOR.label, Sex = clincal.info.mTOR$Sex, Age = clincal.info.mTOR$Age, Prior_Therapies = clincal.info.mTOR$Number_of_Prior_Therapies, Metastasis = clincal.info.mTOR$Tumor_Sample_Primary_or_Metastasis)
                    cox.clincal.data <- na.omit(cox.clincal.data)
                    mTOR.cox.PFS <- Cox.function(time = cox.clincal.data$time, event = cox.clincal.data$event, clinical.data = cox.clincal.data)

                    cox.res <- list(icb.OS = icb.cox.OS, icb.PFS = icb.cox.PFS, mTOR.OS = mTOR.cox.OS, mTOR.PFS = mTOR.cox.PFS)
                }else{

                    icb.clincal.OS <- data.frame(Patient_ID = clincal.info.icb.cohort$RNA_ID, event = clincal.info.icb.cohort$OS_CNSR, time = clincal.info.icb.cohort$OS, sample.label = icb.label)
                    p1 <- plot.surv(icb.clincal.OS, risk.table = T, HR = T, ylab = "Overall Survival", main = paste0(programName, " ", x), surv.median.line = "hv", xlab = "Time (Month)")
                    print(p1)
                    icb.clincal.PFS <- data.frame(Patient_ID = clincal.info.icb.cohort$RNA_ID, event = clincal.info.icb.cohort$PFS_CNSR, time = clincal.info.icb.cohort$PFS, sample.label = icb.label)
                    p1 <- plot.surv(icb.clincal.PFS, risk.table = T, HR = T, ylab = "Progression-Free Survival", main = paste0(programName, " ", x), surv.median.line = "hv", xlab = "Time (Month)")
                    print(p1)

                    cox.clincal.data <- data.frame(Patient_ID = clincal.info.icb.cohort$RNA_ID, time = clincal.info.icb.cohort$OS, event = clincal.info.icb.cohort$OS_CNSR, sample.label = icb.label, Sex = clincal.info.icb.cohort$Sex, Age = clincal.info.icb.cohort$Age)
                    cox.clincal.data <- na.omit(cox.clincal.data)
                    icb.cox.OS <- Cox.function(time = cox.clincal.data$time, event = cox.clincal.data$event, clinical.data = cox.clincal.data)
                    cox.clincal.data <- data.frame(Patient_ID = clincal.info.icb.cohort$RNA_ID, time = clincal.info.icb.cohort$PFS, event = clincal.info.icb.cohort$PFS_CNSR, sample.label = icb.label, Sex = clincal.info.icb.cohort$Sex, Age = clincal.info.icb.cohort$Age)
                    cox.clincal.data <- na.omit(cox.clincal.data)
                    icb.cox.PFS <- Cox.function(time = cox.clincal.data$time, event = cox.clincal.data$event, clinical.data = cox.clincal.data)

                    cox.res <- list(icb.OS = icb.cox.OS, icb.PFS = icb.cox.PFS)
                }
            }     
            return(cox.res)
        })
        
    })
    names(signature.list.score) <- names(signature.list)
    return(signature.list.score)
}

plot.surv <- function(clinical.data, upper.time = NULL, xscale = 1, xlab = "Time", median.time = TRUE, 
                      surv.median.line = "none", HR = FALSE, risk.table = TRUE, pval = TRUE, 
                      conf.int = FALSE, main = NULL, ylab = "Survival probability", colors = c("#D95319", "#3B6793","#EA4335","#4285F4","#34A853","#000000")) {
  
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
