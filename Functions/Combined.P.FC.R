####整合p和FC
#reference: A novel signiﬁcance score for gene selection and ranking. 2012 Bioinformatics

#' @author Zhilin Long
#' @param data data.frame, 第一列为FC，第二列为p value
#' @param log2FC logical；T表示FC已经经过了log2处理
#' @param log10P logical; T表示P已经经过了-log10处理
#' @return data.frame，增加了一列pi值, pi越大表示统计学意义越强，一般用于GSEA基因排序

Combined.P.FC <- function(data, log2FC = T, log10P = T, min.p = 0.00001){

    data <- na.omit(data)
    if(log2FC){
        LFC <- abs(data[,1])
    }else{
        LFC <- abs(log2(as.numeric(data[,1])))
    }

    ####处理p为0的情况：
	#若数据中最小值大于0.00001，则用0.00001替代0，
	#若数据中最小值小于0.00001，则在最小值基础上缩小0.1倍
    data.pro <- data[which(data[,2]>0),]
    min.value <- min(data.pro[,2])
    if(min.value < min.p){
        min.value <- min.value*0.1
    }else{
        min.value <- min.p
    }

    index <- which(data[,2]==0)
    data[index,2] <- min.value

    if(log10P){
        LP <- data[,2]
    }else{
        LP <- -log10(as.numeric(data[,2]))
    }

    pi <- LFC*LP
    data$pi <- pi
    return(data)
}