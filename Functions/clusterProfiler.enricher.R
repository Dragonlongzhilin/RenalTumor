#### 使用clusterProfiler进行通路富集分析

#' @param gene 字符向量，感兴趣的基因集合, 如"FUT7"         "AIRE"         "IFI35"
#' @param geneType 字符型向量，指明输入的geneList中name属性类型，即基因id的类型
#' @param db.type 字符型向量，指明GSEA针对的数据库
#' @param saveDir gsea输入目录
#' @param pAdjustMethod, pvalueCutoff p值矫正方法和阈值
#' @param MsigDB.category 指定特定的MsigDB中的数据库进行GSEA分析
#' MsigDB.pathway, msigdbr(species = "Homo sapiens", category = i)， 指定i为H, C1, C2, C3, C4, C5, C6, C7, C8
# 通过gs_cat gs_subcat，来选择库
#' @param GO.ont 指定GSEA分析GO数据库时使用的条目类型，支持"BP", "MF", and "CC"，以及"ALL"

#' 功能富集分析的原理是基于hypergeometric test
#' clusterProfiler主要基于EntrezID做分析
#' 分析过程中的universe，即背景基因的数目，是跟你选择的db中的subcategory有关系的，而不是固定不变的
cluterProfiler.enricher <- function(gene, geneType = c("ENTREZID", "SYMBOL", "ENSEMBL"),
                                    db.type = c("MsigDB", "GO", "KEGG"),
                                    saveDir, title,
                                    pAdjustMethod = "BH",
                                    qvalueCutoff = 0.2,
                                    pvalueCutoff = 0.05,
                                    MsigDB.category = list(H = c("All"), C2 = c("BIOCARTA", "KEGG", "REACTOME"), C5 = c("BP")),
                                    GO.ont = "BP", simplify = T,
                                    minGSSize = 10, maxGSSize = 500,
                                    showCategory = 20
                                    ){
  require(enrichplot)
  require(clusterProfiler)
  require(tidyverse)

  
  geneType <- match.arg(geneType)
  db.type <- match.arg(db.type)  
  # 如果基因名字类型不是SYMBOL，则需要进行ID转换
  if(geneType == "ENTREZID"){
    cat("不需要转换\n")
    geneList <- gene
  }
  if(geneType %in% c("SYMBOL", "ENSEMBL")){
    genes <- IDConvert(genes = gene, fromType = geneType, toType = "ENTREZID")
    cat("转换前基因数目：", length(gene), "\n")
    index <- match(gene, genes[,1])
    geneList <- as.character(genes[na.omit(index),2])
    cat("转换后基因数目：", length(geneList), "\n")    
  }
  
  ## 第一种富集策略：MsigDB基因集合
  if(db.type == "MsigDB"){
      library(msigdbr)
    
      # 通路提取
      MsigDB.pathway <- tibble(gs_name = NULL, entrez_gene = NULL)
      
      # 将MsigDB所有基因作为背景基因
      #MsigDB.gene <- unique(msigdbr(species = "Homo sapiens") %>% dplyr::select(gs_name,entrez_gene))
      
      if(class(MsigDB.category) == "character"){
        MsigDB.pathway <- msigdbr(species = "Homo sapiens") %>% dplyr::select(gs_name, entrez_gene)
      }else{
        for(i in names(MsigDB.category)){
          if(MsigDB.category[[i]][1] == "All"){
            MsigDB.pathway <- rbind(MsigDB.pathway, msigdbr(species = "Homo sapiens", category = i) %>% dplyr::select(gs_name, entrez_gene))
          }else{
            for(j in MsigDB.category[[i]]){
              MsigDB.pathway <- rbind(MsigDB.pathway, msigdbr(species = "Homo sapiens", category = i, subcategory = j) %>% dplyr::select(gs_name, entrez_gene))
            }
          }
        }
      }
      em.res <- enricher(geneList, TERM2GENE = MsigDB.pathway, qvalueCutoff = qvalueCutoff, pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize)
  }
  
  if(db.type == "GO"){
    em.res <- enrichGO(gene = geneList,
                    OrgDb = org.Hs.eg.db,
                    ont = GO.ont,minGSSize = minGSSize, maxGSSize = maxGSSize,
                    pAdjustMethod = pAdjustMethod,
                    pvalueCutoff  = pvalueCutoff,
                    qvalueCutoff  = qvalueCutoff)
    #基于simplify method, we can remove redundant terms
    if(simplify){
      em.res <- simplify(em.res, cutoff = 0.7, by="p.adjust", select_fun = min)
    }
  }
  
  if(db.type == "KEGG"){
    em.res <- enrichKEGG(gene = geneList, organism = "hsa", keyType = "kegg",
                       pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize,
                       pvalueCutoff  = pvalueCutoff,
                       qvalueCutoff  = qvalueCutoff)
  }
  #转换为geneSymbol
  em.res.geneSymbol <- setReadable(em.res, org.Hs.eg.db, keyType = "ENTREZID")
    
  # if(nrow(em.res@result) > 0){
  #   write.csv(em.res.geneSymbol, file = file.path(saveDir, paste0(title, "_", "enricherResult.csv")))
  #   tryCatch({
  #     p1 <- barplot(em.res.geneSymbol, showCategory = showCategory, main = title) + theme(axis.text.x = element_text(size=4), axis.text.y = element_text(size=4))
  #     print(p1)
  #     p2 <- dotplot(em.res.geneSymbol, showCategory = showCategory, title = title) + theme(axis.text.x = element_text(size=4), axis.text.y = element_text(size=4))
  #     print(p2)
  #     p3 <- heatplot(em.res.geneSymbol, showCategory = showCategory) + theme(axis.text.x = element_text(size=4), axis.text.y = element_text(size=4))
  #     print(p3)
  #   }, warning = function(w){
  #     # 这里是出现warning状态时，应该怎么做，可以用print打印出来，可以执行其它命令
  #     cat("Warning in plot!\n")
  #   }, error = function(e){
  #     # 这里时出现Error状态时，应该怎么做，可以用print打印出来，也可以执行其它命令
  #     cat("Error in plot!\n")
  #   },finally = {
  #     # 这是运行正常时，应该怎么做，可以用print打印出来，也可以执行其它命令
  #     cat("Plot complete!\n")
  #   })
  # }
  return(em.res = list(em.res.entrezid = em.res, em.res.genesymbol = em.res.geneSymbol))
}


##### ID转换
#author: Zhilin long
#time: 2021.1.19

#' @author Zhilin Long
#' @param gene 字符型向量，需要转换的ID向量
#' @param fromType 字符向量，初始的ID类型
#' @param toType 字符向量，需要转换的ID类型
#' 
#' @return data.frame 已经去除了重复和NA，可直接使用
IDConvert <- function(genes, 
                      method = c("clusterProfiler", "org.Hs.eg.db"),
                      fromType = c("ENTREZID", "SYMBOL", "ENSEMBL"), 
                      toType = c("ENTREZID", "SYMBOL", "ENSEMBL")
                      ){
    require(org.Hs.eg.db)
    require(clusterProfiler)

    #规范参数输入
    method <- match.arg(method)
    fromType <- match.arg(fromType)
    toType <- match.arg(toType)

    if(method == "clusterProfiler"){
        gene.df <- bitr(genes, fromType = fromType, toType = toType, OrgDb = "org.Hs.eg.db")
        #默认会过滤没有map到的id，即条目会减少
    }else{
        gene.df <- mapIds(org.Hs.eg.db, 
                          keys = genes, 
                          column = toType, 
                          keytype = fromType,
                          multiVals="first")
        gene.df <- data.frame(ENSEMBL = names(gene.df), SYMBOL = gene.df)
        colnames(gene.df) <- c(fromType, toType)
        #没有map，则为NA，即条目不变
        gene.df <- na.omit(gene.df) #去除NA的行
    }

    #去除重复, 
    gene.df <- gene.df[!duplicated(gene.df[,2]),]
    gene.df <- gene.df[!duplicated(gene.df[,1]),]
    return(gene.df)
}