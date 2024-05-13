#整理各数据集中三个干性量化算法的得分，并对其加权平均，保存结果

GSE <- c("GSE110009","GSE178318","GSE183916","GSE221575","GSE225857")
weight<-read.csv("E:/2021020529/课题/细胞干性课题/result/各算法权重/weight(H1,H9).csv",header = T,row.names = 1)
for (i in GSE) {
  print(i)
  scRNA <- readRDS(paste0("E:/2021020529/课题/细胞干性课题/data/",i,"/scRNA.rds"))
  
  setwd("E:/2021020529/课题/细胞干性课题/result/各算法权重/")
  sc_allcell<-readRDS(paste(i,"/sc_allcell.rds",sep = ""))
  score<-readRDS(paste(i,"/score.rds",sep = ""))
  results<-readRDS(paste(i,"/results.rds",sep = ""))
  
  # 多种干性指数标准化
  slice_normalized <- (sc_allcell@entropies - min(sc_allcell@entropies))/diff(range(sc_allcell@entropies )) # SLICE
  mRNAsi_normalized <- (score - min(score))/diff(range(score)) # mRNAsi
  cytotrace_normalized <- (results$CytoTRACE- min(results$CytoTRACE))/diff(range(results$CytoTRACE)) # CytoTRACE
  # 加权平均值
  slice_weighted <- slice_normalized*weight[1,]
  mRNAsi_weighted <- mRNAsi_normalized*weight[2,]
  cytotrace_weighted <- cytotrace_normalized*weight[3,]
  # 总值  ##############################
  csi <- slice_weighted + mRNAsi_weighted + cytotrace_weighted
  
  stem_score <- data.frame(slice_normalized, mRNAsi_normalized, cytotrace_normalized,
                           slice_weighted, mRNAsi_weighted, cytotrace_weighted,
                           csi)
  colnames(stem_score) <- c("slice_normalized","mRNAsi_normalized","cytotrace_normalized","slice_weighted","mRNAsi_weighted","cytotrace_weighted","csi")
  write.csv(stem_score,paste0(i,"_CSI.csv"))
  
}