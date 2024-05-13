GSE <- c("GSE110009","GSE178318","GSE183916","GSE221575","GSE225857")
for (i in GSE) {
  print(i)
  scRNA <- readRDS(paste0("E:/2021020529/课题/细胞干性课题/data/",i,"/scRNA.rds"))
  meta <- scRNA@meta.data
  CSI <- read.csv(paste0("E:/2021020529/课题/细胞干性课题/result/各算法权重/",i,"_CSI.csv"),header = T,row.names = 1)
  
  csi_befor_ranked <- matrix(0,nrow = 1,ncol = length(unique(as.character(meta$seurat_clusters))))
  rownames(csi_befor_ranked) <- "csi"
  colnames(csi_befor_ranked) <- as.character(0:(length(unique(as.character(meta$seurat_clusters)))-1))
  for (j in unique(as.character(meta$seurat_clusters))) {
    print(j)
    cluster <- subset(scRNA, seurat_clusters == j)
    score <- CSI[which(rownames(CSI) %in% colnames(cluster)),]
    mean_csi <- mean(score$csi)
    
    csi_befor_ranked[,j] <- mean_csi
  }
  csi_ranked <- t(as.data.frame(csi_befor_ranked[,order(csi_befor_ranked,decreasing = T)]))
  rownames(csi_ranked) <- "csi"
  
  write.csv(csi_ranked,paste0("E:/2021020529/课题/细胞干性课题/data/",i,"/cluster_csi_ranked.csv"))
}