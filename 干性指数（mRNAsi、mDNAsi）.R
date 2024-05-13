## https://zhuanlan.zhihu.com/p/398410279

install.packages("synapser", repos = c("http://ran.synapse.org", "http://cran.fhcrc.org"))
install.packages("synapser")

##  在 https://www.synapse.org/ 网站上注册账号，并记住邮箱号密码
library(synapser)
setwd("E:/2021020529/R代码/干性指数（mRNAsi、mDNAsi）/")
## 保证使用的电脑在https://www.synapse.org/的登陆状态
synLogin(email = "heilingyue@gmail.com", password = "cty19980922")

## 登录成功之后，使用 synGet 函数获取 RNA 表达数据
synRNA <- synGet( "syn2701943", downloadLocation = "E:/2021020529/R代码/干性指数（mRNAsi、mDNAsi）/Downloads/PCBC" )

library(tidyverse)
exp <- read_delim(file = synRNA$path) %>%
  # 去除 Ensembl ID 的后缀
  separate(col = "tracking_id", sep = "\\.", into = c("Ensembl_ID", "suffix")) %>%
  dplyr::select(-suffix) %>%
  column_to_rownames("Ensembl_ID") %>%
  as.matrix()

exp[1:3,1:3]
## 将 ENSEMBL 转换为 Symbol
library(org.Hs.eg.db)

unimap <- mapIds(
  org.Hs.eg.db, keys = rownames(exp), keytype = "ENSEMBL", 
  column = "SYMBOL", multiVals = "filter"
)

data.exp <- exp[names(unimap),]
rownames(data.exp) <- unimap

synMeta <- synTableQuery("SELECT UID, Diffname_short FROM syn3156503")

metaInfo <- synMeta$asDataFrame() %>%
  dplyr::select(UID, Diffname_short) %>%
  column_to_rownames("UID") %>%
  filter(!is.na(Diffname_short))

X <- data.exp  
y <- metaInfo[colnames(X), ]
names(y) <- colnames(X)

head(y)


##  构建模型 
# 对数据进行均值中心化
X <- data.exp
m <- apply(X, 1, mean)
X <- X - m
# 将样本分为干细胞组和非干细胞组
sc <- which(y == "SC")
X.sc <- X[, sc]
X.or <- X[, -sc]
#install.packages("gelnet")
library(gelnet)
model.RNA <- gelnet(t(X.sc), NULL, 0, 1)
#得到每个基因的权重
head(model.RNA$w)
# 保存模型
save(X, y, model.RNA, file = "E:/2021020529/R代码/干性指数（mRNAsi、mDNAsi）/Downloads/PCBC/model.rda")

## 模型预测

exp <- read.table("E:/2021020529/R代码/干性指数（mRNAsi、mDNAsi）/TCGA.BRCA.sampleMap_HiSeqV2/HiSeqV2",header = T)  ## 读入TCGA 乳腺癌 RNA表达谱 
rownames(exp) <- exp$sample
exp <- exp[,-1]

## 将 ENSEMBL 转换为 Symbol
library(clusterProfiler)
library(org.Hs.eg.db)
gene_ids <- bitr(rownames(exp), fromType="ENSEMBL", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")
exp<-exp[gene_ids$ENSEMBL,]
rownames(exp)<-gene_ids$SYMBOL

# 导入模型
load("E:/2021020529/R代码/干性指数（mRNAsi、mDNAsi）/Downloads/PCBC/model.rda")
# 提取交叠基因的表达谱及权重
common <- intersect(names(model.RNA$w), rownames(exp))
X <- exp[common, ]
w <- model.RNA$w[common]

# 对于 RNA 表达数据，使用 spearman 计算权重与表达值之间的相关性来衡量样本的干性指数，并进行标准化使其落在 [0,1] 之间
score <- apply(X, 2, function(z) {cor(z, w, method="sp", use="complete.obs")})
score <- score - min(score)
score <- score / max(score)

# 样本的干性指数为
head(score)

# 封装成函数
predict.mRNAsi <- function(exp, modelPath='model.rda') {
  load(modelPath)
  
  common <- intersect(names(model.RNA$w), rownames(exp))
  X <- exp[common, ]
  w <- model.RNA$w[common]
  
  score <- apply(X, 2, function(z) {cor(z, w, method="sp", use="complete.obs")})
  score <- score - min(score)
  score <- score / max(score)
}
score <- predict.mRNAsi(exp, "E:/2021020529/R代码/干性指数（mRNAsi、mDNAsi）/Downloads/PCBC/model.rda")  ## exp 行为基因，列为样本


## mDNAsi----------------------------------------------------------
load('E:/2021020529/R代码/干性指数（mRNAsi、mDNAsi）/Downloads/PCBC/pcbc.data.Rda')
load('E:/2021020529/R代码/干性指数（mRNAsi、mDNAsi）/Downloads/PCBC/pcbc.pd.f.Rda')

m <- apply(pcbc.data, 1, mean)
pcbc.data.norm <- pcbc.data - m

SC <- pcbc.pd.f[pcbc.pd.f$Diffname_short %in% 'SC',]
X <- pcbc.data.norm[, SC$UID]
library(gelnet)
model.DNA <- gelnet(t(X), NULL, 0, 1)

save(model.DNA, model.RNA, file = "E:/2021020529/R代码/干性指数（mRNAsi、mDNAsi）/Downloads/PCBC/model-weight.rda")

length(model.DNA$w)
head(model.DNA$w)

setwd("E:/2021020529/R代码/干性指数（mRNAsi、mDNAsi）")

met<-read.table("E:/2021020529/R代码/干性指数（mRNAsi、mDNAsi）/TCGA.COAD.sampleMap_HumanMethylation27/HumanMethylation27",header = T)  ## 读入TCGA 结肠癌 甲基化谱
gene<-met$sample
met<-met[-1]
rownames(met)<-gene
## 其实这里还涉及到补缺失值的问题，补缺失值的方法有很多，为了方便，直接将缺失值替换为0
met[is.na(met)]<-0

## 对于 DNA 数据，使用线性预测模型来计算样本的干性指数

predict.mDNAsi <- function(met, modelPath='model-weight.rda') {
  load(modelPath)
  
  common <- intersect(names(model.DNA$w), rownames(met))
  X <- met[common, ]
  w <- model.DNA$w[common]
  
  score <- t(w) %*% X
  score <- score - min(score)
  score <- score / max(score)
}

score.meth <- predict.mDNAsi(as.matrix(met), "E:/2021020529/R代码/干性指数（mRNAsi、mDNAsi）/Downloads/PCBC/model-weight.rda")
colnames(score.meth)<-colnames(met)
head(t(score.meth))


