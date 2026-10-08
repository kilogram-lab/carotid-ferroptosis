#引用包
# rm(list = ls())
# 机器学习交集后 c("Pdhx",  "Hmgb1",  "Slc25a3","Fdx1")
# 验证后 ("Pdhx",  "Fdx1")在Neo组下调


library(limma)
library(ConsensusClusterPlus)
library(tidyverse)
## 前接 03 gene_survival   exprset矩阵
## 前接 03_0
#表达输入文件

# exprSet <- normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T, las=2)
## 应该行为基因，必要时转置，用TCGA数据
load(file = "exprSet_rmdup.Rdata")
exprSet <- exprSet[,5:8]
load(file = "exprSet_gse126neo.Rdata")  
boxplot(exprSet_gse126neo,outline=FALSE, notch=T, las=2)


exprSet <- rownames_to_column(exprSet, var = "symbol")
exprSet_gse126neo <- rownames_to_column(exprSet_gse126neo, var = "symbol")
exprSet02 <- inner_join(exprSet, exprSet_gse126neo, by = "symbol")
# exprSet <- as.data.frame(t(exprSet))
exprSet[1:4,1:4]
exprSet_gse126neo[1:4,1:4]
exprSet02[1:4,1:4]
exprSet02 <- column_to_rownames(exprSet02, var = "symbol")

exprSet <- exprSet02
dim(exprSet)
# dat03$Cuproptosis %in% rownames(exprSet)
## 读取铁死亡基因
cup_list <- read.csv(file = "./01input/db_ferdb.csv")


genelist <- cup_list$Symbol
genelist <- gsub("(\\w)(\\w+)", "\\U\\1\\L\\2", genelist, perl = TRUE)

intersec01 <- intersect(genelist, rownames(exprSet))


exprSet[1:4,1:4]

rt <- exprSet
rt=as.matrix(rt)
exp <- rt

# dat03 <- read.csv(file = "01-05Mitophagy-GeneCards-SearchResults.csv")

# 机器学习交集后 c("Pdhx",  "Hmgb1",  "Slc25a3","Fdx1")
# 验证后 ("Pdhx",  "Fdx1")在Neo组下调

dat03 <- data.frame(Cuproptosis = intersec01)

# dat03 <- read.csv(file = "./lasso/lasso_gene.csv" , header = T)
# colnames(dat03)[1] <- "Cuproptosis"
#  write.csv(dat03, col.names = T, row.names = F, file = "dat03.csv")


# rt=read.table('exp.txt', header=T, sep="\t", check.names=F)
# rt=as.matrix(rt)
# rownames(rt)=rt[,1]
# exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
###  这句不知道为啥要注释掉，因为删除了很多数据，不得不注释
# data=data[rowMeans(data)>0,]

table(dat03$Cuproptosis %in% rownames(data))

# dat03$Cuproptosis[1:3]
# substring(dat03$Cuproptosis,1,1)
# 
# tolower(substring(dat03$Cuproptosis,2))
# nogene <- paste0(substring(dat03$Cuproptosis,1,1),tolower(substring(dat03$Cuproptosis,2)))
# genelist <- row.names(data)[row.names(data) %in% nogene] 


# nogene <- c("FUNDC1","PGAM5","TOMM22","UBB","ULK1")
# genelist <- row.names(data)[!(row.names(data) %in% nogene)] 
genelist <- dat03$Cuproptosis
dim(data)
rownames(data)[1:3]
dat03$Cuproptosis %in% rownames(data)


data=data[genelist,]
dim(data)
#聚类
maxK=6
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=50,
                             pItem=0.8,
                             pFeature=1,
                             title="plot/",
                             clusterAlg="km",
                             distance="euclidean",
                             seed=123456,
                             plot="png")

#基于PAC方法计算最佳K值
Kvec = 2:maxK
x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec))
names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)}#end for i# The optimal K
optK = Kvec[which.min(PAC)]
print(optK)
###optk变量为PAC方法确定的最佳K值
#输出分型结果
clusterNum=2      #分几类，根据判断标准判断
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("geneCluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$geneCluster))
cluster$geneCluster=letter[match(cluster$geneCluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="geneCluster3.txt", sep="\t", quote=F, col.names=F)

save(clusterOut, file = "clusterOut.Rdata")




