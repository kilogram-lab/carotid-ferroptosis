################################################
################################################
### 作者：果子
### 更新时间：2020-01-04
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人博客: httpLOCAL_PROJECT_PATH

## 练习GEO数据的处理流程

rm(list = ls())

#更换工作目录，可以在Rstudio中完成
# Session,Set working directory,Choose working directory,选择03_GEO

getwd()
#安装bioconductor包
#options(BioC_mirror="httpLOCAL_PROJECT_PATH")
#if(!require("limma")) BiocManager::install("limma",update = F,ask = F)

##本次处理的GEO数据编号：GSE42872
##如果网络不是很通畅，手工获取上一步matrix的链接，下载到目的文件夹
## 解压缩
#在浏览器中打开
#ftLOCAL_PROJECT_PATH
#或者先登录
#httpLOCAL_PROJECT_PATH
#在GEO accession 中输入GSE42872 也可以找到

###########
###########
#正式的练习从这里开始
###########
#先解压GSE42872_series_matrix.txt.gz，注意解压到当前目录，再读入
exprSet <- read.table("GSE164050_All.counts.txt",comment.char="!",stringsAsFactors=F,
                      header=T)
exprSet <- as.data.frame(exprSet)
exprSet <- as.numeric(exprSet)
class(exprSet)
str(exprSet)
library(tidyverse)
colnames(exprSet)
exprSet <- exprSet %>%
  column_to_rownames("AccID")
no_zero <- as.logical(apply(exprSet, 1, sum)>4)
table(no_zero)

exprSet_logical <- as.data.frame(lapply(exprSet, as.logical))
no_zero <- as.logical(apply(exprSet_logical, 1, sum)>4)

exprSet <- exprSet[no_zero,]
exprSet[1:4,1:4]


# #,fill = T, na.strings = ""
# #comment.char="!" 意思是！后面的内容不要读取，可以打开文件看一下?read.table
# apply(exprSet,1,mean)
# 
# # 因为这个GEO数据集只有一个GPL平台，所以下载到的是一个含有一个元素的list
# exprSet[is.na(exprSet)]
# ind<-which(apply(exprSet, 1, function(row) any(is.na(row)))) ##哪一行有缺失值
# 
# exprSet<-exprSet[-ind,]


# 其实不解压也可以读的
#exprSet <- read.table("GSE42872_series_matrix.txt.gz",comment.char="!",stringsAsFactors=F,header=T)

#本例暂不需要
# rownames(exprSet) <- exprSet[,1]
# exprSet <- exprSet[,-1] 
save(exprSet,cli01,cli01_GSE164050,  file = "exprSet_readGSE_GSE164050_All.counts.Rdata")
## GEO教程长期更新的链接是这个:
## httpLOCAL_PROJECT_PATH