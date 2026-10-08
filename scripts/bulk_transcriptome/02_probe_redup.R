################################################
################################################
### 作者：果子
### 更新时间：2020-01-04
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人博客: https://codingsoeasy.com/

rm(list = ls())
load(file = "exprSet_readGSE_GSE164050_All.counts.Rdata")
library(tidyverse)
colnames(exprSet)
exprSet <- column_to_rownames(exprSet,"AccID")

boxplot(exprSet,outline=FALSE, notch=T, las=2)
## 组间校正
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T, las=2)

# 自动log化
ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { 
  ex[which(ex <= 0)] <- NaN
  exprSet <- log2(ex)
  print("log2 transform finished")
  }else{
    print("log2 transform not needed")
    }

exprSet <- as.data.frame(exprSet)

#################################################################
## 探针基因名转换
##platformMap 中有常见的平台个R注释包的对于关系，这是我整理的。
## 读取，这都是我们已经讲过的
platformMap <- data.table::fread("platformMap.txt")

## 平台的名称如何知道?
index <- "GPL24688"  ##悲催，本例没有R包
#paste0(platformMap$bioc_package[grep(index,platformMap$gpl)],".db")

## 安装R包
#options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
#if(!require("hugene10sttranscriptcluster.db")) BiocManager::install("hugene10sttranscriptcluster.db",update = F,ask = F)

#获取探针
#probe2symbol_df <- toTable(get("hugene10sttranscriptclusterSYMBOL"))

##如果没有安装好"hugene10sttranscriptclusterSYMBOL"，说明你不是很听话
#没有安装就load吧，没有关系的，我们就是要练习
#save(probe2symbol_df,file = "probe2symbol_df.Rdata")
#load(file = "probe2symbol_df.Rdata")

##悲催，本例没有R包 读取sofa,没有，发现数据矩阵好像是注释好的，所以下面都不要
if F {
library(tidyverse)
## 读取GEO平台注释信息soft文件,skip根据具体情况改
GPL10332_anno <-data.table::fread("GPL24688_family.soft",skip ="#",stringsAsFactors=F, header=T, data.table = F)

probe2symbol_df<-data.frame(probe_id=GPL10332_anno$ID,symbol=GPL10332_anno$GENE_SYMBOL)
c("probe_id", "symbol" ) 
#看一下symbol有没有重复，发现只有18834个，
length(unique(probe2symbol_df$symbol))
## 而探针和基因的对应关系要更多
nrow(probe2symbol_df)
# 所以需要去重，多个探针对应一个基因
class(probe2symbol_df$probe_id)
probe2symbol_df$probe_id<-as.character(probe2symbol_df$probe_id)
probe2symbol_df$symbol<-as.character(probe2symbol_df$symbol)
###探针转换以及去重，获得最终的表达矩阵
### 拆分体会
library(dplyr)
library(tibble)
exprSet <- exprSet %>%   #确保正确后把左边改为exprSet
  ## 行名转列名,因为只有变成数据框的列,才可以用inner_join
  rownames_to_column("probe_id") %>% 
  #合并探针的信息
  inner_join(probe2symbol_df,by="probe_id") %>% 
  #去掉多余信息
  dplyr::select(-probe_id) %>%  
  #重新排列
  dplyr::select(symbol,everything()) %>%  
  #求出平均数(这边的.代表上面传入的数据)
  ## .[,-1]表示去掉出入数据的第一列，然后求行的平均值
  mutate(rowMean =rowMeans(.[,-1])) %>% 
  #把表达量的平均值按从大到小排序
  arrange(desc(rowMean)) %>% 
  # 去重，symbol留下第一个
  distinct(symbol,.keep_all = T) %>% 
  #反向选择去除rowMean这一列
  dplyr::select(-rowMean) %>% 
  ## 列名转行名
  column_to_rownames("symbol")
exprSet["Gapdh",]
probe2symbol_df %>%
  filter(symbol=="OR4K1")
}


## 表达矩阵去重复 ##
colnames(exprSet)
exprSet <- exprSet %>%   #确保正确后把左边改为exprSet
  rownames_to_column("symbol")  %>%  
  #求出平均数(这边的.代表上面传入的数据)
  ## .[,-1]表示去掉出入数据的第一列，然后求行的平均值
  mutate(rowMean =rowMeans(.[,-1])) %>% 
  #把表达量的平均值按从大到小排序
  arrange(desc(rowMean)) %>% 
  # 去重，symbol留下第一个
  distinct(symbol,.keep_all = T) %>% 
  #反向选择去除rowMean这一列
  dplyr::select(-rowMean) %>% 
  ## 列名转行名
  column_to_rownames("symbol")
exprSet <- na.omit(exprSet)

exprSet["Gapdh",]



save(exprSet, cli01,cli01_GSE164050, file = "exprSet_rmdup.Rdata")
## 补充阅读部分:
## 探针对应的信息可以从平台文件获取
## https://mp.weixin.qq.com/s/nWbMO4mULgN__nPjooRDlg
## https://mp.weixin.qq.com/s/CSHdvRK6xoNJU91tpper_w
## https://mp.weixin.qq.com/s/DlioHHXQd-W-96tXLWrQvA
## NM_，NR_开头的识别号如何转换成基因名称
## https://mp.weixin.qq.com/s/FdCcliMCYj4Yb4grzIQMaA
## 非编码序列如何转换
## https://mp.weixin.qq.com/s/X8rUnEasKy3Dk-EoUAvC2A
## 如何让基因名称在多个数据库间随意转换？
## https://mp.weixin.qq.com/s/wsiceQmNVveoggiqeDSlmQ

## GEO教程长期更新的链接是这个:
## https://codingsoeasy.com/archives/geo