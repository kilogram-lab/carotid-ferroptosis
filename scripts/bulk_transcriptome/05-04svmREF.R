# 官方文档
# https://github.com/johncolby/SVM-RFE/blob/master/README.md

# 104机器学习基于r包mlr35--分类--svm
# https://lishensuo.github.io/posts/bioinfo/104%E6%9C%BA%E5%99%A8%E5%AD%A6%E4%B9%A0%E5%9F%BA%E4%BA%8Er%E5%8C%85mlr35--%E5%88%86%E7%B1%BB--svm/

# https://zhuanlan.zhihu.com/p/666557017
#  https://zhuanlan.zhihu.com/p/623191566

# https://blog.csdn.net/weixin_43216017/article/details/87898559

# https://www.jianshu.com/p/e3d3483edaa0

# https://www.jiqizhixin.com/articles/2018-10-17-20
# #安装需要的R包
# install.packages("tidyverse")
# install.packages("randomForest")
# install.packages('e1071')
# install.packages("glmnet")
# install.packages("VennDiagram")
# install.packages("ggplot2")
#加载需要的R包
rm(list = ls())

library(tidyverse)
library(glmnet)
source('msvmRFE.R')  
library(VennDiagram)
library(e1071)
library(caret)

library(mlbench)
library(caret)
## 加载表达数据
load(file = "exprSet_rmdup.Rdata")
## 加载差异列表
load(file = "diffLab.Rda")
load(file = "group.Rda")



heatdata <- exprSet[intersec01,]
data <- t(heatdata)

data <- as.data.frame(data)
str(data)

data <- cbind(group = factor(cli01$group),data )
input <- data
# 构建支持向量机模型
wine_svm <- svm(group ~ ., 
                data = data,
                type = 'C',
                kernel = 'radial')





result <- svmRFE(input, k=3, halve.above=100)
print(result)
# ??svmREF
dd <- data.frame(gene = intersec01, ranking = result)
# 按照 ranking 列升序排列
dd_sorted <- dd[order(dd$ranking), ]

# 打印排序后的数据框
print(dd_sorted)
dd_sorted$gene[1:5]
