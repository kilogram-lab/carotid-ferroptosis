# 官方文档
# httpLOCAL_PROJECT_PATH

# 104机器学习基于r包mlr35--分类--svm
# httpLOCAL_PROJECT_PATH

# httpLOCAL_PROJECT_PATH
#  httpLOCAL_PROJECT_PATH

# httpLOCAL_PROJECT_PATH

# httpLOCAL_PROJECT_PATH

# httpLOCAL_PROJECT_PATH
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
