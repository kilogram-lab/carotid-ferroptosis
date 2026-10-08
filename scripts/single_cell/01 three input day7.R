
# 文献介绍

# 介绍Seurat包

# GEO数据库下载单细胞数据介绍

# 构建Seurat对象及讲解

# BiocManager::install('Seurat')
# BiocManager::install('scater')
# BiocManager::install('cowplot')
rm(list = ls())
library(Seurat)
library(tidyverse)
##########单细胞数据第一种数据类型构建Seurat对象########、

scRNA_data1 <- Read10X(data.dir = "data\\scType1\\")

# 构建Seurat对象
srt1 <- CreateSeuratObject(counts = scRNA_data1, 
                            # 表达矩阵，可以是稀疏矩阵，也可以是普通矩阵
                            min.cells = 3, # 去除在小于3个细胞中表达的基因
                            min.features = 200 # 去除只有200个以下基因表达 的细胞
)

scRNA_data1[c("CD3D","TCL1A","MS4A1"),1:30]

# c查看数据类型
typeof(srt1) # "S4"
class(srt1) # "dgCMatrix"

# 计算稀疏矩阵的内存大小
sparse.size <- object.size(x = srt1)
sparse.size

# 查看稀疏矩阵数据
srt1@assays[["RNA"]][1:10, 1:4]
head(srt1@meta.data)
rownames(srt1@assays[["RNA"]])[1:10]
srt1@assays[["RNA"]][c("Raet1e","Ulbp1"), 1:4]
srt3_all <- srt1
save(srt3_all, file = "srt3_SCT_final.Rda")

