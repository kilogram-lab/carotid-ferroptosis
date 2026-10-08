###单细胞文件读取
rm(list = ls())
library(Seurat)
library(tidyverse)
temp1 <- readRDS(file = "GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
temp1[1:4,1:4]

dd <- data.frame(colname01 = colnames(temp1))
rownames(dd) <- dd$colname01
dd1 <- dd %>% separate(colname01, into = c("d1","d2","d3"), sep = "_") %>% 
  mutate(group = substring(d3,1,1)) %>%
  mutate(patient = substring(d3,2)) %>%
  rownames_to_column(var = "colname01") %>% 
  filter(d2== "LUNG")
table(dd1$group)

temp1 <- temp1[,colnames(temp1) %in% dd1$colname01]
identical(colnames(temp1), dd1$colname01)



library(tidyverse)
library(data.table)
# temp1 <- fread('GSE131907_Lung_Cancer_raw_UMI_matrix.rds') %>% column_to_rownames(var ='V1')
library(Seurat)
scRNA1 <- CreateSeuratObject(counts = temp1,min.cells = 3,  
                             min.features = 200)
scRNA <- scRNA1
dim(scRNA)
dim(dd1)
scRNA@meta.data$group <- dd1$group
scRNA@meta.data$colname01 <- dd1$colname01
scRNA@meta.data$sample <- dd1$group
scRNA@meta.data$patient <- dd1$patient
table(scRNA@meta.data$sample) 
scRNA@meta.data$sample <- ifelse(scRNA@meta.data$sample == "N", "Normal", "Tumor")
scRNA@meta.data$group <- scRNA@meta.data$sample

save(scRNA , file = "scRNA.Rdat")



