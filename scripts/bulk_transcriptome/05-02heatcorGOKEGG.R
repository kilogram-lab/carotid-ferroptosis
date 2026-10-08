rm(list = ls())
##用ggplot2
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyverse)
load(file = "diffLab.Rda")
## 加载表达数据
load(file = "exprSet_rmdup.Rdata")


cup_list <- read.csv(file = "GeneCards-Cuproptosis-SearchResults.csv")
genelist <- cup_list$Gene.Symbol
genelist <- gsub("(\\w)(\\w+)", "\\U\\1\\L\\2", genelist, perl = TRUE)
intersec01 <- intersect(genelist, rownames(diffLab))

########   相关系数热图
corheat <- exprSet[intersec01,]
corheat <- t(corheat)
# corheat <- rownames_to_column(corheat, var = "genes")
# 
# corheat02 <- rbind(c("#group", cli01$group), c("id", cli01$sample),corheat)

# write.csv(corheat02,row.names = F,file = "./02output/02corheat.csv")

# 定义输出文件路径
output_file <- "./02output/02corheat.csv"

# 将数据框写入CSV文件，不包含行名和列名
write.table(corheat, row.names = FALSE, col.names = T, file = output_file, sep = ",")

#############    GO KEGG

# 定义输出文件路径
output_file <- "./02output/03difflabGOKEGG.csv"

# 将数据框写入CSV文件，不包含行名和列名
write.table(diffLab, row.names = T, col.names = T, file = output_file, sep = ",")







############  GO KEGG 用difflab的基因

## 读取铁死亡基因
cup_list <- read.csv(file = "./01input/db_ferdb.csv")


genelist <- cup_list$Symbol
genelist <- gsub("(\\w)(\\w+)", "\\U\\1\\L\\2", genelist, perl = TRUE)

diffLab <- allDiff %>% 
  rownames_to_column() %>% 
  filter(P.Value < 0.03) %>% 
  filter(abs(logFC) >0.3) %>% 
  column_to_rownames()

intersec01 <- intersect(genelist, rownames(diffLab))
intersec01
write.csv(intersec01, file = "./02output/03gokegg.csv")



############ GSEA data
library(tidyverse)
colnames(allDiff)
gsea_dat <- allDiff %>%
  rownames_to_column(var = "gene_name") %>%
  select(c("gene_name" ,"logFC" )) %>%
  arrange(desc(logFC))


colnames(gsea_dat) <- c("id","value")
write.csv(gsea_dat, file = "./02output/gsea_dat.csv", row.names = F,col.names = T)

