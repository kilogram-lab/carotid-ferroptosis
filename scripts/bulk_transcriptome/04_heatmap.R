################################################
################################################
### 作者：果子
### 更新时间：2020-01-04
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人博客: https://codingsoeasy.com/

###############
####heatmap热图
#用行名提取数据
rm(list = ls())
## 加载表达数据
load(file = "exprSet_rmdup.Rdata")
## 加载差异列表
load(file = "diffLab.Rda")
load(file = "group.Rda")
library(dplyr)
library(tibble)
# diffLab <- allDiff %>% 
#   rownames_to_column() %>% 
#   filter(adj.P.Val < 0.01) %>% 
#   filter(abs(logFC) >4) %>% 
#   column_to_rownames()
## 加载热图的R包
if(!require("pheatmap")) install.packages("pheatmap")
library(pheatmap)
##用名称提取部分数据用作热图绘制

## 读取铁死亡基因
cup_list <- read.csv(file = "./01input/db_ferdb.csv")


genelist <- cup_list$Symbol
genelist <- gsub("(\\w)(\\w+)", "\\U\\1\\L\\2", genelist, perl = TRUE)
intersec01 <- intersect(genelist, rownames(diffLab))

heatdata <- exprSet[intersec01,]



##制作一个分组信息用于注释
#group <- c(rep("con",3),rep("treat",3)) 
annotation_col <- data.frame(group)
rownames(annotation_col) <- colnames(heatdata)
#annotation_col$group <- factor(annotation_col$group,levels = c("con","treat"),ordered = F)

## 直接作图
pheatmap(heatdata)

#如果注释出界，可以通过调整格子比例和字体修正
pheatmap(heatdata, #热图的数据
         cluster_rows = TRUE,#行聚类
         cluster_cols = TRUE,#列聚类，可以看出样本之间的区分度
         annotation_col =annotation_col, #标注样本分类
         annotation_legend=TRUE, # 显示注释
         show_rownames = F,# 显示行名
         scale = "row", #以行来标准化，这个功能很不错
         color =colorRampPalette(c("blue", "white","red"))(100),#调色
         #filename = "heatmap_F.pdf",#是否保存
         cellwidth = 40, cellheight = 2,# 格子比例
         fontsize = 10)

### 可以重新筛选，阈值设大一点
diffLab <- allDiff %>% 
  rownames_to_column() %>% 
  filter(adj.P.Val < 0.05) %>% 
  filter((logFC >1) | (logFC < -1)  ) %>% 
  column_to_rownames()

intersec01 <- intersect(genelist, rownames(diffLab))
heatdata <- exprSet[intersec01,]



##制作一个分组信息用于注释
#group <- c(rep("con",3),rep("treat",3)) 
annotation_col <- data.frame(group)
rownames(annotation_col) <- colnames(heatdata)

#如果注释出界，可以通过调整格子比例和字体修正
p_heat <- pheatmap(heatdata, #热图的数据
         cluster_rows = TRUE,#行聚类
         cluster_cols = TRUE,#列聚类，可以看出样本之间的区分度
         annotation_col =annotation_col, #标注样本分类
         annotation_legend=TRUE, # 显示注释
         show_rownames = T,# 显示行名
         scale = "row", #以行来标准化，这个功能很不错
         color =colorRampPalette(c("blue", "white","red"))(100),#调色
         #filename = "heatmap_F.pdf",#是否保存
         cellwidth = 30, cellheight = 10,# 格子比例
         fontsize = 10)

png(filename = "dif_heat_carotid.png", height = 2500, width =  3500, res = 600)
p_heat
dev.off()

## 记住保存

## GEO教程长期更新的链接是这个:
## https://codingsoeasy.com/archives/geo