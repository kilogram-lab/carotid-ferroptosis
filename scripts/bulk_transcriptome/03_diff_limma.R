################################################
################################################
### 作者：果子
### 更新时间：2020-01-04
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人博客: https://codingsoeasy.com/

##############
##############
####差异分析
#############
#加载limma包，用于校正和比较差异
rm(list = ls())
library(limma)
load(file = "exprSet_rmdup.Rdata")

## 记住保存
cli01<-data.table::fread("sample.csv",header = F)
colnames(cli01)
colnames(cli01) <- c("sample","group")
cli01$group <- c(rep("Con",4), rep("Neo",4))
cli01$sample <- c(paste0("C",1:4),paste0("N",1:4))
identical(cli01$sample , colnames(exprSet))

cli01<-as.data.frame(cli01)
a<-colnames(exprSet) ## a为样本，一个向量

#  rownames(cli)<-cli[,1]
rownames(cli01)<-cli01[,1]
colnames(cli01)

cli01_GSE164050 <- cli01
save(cli01_GSE164050, file = "cli01_GSE164050.Rdata")

group<-cli01$group
metadata<-data.frame(sample_id=cli01[a,1],
                     group=cli01[a,]$group,   
                     sample=cli01[a,]$group) ##每次要变，看看group是什么,sample为分组
#####此时为所需metadata
####
####  此处确定那两组比较，
####
# table(cli$path)
# table(cli$sample)
# cli$sample[1]

##每次要变，看看group是什么,sample为分组
if (T) {
#cli01<- cli01[,1:2]
colnames(cli01)<-c("sample_id","sample")

group<-cli01$sample
#group <- c(rep("stenosis",3),rep("normal",6))   #根据cli01$path 确定分组
group <- factor(group,levels = c("Con","Neo"),ordered = F)

}   #筛选静脉，比较正常组合狭窄组

#save(group,file = "group.Rda")

if (T) {
exprSet<-t(exprSet) #行基因->列基因； 列样本->行样本
exprSet<-as.data.frame(exprSet)
exprSet[1:4,1:4]
exprSet<-exprSet[cli01[,1],] # 行样本，根据行排序
exprSet[1:4,1:4]
exprSet<-t(exprSet)
exprSet<-as.data.frame(exprSet)
}    #根据确定分组，之前定义所需那几组数据，筛选所需矩阵
                #按照cli01排序





## 构建比较矩阵
design <- model.matrix(~group)
## 比较矩阵命名
colnames(design) <- levels(group)
#design



#differential差异分析
##1.构建分组矩阵,有两种方式
##https://dwz.cn/wR2qU7s9
#这一段有两种方法，但是初学的时候特别容易误解，这一步分完全就是独立的
## 创建分组
#group <- c(rep("con",3),rep("treat",3)) 
## 分组变成向量，并且限定leves的顺序
## levels里面，把对照组放在前面
#group <- factor(group,levels = c("con","treat"),ordered = F)
## 构建比较矩阵
design <- model.matrix(~group)
## 比较矩阵命名
colnames(design) <- levels(group)
#design

#2.线性模型拟合
fit <- lmFit(exprSet,design)
#3.贝叶斯检验
fit2 <- eBayes(fit)
#4.输出差异分析结果,其中coef的数目不能操过design的列数
# 此处的2代表的是第二列和第一列的比较
allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf) 
#save(allDiff,file = "allDiff.Rda")
write.csv(allDiff,file = "allDiff_GSE164050.csv")
#找出差异两倍以上，pvalue小于0.05，1078个
diffLab <- subset(allDiff,abs(logFC) >0.3 & P.Value < 0.05)
write.csv(diffLab,file = "./01input/diffLab_GSE164050.csv")
# diffLab <- subset(allDiff,abs(logFC) >1 & P.Value < 0.05)
# write.csv(diffLab,"DEG_sig_F1P0.5_GSE164050.csv")


diffLabUP <- subset(allDiff,logFC >0.3 & P.Value < 0.05)
diffLabDOWN <- subset(allDiff,logFC <(-0.3) & P.Value < 0.05)
#
###################################################################################
## 也可以通过dplyr来完成，更加直观
# library(dplyr)
# diffLab <- allDiff %>% 
#   filter(adj.P.Val < 0.05) %>% 
#   filter(abs(logFC) >1)
## 不过你会发现，这时候会丢失行名，这是坑，要保护一下
library(tibble)
diffLab <- allDiff %>% 
  rownames_to_column() %>% 
  filter(P.Value < 0.03) %>% 
  filter(abs(logFC) >0.3) %>% 
  column_to_rownames()
#save(allDiff,diffLab,group,diffLabUP,diffLabDOWN,file = "diffLab.Rda")

## 作图环节
##############################################
## 1.把现在数据调整成可以作图的格式

exprSet <- t(exprSet)  # 行gene ---> 列gene
## 矩阵转数据框
exprSet <- as.data.frame(exprSet)
## 增加一列为第一列
dd <- cbind(group=group,exprSet)

## 2.作图展示
library(ggplot2)
ggplot(data = dd,aes(x=group,y=Cd36,fill=group))+
  geom_boxplot()+
  geom_point()+
  theme_bw()

## 3.steal plot
my_comparisons <- list(
  c("Con", "Neo")
)
library(ggpubr)
ggboxplot(
  dd, x = "group", y = "Cd36",
  color = "group", palette = c("#00AFBB", "#E7B800"),
  add = "jitter"
)+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

## 改写成函数
diffplot <- function(gene){
  my_comparisons <- list(
    c("Con", "Neo")
  )
  library(ggpubr)
  ggboxplot(
    dd, x = "group", y = gene,
    color = "group", palette = c("#00AFBB", "#E7B800"),
    add = "jitter"
  )+
    stat_compare_means(comparisons = my_comparisons, method = "t.test")
}

diffplot("Fdx1")
diffLab["Fdx1",]
?stat_compare_means


## 4.多个基因作图查看
## 先把基因提取出来
genelist_diff <- rownames(diffLab)

## 读取铁死亡基因
cup_list <- read.csv(file = "./01input/db_ferdb.csv")


genelist <- cup_list$Symbol
genelist <- gsub("(\\w)(\\w+)", "\\U\\1\\L\\2", genelist, perl = TRUE)

intersec01 <- intersect(genelist, rownames(diffLab))

library(tools)
genelist <- gsub("(\\w)(\\w+)", "\\U\\1\\L\\2", genelist, perl = TRUE)

intersect(genelist, colnames(dd))
intersect(genelist, rownames(diffLab))
intersect(genelist, rownames(diffLabUP))
intersect(genelist, rownames(diffLabDOWN))
exprSet[,"Cbs"]

difcup <- allDiff[intersect(genelist, rownames(diffLab))[1:10],]
difcup[order(difcup$adj.P.Val, decreasing = F),]



# 输出结果
print(genelist)
intersect(genelist, colnames(dd))
genelist <- intersect(genelist, genelist_diff)[1:10]
## 再提取表达量

# genelist <- dput(gene_used)
genelist <- c( "Plin2", "Atg13",  "Cdo1")
data <- dd[,c("group",genelist)]
## 用gather调整数据
library(tidyr)
data <- data %>% 
  pivot_longer(cols=-1,
               names_to= "gene",
               values_to = "expression")
## 分面作图
library(ggsci)
ggplot(data = data,aes(x=group,y=expression,fill=group))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  facet_grid(.~gene)+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")+  scale_fill_npg()


###################################################################

saveRDS(genelist, file = "genelist.RDS")


p3 <- ggplot(data = data,aes(x=group,y=expression,fill=group))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  facet_grid(.~gene)+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

ggsave("./02output/01pic/pic02-a-box.pdf", p3, width=15 ,height=4)






# GSE126627
# genelist[!(genelist %in% rownames(exprSet02))]
# [1] "Slc31a1" "Lipt1"   "Arid1a"  "Nlrp3"  

# GSE174098
# genelist <- readRDS(file = 'genelist.RDS')
# table(genelist %in% rownames(scRNA))





## 推荐阅读
## GEO芯片分析中的大坑，差异基因完全相反！
## https://dwz.cn/TuaCFpNl
## GEO芯片如果超过了两组，也可以一次搞定差异分析
## https://dwz.cn/l6ocOQHN
## GEO芯片中配对样本如何做差异分析
## https://dwz.cn/YpniLsiP
## 因子(factor)就像贤内助，让你始终分清主次，拨开云雾。
## https://dwz.cn/KMo5SV0L

## GEO教程长期更新的链接是这个:
## https://codingsoeasy.com/archives/geo