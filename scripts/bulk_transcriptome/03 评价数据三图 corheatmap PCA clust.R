
#矩阵相关热图 ，矩阵热图
#行是基因，列是样本，纯表达矩阵
#用行名提取数据
rm(list = ls())
## 加载表达数据
load(file = "exprSet_rmdup.Rdata")
#############  



#############  热图
library(pheatmap)
pheatmap::pheatmap(cor(exprSet))
# install.packages("export")
## 加载R包
library(export)
## 导成PPT可编辑的格式
graph2ppt(file="cor_heatmap.pptx")

dev.off()
png(filename =  "cor_heatmap.png", height = 3000, width = 3000, res = 600)
pheatmap::pheatmap(cor(exprSet))
dev.off()

## 记住保存
cli01<-data.table::fread("sample.csv",header = F)
colnames(cli01)
colnames(cli01) <- c("sample","group")
cli01$group <- c(rep("neo",4),rep("con",4))
cli01$sample <- c(paste0("C",1:4),paste0("N",1:4))
identical(cli01$sample , colnames(exprSet))

# cli02 <- cli01 %>%
#   tidyr::separate(V2,into = c("group","num"),sep="patient") 
#   colnames(cli02)[1]<-"sample"
# cli01<-cli02

if(F){

if (T) {
  library(dplyr)
  # cli02<-cli02 %>%
  #  arrange(sample)
  cli02[,1]
  
  exprSet<-t(exprSet)
  exprSet<-as.data.frame(exprSet)
  exprSet[1:4,1:4]
  exprSet<-exprSet[cli02$sample,]  # 按照sample名称排序
  exprSet[1:4,1:4]  ## 行样本，列基因
  exprSet<-t(exprSet) ## 行基因，列样本
  exprSet<-as.data.frame(exprSet)
  
 
  
}     #按照cli02排序  
  ################ 
  pheatmap(cor(exprSet), #热图的数据
           cluster_rows = F,#行聚类
           cluster_cols = F,#列聚类，可以看出样本之间的区分度
           #  annotation_col =annotation_col, #标注样本分类
           annotation_legend=TRUE, # 显示注释
           show_rownames = T,# 显示行名
           show_colnames = T,# 显示行名
           scale = "row", #以行来标准化，这个功能很不错
           color =colorRampPalette(c("blue", "white","red"))(100),#调色
           #filename = "heatmap_F.pdf",#是否保存
     #      cellwidth = 25, cellheight = 10,# 格子比例
      #     fontsize = 10
     )   # 包括参数的热图
} #   按照group排序，再做热图，似乎没啥用

if (F) {
 
  ### 查看TCGA_id分组意义
  ### https://dwz.cn/WVgQUqfw
  ### 样本名称
  TCGA_id <- colnames(exprSet)
  
  
  ### 创建分组信息
  sample <- group
#  sample <- factor(sample,levels = c("recurr","cancer"),ordered = F)
  ### 获取配对信息，如果不是配对样本，就不需要这个信息
  # paire_info <- as.factor(as.numeric(as.factor(substring(TCGA_id,1,12))))
  # ### 创建metadata
  # metadata <- data.frame(TCGA_id,sample) 
  ##save(metadata,file = "metadata.Rdata")
  
  #######################################################################
  ### 核心环节，构建dds对象，前面的操作都是铺垫
  ### 要有数据countData，
  ### 要有分组信息，在colData中
  ### design部分，本次是配对
  ### 如果不是配对的是这样:design=~sample
  ### 第一列如果有基因名称，需要处理，所以tidy=TRUE
  library(DESeq2)
  dds <-DESeqDataSetFromMatrix(countData=exprSet, 
                               colData=metadata, 
                               design=~paire_info+sample,
                               tidy=F)
  ### 获取行数
  nrow(dds)
  ### 过滤
  ### 如果一个基因在所有样本中的counts数小于等于1，我们就把他删掉
  dds <- dds[rowSums(counts(dds))>1,]
  ### 获取行数
  nrow(dds)
  
  ### 数据标准化用于看聚类，下面的链接中解释了为什么使用vst方法
  ### https://dwz.cn/xJTuI4aO
  ### 第一个函数vst,用来把数据标准化，类似于取log
  vsd <- vst(dds, blind = FALSE)
  
  ### PCA，主成分分析，看分类
  plotPCA(vsd, "sample")
  
  ### PCA，主成分分析，看分类
  plotPCA(vsd, "sample")  ## exprSet->dds->vsd
  ## dds需要metadata，design=~sample。
  ###############  导出PPT
  library(export)
  ## 导成PPT可编辑的格式
  graph2ppt(file="pic_PCA_V-AV20200213.pptx")
  
  
  
} ## 此次分析没用，似乎count才可以
if (F) {
  
#  install.packages("DESeq2")

  
    library(DESeq2)
  
  library("FactoMineR")# 计算PCA
  library("factoextra")# 画图展示
  ### 第二个函数是assay
  ### assay函数提取vst标准化后的数据，保存数据用于热图
  exprSet_vst <- as.data.frame(assay(vsd))
  exprSet03<-exprSet_vst
  
  exprSet03<-t(exprSet03)  #此时exprSet03为纯表达矩阵，行为样本，列为基因
  exprSet03<-as.data.frame(exprSet03)
  Exprset.pca <- PCA(exprSet03, graph = F) #
  
  
  #######################
  group<-sample   ##此处group定义要注意，插入可能不兼容
  
  exprSet03=cbind(exprSet03,Groups=group) #cbind横向追加，即将样本分组信息追加到最后一列
  
  exprSet03[,ncol(exprSet03)]  ##此处插入可能不兼容
  
  
}

###############  PCA图
if (T) {
  #######################
  library(dplyr)
  
  #  exprSet  ###此时 exprSet行为基因，列为样本，基因名称变成行名，
  ###纯表达矩阵
  # cli01<-cli02
  cli01<-as.data.frame(cli01)
  a<-colnames(exprSet) ## a为样本，一个向量
  
  #  rownames(cli)<-cli[,1]
  rownames(cli01)<-cli01[,1]
  colnames(cli01)
  group<-cli01$group
  metadata<-data.frame(sample_id=cli01[a,1],
                       group=cli01[a,]$group,   
                       sample=cli01[a,]$group) ##每次要变，看看group是什么,sample为分组
  #####此时为所需metadata
  
  class(exprSet)   ###此时行基因，列样本
  group<-cli01$group ##此处group定义要注意，插入可能不兼容
  
group_list<-group
exprSet[1:4,1:4] ##此时 exprSet行为基因，列为样本，基因名称变成行名，
## 下面是画PCA的必须操作，需要看说明书。
exprSet=t(exprSet)#画PCA图时要求是行名时样本名，列名时探针名，因此此时需要转换
exprSet=as.data.frame(exprSet)#将matrix转换为data.frame
exprSet=cbind(exprSet,group_list) #cbind横向追加，即将分组信息追加到最后一列
library("FactoMineR")#画主成分分析图需要加载这两个包
library("factoextra") 
# The variable group_list (index = 54676) is removed
# before PCA analysis
exprSet.pca <- PCA(exprSet[,-ncol(exprSet)], graph = FALSE)#现在dat最后一列是group_list，需要重新赋值给一个dat.pca,这个矩阵是不含有分组信息的
fviz_pca_ind(exprSet.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = exprSet$group_list, # color by groups
             # palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)

png(filename = "pic_PCA_carotid.png", height = 3000, width =  3000, res = 600)
fviz_pca_ind(exprSet.pca,
             geom.ind = c("point", "text"), # show points only (nbut not "text")
             col.ind = exprSet$group_list, # color by groups
             # palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
dev.off()


#############  CLUST图
if (F) {
  sample_name<-cli01$group
  table(cli01$group)
  av<-1;
  con<-1;
  for (i in 1:length(sample_name)){
    if (sample_name[i]=="AV fistula") {
      sample_name[i]<-paste0(sample_name[i],av);
      av<-av+1
    } else {
      sample_name[i]<-paste0(sample_name[i],con);
      con<-con+1
    }
  }
  exprSet<-t(exprSet)  ### 此时行是基因，列是样本
  colnames(exprSet)<-sample_name
  
}  ##有更简单的，此处是更改行名，好确定每个标本是哪组哪个

if (T) {

  exprSet[1:4,1:4]  # 行是基因，列是样本

hc<-hclust(dist(t(exprSet)))  #hclust 输入需要 行为样本，列为基因，此处t转置完才是
plot(hclust(dist(t(exprSet))))
#colnames(exprSet)<-paste(group_list,1:22,sep='')
colnames(exprSet)
plot(hclust(dist(t(exprSet))))

}  ##基本clust


#colnames(exprSet)<-paste(group,1:22,sep='')
plot(hclust(dist(t(exprSet))))

nodePar<-list(lab.cex=0.6,pch=c(NA,19),
              cex=0.7,col='blue')
hc<-hclust(dist(t(exprSet)))
par(mar=c(5,5,5,10))
plot(as.dendrogram(hc),nodePar = nodePar,horiz = T)


#### 好看的clust图

sampleDists<-dist(t(exprSet))
hc<-hclust(sampleDists,method = "ward.D2")
plot(hc,hang=-1)

library(factoextra)
res <- hcut(sampleDists, k = 2, stand = TRUE)
# Visualize
fviz_dend(res,
          # 加边框
          rect = TRUE,
          # 边框颜色
          rect_border="cluster",
          # 边框线条类型
          rect_lty=2,
          # 边框线条粗细
          lwd=1.2,
          # 边框填充
          rect_fill = T,
          # 字体大小
          cex = 1,
          # 字体颜色
          color_labels_by_k=T,
          # 平行放置
          horiz=T)



