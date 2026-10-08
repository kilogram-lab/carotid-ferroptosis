#   rm(list = ls())
library(Seurat)
library(patchwork)

##  接02 PCA marker day9
load(file = "srt3_all.Rda")

###1.文件读取
###2.创建Seurat对象
###3.质控QC
###4.标准化-归一化
###5.降维，聚类
###6.细胞注释

# lung_T <- readRDS('scRNA_0724.RDS')
# #scRNA <- lung_T[,substr(colnames(lung_T),23,23) %in% 'T']
# #创建seurat
# scRNA <- CreateSeuratObject(counts = lung_T, project = "scRNA", min.cells = 3, min.features = 200)
#  load(file = "scRNA.Rdat")

library(Seurat)
library(tidyverse)
library(stringr)
library(patchwork)


scRNA <- srt3_all

dim(scRNA)

table(scRNA@meta.data[["Group"]])
scRNA@meta.data[["Group"]] <- factor(scRNA@meta.data[["Group"]], levels = c("Control" , "CAI"))


scRNA@meta.data[["group"]] <- scRNA@meta.data[["Group"]]
scRNA@meta.data[["sample"]] <- scRNA@meta.data[["Group"]]
###查看基本信息
as.data.frame(scRNA@assays$RNA@counts)
scRNA@meta.data

library(tidyverse)
# scRNA@meta.data$sample=substr(rownames(scRNA@meta.data),23,25)
# scRNA@meta.data$patient_total=substring(scRNA@meta.data$colname01,nchar(scRNA@meta.data$colname01)-2,nchar(scRNA@meta.data$colname01))
# scRNA@meta.data$colname01[1:3]
# table(substring(scRNA@meta.data$colname01,nchar(scRNA@meta.data$colname01)-2,nchar(scRNA@meta.data$colname01)))
# 
# table(substring(scRNA@meta.data$colname01,nchar(scRNA@meta.data$colname01)-2,nchar(scRNA@meta.data$colname01)-2))
# 
# sample_dat <- read.csv(file = "sample.csv", header = T)
# sample_dat02 <- sample_dat %>% 
#     separate(sample ,sep = "_", into = c("d1","d2"))
# sample_dat03 <- sample_dat02 %>% 
#   filter(d1 == "LUNG") %>%
#    mutate(group = substring(.$d2, 1,1))
# table(sample_dat03$group)



dir.create('QC')
##计算质控指标
#计算细胞中核糖体基因比例，无脑继续
Sys.setenv(LANGUAGE = "en")
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")

#计算红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")


HB.genes <- paste0(substring(HB.genes,1,1),tolower(substring(HB.genes,2)))
HB.genes
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 
Idents(scRNA) <- scRNA@meta.data[["group"]]
col.num <- length(levels(scRNA@active.ident))
library(ggplot2)
#质控前

violin <- VlnPlot(scRNA,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), ###颜色调整
                #  pt.size = 0.01,
                  group.by = "group",
                  pt.size = 0, #不需要显示点，可以设置
                  ncol = 4) #+ 
#   theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
violin

###展示多个样本的质控
VlnPlot(object = scRNA, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
        group.by  = "group",
        #cols =pal[1:2], ###颜色调整
        log = T,
        pt.size = 0)

ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 12, height = 6) 
plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot <- CombinePlots(plots = list(plot2, plot1, plot3), nrow=1, legend="none") 
ggsave("QC/pearplot_before_qc.pdf", plot = pearplot, width = 12, height = 6) 


##设置质控标准，按照文献来定，如果发现后面结果不太满意，可以调整一下阈值，但是影响不大的
print(c("请输入允许基因数和核糖体比例，示例如下：", "minGene=200", "maxGene=2500", "pctMT=10",'pctHB=3'))
minGene=200
maxGene=3000
pctMT=20
pctHB=3

##1.数据质控
dim(scRNA)
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT )
dim(scRNA)
col.num <- length(levels(scRNA@active.ident))
violin <-VlnPlot(scRNA,
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                 cols =rainbow(col.num), 
                 group.by = "group",
                 pt.size = 0, 
                 ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
violin
#  cols =rainbow(col.num)
#  cols =rainbow(30)[c(2,4)]   02
#  cols =rainbow(30)[c(3,15)]   03
#  cols =rainbow(30)[c(8,22)]   04

ggsave("QC/vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6)


###展示多个样本的质控
VlnPlot(object = scRNA, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
        group.by  = "group",
        log = T,
        pt.size = 0)

# ggsave("QC/vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6) 



# 标准化
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)###SCT标准化


#2.降维和聚类,无脑运行#########################
library(Seurat)
library(tidyverse)
#install.packages('patchwork')
library(patchwork)
dir.create("cluster")
#高变基因3000个
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 5000) 
top10 <- head(VariableFeatures(scRNA), 10) 
plot1 <- VariableFeaturePlot(scRNA) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5) +NoLegend()
plot <- CombinePlots(plots = list(plot1, plot2),legend="bottom") 
ggsave("cluster/VariableFeatures.pdf", plot = plot2, width = 6, height = 6) 

#如果内存足够最好对所有基因进行中心化
scale.genes <-  rownames(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)
##如果内存不够，可以只对高变基因进行标准化!!!!!
# scale.genes <-  VariableFeatures(scRNA)
# scRNA <- ScaleData(scRNA, features = scale.genes)

#### PCA，等一会会儿
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))###仅用高变基因来做PCA 
plot1 <- DimPlot(scRNA, reduction = "pca", group.by="sample", cols =rainbow(30)[c(1,4)] ) 
plot1

#  cols =rainbow(col.num)  cols =rainbow(30)[c(1,4)]
#  cols =rainbow(30)[c(2,25)]   02
#  cols =rainbow(30)[c(3,17)]   03
#  cols =rainbow(30)[c(8,22)]   04
ggsave("cluster/pca.pdf", plot = plot1, width = 6, height = 6)


plot2 <- ElbowPlot(scRNA, ndims=20, reduction="pca") 
plotc <- plot1+plot2
plotc
ggsave("cluster/pca3.pdf", plot = plot1, width = 6, height = 6) 

# 选取平缓的elbow，不用更改
pc.num=1:15

scRNA <- FindNeighbors(scRNA, dims = pc.num) 

# 聚类
# scRNA <- FindClusters(scRNA)
scRNA <- FindClusters(scRNA,
                         resolution = 0.2, # 最重要参数，该值越大，cluster 越多
                         method = "igraph", # 根据结果调整
                         algorithm = 1, # 根据结果调整'
                         random.seed = 2021) 
table(scRNA@meta.data$seurat_clusters)



metadata <- scRNA@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'cluster/cell_cluster.csv',row.names = F)


######################  每个cluster的marker的热图


###5.多核运算
library(future)
options(future.globals.maxSize = 20 * 1024^3)
plan(multisession,workers = 6) #开启多核运算

plan('sequential')#终止多核运算
# —————————————————————————————————————————————————————————————————————
Idents(scRNA) <- scRNA@meta.data[["seurat_clusters"]]
table(scRNA@meta.data[["seurat_clusters"]])
###6.查看所有亚群的差异基因并输出
scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers <- scRNA.markers %>% 
  select(gene , everything()) %>%
  subset(p_val < 0.05)
write.csv(all.markers , "cluster/all_cluster_marker.csv", row.names = F)
mmarkers <- scRNA.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
features <- mmarkers$gene
write.csv(mmarkers ,file="cluster/all_markers_top10.csv",row.names=F)
# —————————————————————————————————————————————————————————————————————

##top10基因绘制热图
top10_genes <- read.csv("cluster/all_markers_top10.csv")
top10_genes = CaseMatch(search = as.vector(top10_genes$gene), match = rownames(scRNA)) 
plot1 = DoHeatmap(scRNA, features = top10_genes, group.by = "seurat_clusters", group.bar = T, size = 4)
plot2 = DoHeatmap(scRNA, features = top10_genes, group.by = "seurat_clusters", group.bar = T, size = 8)

ggsave("cluster/top10_markers03.pdf", plot=plot1, width=25, height=20)
ggsave("cluster/top10_markers04.pdf", plot=plot2, width=25, height=12)

ggsave("cluster/top10_markers.png", plot=plot1, width=20, height=20)
# ————————————————————














#UMAP可视化
scRNA <- RunUMAP(scRNA, dims = pc.num)
scRNA <- RunTSNE(scRNA, dims = pc.num)
embed_umap <- Embeddings(scRNA, 'umap')
write.csv(embed_umap,'cluster/embed_umap.csv') 
plot2 = DimPlot(scRNA, reduction = "umap",label = T) 
plot2

p1 <- PCAPlot(object = scRNA,
              group.by = "seurat_clusters",
              pt.size = 2, 
              label = TRUE)

p2 <- TSNEPlot(object = scRNA,
               group.by = "seurat_clusters",
               split.by = "group",
               pt.size = 2, 
               label = TRUE)    

p3 <- UMAPPlot(object = scRNA,
               group.by = "seurat_clusters",
               split.by = "group",
               pt.size = 2, 
               label = TRUE)
p1 + p2 + p3 +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")


ggsave("cluster/UMAP.pdf", plot = p3, width = 8, height = 7)
ggsave("cluster/TSNE.pdf", plot = p2, width = 8, height = 7)

#合并tSNE与UMAP
library(patchwork)
##保存数据
Idents(scRNA)=scRNA$seurat_clusters


dir.create('cell_identify')
###3.细胞类型鉴定

# save(scRNA, file = "scRNA_allcluster.Rdat")
# load(file = "scRNA_allcluster.Rdat")


# 
habermann_imm <- c('CD274',"CD3E", "CD4", "FOXP3", "IL7R", "IL2RA", "CD40LG", "CD8A", "CCL5", "NCR1", "KLRB1", "NKG7", "LYZ", "CD68", "ITGAX", "MARCO", "FCGR1A", "FCGR3A", "C1QA", "APOC1", "S100A12", "FCN1", "S100A9", "CD14", "FCER1A", "CD1C", "CD16", "CLEC9A", "LILRA4", "CLEC4C", "JCHAIN", "IGHG1", "IGLL5", "MS4A1", "CD19", "CD79A", "CPA3", "KIT", "MKI67", "CDK1", "EPCAM")

habermann_oth <- c("VWF", "PECAM1", "CCL21", "PROX1", "ACTA2", "MYH11", "PDGFRB", "WT1", "UPK3B", "LUM", "PDGFRA", "MYLK", "HAS1", "PLIN2", "FAP", "PTPRC", "EPCAM")

habermann_imm <- paste0(substring(habermann_imm,1,1),tolower(substring(habermann_imm,2)))

habermann_oth <- paste0(substring(habermann_oth,1,1),tolower(substring(habermann_oth,2)))

DimPlot(scRNA,label = T)
DotPlot(scRNA, features = habermann_oth,group.by = "seurat_clusters") + coord_flip()
DotPlot(scRNA, features = habermann_imm,group.by = "seurat_clusters") + coord_flip()

### 载入注释
library("readxl")
## 
# new_ids_imm <- c('Conventional T cells','T cells','B cells','Regulatory T cells',
#                  'Macrophages','Mast cells','Epithelial','Plasma cells','NK cells','Myofibroblasts',
#                  'Macrophages','T cells','Epithelial','Epithelial','Regulatory T cells',
#                  'Monocytes','Monocytes','Epithelial','Stromal cells','T cells',
#                  'Regulatory T cells','Plasmacytoid dendritic cells','Regulatory T cells','Epithelial')
#0 1 12 T ;2 3 15 Macrophages; 4 11 13 14 Epithelial;
#5 B cells; 6 Fibroblasts; 7 Mast cells; 8 Monocytes
#9 Endothelial cells; 10 Plasma cells;16 Plasmacytoid dendritic cells


dim(scRNA)
cluster_celltype <- c("0"="VSMCs",
                      "1"="Fibroblasts",
                      "2"="VSMCs", 
                      "3"= "VSMCs", 
                      "4"= "ECs", 
                      "5"= "Macrophages",
                      "6"= "VSMCs", 
                      "7"= "VSMCs", 
                      "8"= "Other Cells"
                      # ,
                      # "9"="Endothelial cells"
                      
                      # ,
                      # "10"="Plasma cells",
                      # "11"="Epithelial",
                      # "12"="T cells",
                      # "13"="Epithelial",
                      # "14"="Epithelial",
                      # "15"="Macrophages",
                      # "16"="Dendritic cells"
                      )
                      
                      # "17"="Epithelial",
                      # "18"="Stromal cells",
                      # "19"="T cells",
                      # "20"="Regulatory T cells",
                      # "21"="Plasmacytoid dendritic cells",
                      # "22"="Regulatory T cells",
                      # "23"="Epithelial")
# #这里对seurat的写入方式也进行了更换，可见我们是将结果直接创建并写入了对象meta.data下的cell_type
# 并没用上一方法中的RenameIdents函数
scRNA[['celltype']] = unname(cluster_celltype[scRNA@meta.data$seurat_clusters])
celltype <- DimPlot(scRNA, reduction = 'umap', group.by = 'celltype',
split.by = "group",
        label = TRUE, pt.size = 0.5) + NoLegend()
ggsave("cluster/cell_type_normal_UMAP.pdf", plot = celltype, width = 10, height = 7)


#### 细胞构成比，转03-02 compare precent
barplot(table(scRNA[['celltype']]))
barplot(unname(table(scRNA[['celltype']])))





#Dimplot中是可以通过group.by参数指定可视化内容的。为了演示，进行了参数的补全。

# Idents(scRNA) <- scRNA@meta.data$seurat_clusters
# names(new_ids_imm) <- levels(scRNA)
# scRNA<- RenameIdents(scRNA, new_ids_imm)
# 
# scRNA@meta.data$celltype <- Idents(scRNA)

###调整配色
library(RColorBrewer) 
library(viridis)
library(wesanderson)
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
pie(rep(6,n), col=sample(color, n))
col_vector
col_vector =c(wes_palette("Darjeeling1"), wes_palette("GrandBudapest1"), wes_palette("Cavalcanti1"), wes_palette("GrandBudapest2"), wes_palette("FantasticFox1"))
pal <- wes_palette("Zissou1", 10, type = "continuous")
pal2 <- wes_palette("Zissou1", 5, type = "continuous")
pal[1:10]


p2 = DimPlot(scRNA, group.by="celltype", label=T, label.size=4.5, 
             reduction='umap',cols = col_vector,
             split.by = "group")+NoLegend()
p2
ggsave("cluster/cell_type_UMAP.pdf", plot = p2, width = 10, height = 7)
###基因的细胞定位，输入感兴趣的基因，也可以是marker基因
VlnPlot(scRNA, features = c("MS4A1", "MS4A6A",'SLC1A5'),group.by  = "celltype",
        split.by = "group") #小提琴图

FeaturePlot(scRNA, features = c("MS4A1", "MS4A6A"))#坐标映射图

RidgePlot(scRNA, features = c("MS4A1", "MS4A6A"), ncol = 1)#峰峦图

features= c('IL7R', 'CCR7','CD14', 'LYZ',  'IL7R', 'S100A4',"MS4A1", "CD8A",'FCGR3A', 'MS4A7', 'GNLY', 'NKG7','FCER1A', 'CST3','PPBP')
DoHeatmap(subset(scRNA, downsample = 100), features = features, size = 3)#热图
DoHeatmap(subset(scRNA), features = features, size = 3)#热图
###保存图片
ggsave("cell_identify/UMAP_celltype.pdf", p2, width=6 ,height=4)
Idents(scRNA)=scRNA$celltype 


dir.create("03fig03")

p_celltype <- DimPlot(scRNA,
        group.by="celltype", label=T, label.size=4.5, 
        reduction='umap',repel = T)+ 
  scale_color_npg()


ggsave("03fig03/01_celltype.pdf", p_celltype, width=6 ,height=4)


p_cellgroup <- DimPlot(scRNA,
                      group.by="group", label=T, label.size=4.5, 
                      reduction='umap',repel = T)+ 
  scale_color_npg()

# 计算每个group的细胞数量
Idents(scRNA)
Idents(scRNA)=scRNA$celltype 
Idents(scRNA)=scRNA$group
cell_counts <- table(Idents(scRNA))

# 打印每个group的细胞数量
print(cell_counts)
Idents(scRNA)=scRNA$celltype

ggsave("03fig03/02_cellgroup.pdf", p_cellgroup, width=6 ,height=4)



# 
# install.packages("patchwork")
# library(patchwork)
# ??plot_layout
# 
# # 重要参数解释
# plot_layout(
#   ncol = NULL, # 设置列数
#   nrow = NULL, # 设置行数
#   byrow = NULL, # 设置案列输出方式 行 or 列
#   widths = NULL, # 设置宽度
#   heights = NULL, # 设置高度
#   guides = NULL,
#   tag_level = NULL,
#   design = NULL
# )

#### httpLOCAL_PROJECT_PATH
library(patchwork)
p_cell_combine01 <- p_celltype + p_cellgroup +
  plot_layout(guides = "collect") 


ggsave("01fig01/01_cellcombine.pdf", p_cell_combine01, width=8 ,height=4)



?DimPlot

#saveRDS(scRNA,file ='scRNA_anno.RDS')

# save(scRNA, file = "scRNA_allcluster.Rdat")
# save(scRNA, file = "scRNA_celltype.Rdat")
# load(file = "scRNA_allcluster.Rdat")
## 每次可以用这个
# setwd('./scRNA/')
# scRNA=readRDS('scRNA_anno.RDS')




#################
#################  铁死亡相关基因



load(file = "scRNA_allcluster.Rdat")
#gene=read.table('Fxgene.txt',header = T)
gene <- read.csv(file = "db_ferdb.csv")
gene=gene$Symbol

gene <- paste0(substring(gene,1,1),tolower(substring(gene,2)) )
gene <- gene[gene %in% rownames(scRNA)]
table(gene %in% rownames(scRNA))

length(unique(gene))
gene <- unique(gene)


Ferroptosis_score=list(gene)
names(Ferroptosis_score)='Ferroptosis_score'

library(Seurat)
scRNA <- AddModuleScore(object = scRNA, features = Ferroptosis_score, name ='Ferroptosis_score',replace = T)




library(ggplot2)
p1=VlnPlot(scRNA,features ='Ferroptosis_score1', pt.size = 0,
           split.by = "group") +  theme_classic() +
    theme(axis.text.x = element_text(angle=45,
                                   hjust=1,size =10)) + 
  labs(title = "", y = " Ferroptosis_score", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")
p1

ggsave("cell_identify/ferr_celltype.pdf", p1, width=6 ,height=4)

# + 
#   theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
#   theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())
# p1

library(viridis)
p2 <- FeaturePlot(scRNA,features = 'Ferroptosis_score1',reduction = 'umap',label=F,pt.size = 1,cols = c("purple","orange","darkgreen"))
FeaturePlot(scRNA,features = 'Ferroptosis_score1',reduction = 'umap')


?FeaturePlot
p2
ggsave("cell_identify/ferr_density_celltype.pdf", p2, width=6 ,height=4)
ggsave("03fig03/03ferr_density_celltype.pdf", p2, width=6 ,height=4)

scRNA@meta.data$ferr_group <- ifelse(scRNA@meta.data$Ferroptosis_score1>median(scRNA@meta.data$Ferroptosis_score1),'high ferr','low ferr')
p3 = DimPlot(scRNA, group.by="ferr_group", label=F, label.size=4.5, 
             reduction='umap',cols = col_vector)+theme(title = element_blank())
p3

ggsave("03fig03/04ferr_group.pdf", p3, width=6 ,height=4)


# ggsave("cell_identify/ferr_group_celltype.pdf", p3, width=8 ,height=4)



p_cell_combine02 <- p_celltype + p_cellgroup + p2 + p3 +
  plot_layout(guides = "collect") 


ggsave("03fig03/02_cellcombine.pdf", p_cell_combine02, width=12 ,height=12)















# save(scRNA, file = "scRNA_allcluster.Rdat")
# load(file = "scRNA_allcluster.Rdat")
