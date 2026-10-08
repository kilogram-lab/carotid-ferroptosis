#################
#################  铁死亡相关基因
library(tidyverse)
library(Seurat)
library(ggsci)
load(file = "scRNA_allcluster.Rdat")
#gene=read.table('Fxgene.txt',header = T)

scRNA_DEG <- subset(scRNA, celltype == "VSMCs")
dim(scRNA)
dim(scRNA_DEG)

scRNA <- scRNA_DEG


# rownames(scRNA_DEG) <- toupper(rownames(scRNA_DEG))

# rownames(scRNA@assays[["RNA"]]@data) <- toupper(rownames(scRNA@assays[["RNA"]]@data))


# 选取平缓的elbow，不用更改
pc.num=1:20

scRNA <- FindNeighbors(scRNA, dims = pc.num) 

# 聚类
# scRNA <- FindClusters(scRNA)
scRNA <- FindClusters(scRNA,
                      resolution = 0.5, # 最重要参数，该值越大，cluster 越多
                      method = "igraph", # 根据结果调整
                      algorithm = 1, # 根据结果调整'
                      random.seed = 2021) 
table(scRNA@meta.data$seurat_clusters)



p_celltype <- DimPlot(scRNA,
                      group.by="seurat_clusters", label=T, label.size=4.5, 
                      reduction='umap',repel = T)+
  scale_color_npg()
p_celltype

ggsave("03fig03/03_smcs_cluster.pdf", p_celltype, width=6 ,height=4)


p_celltype <- DimPlot(scRNA,
                      group.by="group", label=T, label.size=4.5, 
                      reduction='umap',repel = T)+
  scale_color_npg()
p_celltype

ggsave("03fig03/04_smcs_group.pdf", p_celltype, width=6 ,height=4)





gene <- read.csv(file = "db_ferdb.csv")
gene=gene$Symbol

# Convert the first letter to uppercase and the rest to lowercase
gene_modified <- paste0(toupper(substr(gene, 1, 1)), tolower(substr(gene, 2, nchar(gene))))


# 打印结果
print(gene_modified)
gene <- gene_modified



######################
dge.celltype <- FindMarkers(scRNA_DEG, ident.1 = 'CAI', ident.2 = 'Control', group.by = 'group')
write.csv(dge.celltype, file = "dge.celltype.csv")

length(intersect(gene, rownames(dge.celltype)))
gene <- intersect(gene, rownames(dge.celltype))



Ferroptosis_score=list(gene)
names(Ferroptosis_score)='Ferroptosis_score'

library(Seurat)
# scRNA <- AddModuleScore(object = scRNA, features = Ferroptosis_score, name ='Ferroptosis_score')

library(ggplot2)
p1=VlnPlot(scRNA,features ='Ferroptosis_score1', pt.size = 0.1,
           group.by = "group",
           )+
  scale_fill_npg()
p1
# split.by = "group",


ggsave("03fig03/05_smcs_group_fer.pdf", p1, width=6 ,height=4)



p1=VlnPlot(scRNA,features ='Ferroptosis_score1', pt.size = 0.1,
           
           group.by="seurat_clusters")+
  scale_fill_npg()
p1

ggsave("03fig03/06_smcs_cluster_fer.pdf", p1, width=6 ,height=4)



# + 
#   theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
#   theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())
# p1

library(viridis)
p2 <- FeaturePlot(scRNA,features = 'Ferroptosis_score1',reduction = 'umap',label=F,pt.size = 1,cols = c("purple","orange","darkgreen"))
FeaturePlot(scRNA,features = 'Ferroptosis_score1',reduction = 'tsne')
p2 <-FeaturePlot(scRNA,features = 'Ferroptosis_score1',reduction = 'umap',label=F,pt.size = 1,cols = pal_npg("nrc")(3)[3:1])
ggsave("03fig03/07_smcs_fer.pdf", p2, width=6 ,height=4)
# httpLOCAL_PROJECT_PATH
# ggsci

pal_npg("nrc")(3)



scRNA@meta.data$Ferroptosis_group <- ifelse(scRNA@meta.data$Ferroptosis_score1>median(scRNA@meta.data$Ferroptosis_score1),'high Ferroptosis','low Ferroptosis')
p3 = DimPlot(scRNA, group.by="Ferroptosis_group", label=F, label.size=4.5, 
             reduction='umap',cols = col_vector)+theme(title = element_blank())+ 
  scale_color_npg()
p3

ggsave("03fig03/08Ferroptosis_group.pdf", p3, width=6 ,height=4)


# ggsave("cell_identify/Ferroptosis_group_celltype.pdf", p3, width=8 ,height=4)



