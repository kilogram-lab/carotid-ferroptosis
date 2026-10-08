rm(list = ls())
library(Seurat)
library(patchwork)
library(tidyverse)
library(msigdbr)
#install.packages("msigdbr")
library(clusterProfiler)
library(GSVA)
library(ggsci)



#load(file = "scRNA_allcluster.Rdat")
# 02_02_DEGs_GSEA.R


table(scRNA@meta.data$group)
# scRNA_DEG <- subset(scRNA, group == "Tumor")
table(scRNA_DEG@meta.data$celltype)
# scRNA_DEG <- subset(scRNA_DEG, celltype == "Epithelial cell")
dim(scRNA)
dim(scRNA_DEG)


dir.create("GSVA")
genesets <- msigdbr(species = "Homo sapiens", category = "H")
genesets <- subset(genesets, select = c("gs_name", "gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)

# scRNA_DEG@meta.data[["ferr_group"]]
 
#Idents 设置分组
Idents(scRNA_DEG) <- scRNA_DEG@meta.data[["ferr_group"]]
expr <- AverageExpression(scRNA_DEG, assays = "RNA", slot = "data")[[1]]
dim(expr)
# expr[1:4,1:4]
expr <- expr[rowSums(expr)>0,]
expr <- as.matrix(expr)

gsva.res <- gsva(expr,genesets,method = "ssgsea")
saveRDS(gsva.res , "gsva.res.rds")
gsva.df <- data.frame(Genesets = rownames(gsva.res) , gsva.res, check.names = F)
write.csv(gsva.df , "GSVA/gsva_res.csv" , row.names = F)

dev.off()



gsva_pheat <-pheatmap::pheatmap(gsva.res ,show_colnames = T, scale = "row"  ,color = pal_npg("nrc")(3)[1:2]) 


summary(gsva.res)
str(gsva.res)

# p2 <- FeaturePlot(scRNA,features = 'ferroptosis_score1',reduction = 'umap',label=F,pt.size = 1,cols = c("purple","orange","darkgreen"))
#################
###########  col_vector ;  cols = col_vector
# ######  p3 = DimPlot(scRNA, group.by="ferr_group", label=F, label.size=4.5, 
# reduction='umap',cols = col_vector)+theme(title = element_blank())
###########################


ggsave("GSVA/gsva_pheat.pdf", plot = gsva_pheat, width = 8, height = 10) 




###########  单细胞水平GSVA 多个基因集  ############
genesets <- msigdbr(species = "Homo sapiens", category = "C2")
genesets <- subset(genesets, gs_subcat == "CP:KEGG",
                   select = c("gs_name", "gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)


expr <- as.matrix(scRNA_DEG@assays$RNA@data)
dim(expr)
expr[1:4,1:4]
expr <- expr[rowSums(expr)>0,]


gsva.res <- gsva(expr , genesets , method == "ssgsea" , parallel.sz = 8)
# saveRDS(gsva.res,"gsva.res.KEGG.rds")
# readRDS(file = "gsva.res.KEGG.rds")
# readRDS("gsva.res.KEGG.rds")

gsva.df <- data.frame(Genesets = rownames(gsva.res) , gsva.res , checknames = F)
write.csv(gsva.df , "gsva_res_KEGG.csv" , row.names = F)

#CreateAssayObject()
scRNA_DEG[["gsva"]] <- CreateAssayObject(gsva.res)
DefaultAssay(scRNA_DEG) <- "gsva"
## 基因集名字_ 改成 -
VlnPlot(scRNA_DEG , features = "KEGG-ABC-TRANSPORTERS" , group.by = "ferr_group" , slot = "counts" , pt.size=0)
pheatmap::pheatmap(gsva.res ,show_colnames = T, scale = "row")

gsva_pheat_KEGG <- pheatmap::pheatmap(gsva.res ,show_colnames = T, scale = "row")
ggsave("GSVA/gsva_pheat_KEGG.pdf", plot = gsva_pheat_KEGG, width = 8, height = 10) 
