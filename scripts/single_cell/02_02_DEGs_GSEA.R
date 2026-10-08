library(Seurat)
library(tidyverse)
load(file = "scRNA_allcluster.Rdat")
# scRNA_DEG <- subset(scRNA, group == "Tumor")
dim(scRNA)

load(file = "scRNA_allcluster.Rdat")
table(scRNA@meta.data[["group"]])
table(scRNA@meta.data[["celltype"]])

# scRNA_DEG <- subset(scRNA, group == "Tumor")
scRNA_DEG <- subset(scRNA, celltype == "VSMCs")
dim(scRNA)
dim(scRNA_DEG)



# dge.celltype <- FindMarkers(scRNA_DEG, ident.1 = 'high ferr', ident.2 = 'low ferr', group.by = 'ferr_group')

dge.celltype <- FindMarkers(scRNA_DEG, ident.1 = 'CAI', ident.2 = 'Control', group.by = 'group')
write.csv(dge.celltype, file = "dge.VSMCs.csv")
write.csv(rownames(scRNA_DEG), file = "scrna_allgene.csv")
# rownames(scRNA_DEG)


############ GSEA data

colnames(dge.celltype)
gsea_dat <- dge.celltype %>%
  rownames_to_column(var = "gene_name") %>%
  select(c("gene_name" ,"avg_log2FC" )) %>%
  arrange(desc(avg_log2FC))


colnames(gsea_dat) <- c("id","value")
write.csv(gsea_dat, file = "gsea_dat.csv", row.names = F,col.names = T)

#####################



colnames(dge.celltype)
sig_dge.celltype <- subset(dge.celltype, p_val_adj<0.05&abs(avg_log2FC)>0.2)
write.csv(sig_dge.celltype, file = "sig_dge.celltype.csv")
sig_dge.celltype_SMCS <- sig_dge.celltype

sig_dge.celltype_SMCS_UP <- subset(dge.celltype, p_val_adj<0.05& avg_log2FC >0.2)
sig_dge.celltype_SMCS_DOWN <- subset(dge.celltype, p_val_adj<0.05& avg_log2FC < (-0.2))





cyan_db <- read.table(file = "cyan.txt", sep = '\t', header = F)
darkgreen_db <- read.table(file = "darkgreen.txt", sep = '\t', header = F)
magenta_db <- read.table(file = "magenta.txt", sep = '\t', header = F)
orange_db <- read.table(file = "orange.txt", sep = '\t', header = F)
wg_db <- rbind(cyan_db,darkgreen_db,magenta_db,orange_db)
write.csv(wg_db, file = "wg_db.csv")