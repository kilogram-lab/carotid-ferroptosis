
rm(list = ls())
library(Seurat)

library(tidyverse)

load(file = "scRNA_allcluster.Rdat")
dir.create("SCENIC")
#############################################
##############################################


# 细胞通讯分析
# 先定义细胞
table(scRNA@meta.data$celltype)
scRNACD8=subset(scRNA,celltype %in% 'VSMCs')
summary(scRNACD8$Ferroptosis_score1)
scRNACD8$ferr_group='Ferroptosis_median'

summary(scRNACD8$Ferroptosis_score1)["1st Qu."]
summary(scRNACD8$Ferroptosis_score1)["3rd Qu."]
scRNACD8$Ferroptosis_score1[1:4]


for (i in 1:nrow(scRNACD8@meta.data)) {
  if (scRNACD8$Ferroptosis_score1[i]>summary(scRNACD8$Ferroptosis_score1)["3rd Qu."]) {
    scRNACD8$ferr_group[i]='Ferroptosis_high'
  }
  if (scRNACD8$Ferroptosis_score1[i]<summary(scRNACD8$Ferroptosis_score1)["1st Qu."]) {
    scRNACD8$ferr_group[i]='Ferroptosis_low'
  }
}



scRNA$ferr_group=''
scRNAother=subset(scRNA, celltype != 'VSMCs')

scRNA_chat=merge(scRNACD8,scRNAother)

scRNA_chat$Ferroptosiscell=paste0(scRNA_chat$ferr_group," ",scRNA_chat$celltype)
scRNA_chat$Ferroptosiscell[1:3]
table(scRNA_chat$Ferroptosiscell)

#devtools::install_github('sqjin/CellChat')
library(CellChat)
# 选取第一个病人，如果内存不够，再进一步减少细胞数，例如随机抽1000个
# scRNA_chat <- subset(scRNA_chat, orig.ident=='Tumor1')

meta =scRNA_chat@meta.data # a dataframe with rownames containing cell mata data

rownames(scRNA_chat@assays$RNA@data) <- toupper(rownames(scRNA_chat@assays$RNA@data))

data_input <- as.matrix(scRNA_chat@assays$RNA@data)


#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))


cellchat <- createCellChat(object = data_input, meta = meta, group.by = "Ferroptosiscell")

CellChatDB <- CellChatDB.human
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

dplyr::glimpse(CellChatDB$interaction)
##配体-受体分析
# 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)
unique(cellchat@idents)
cellchat@idents[1:3]
table(cellchat@idents)

cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

df.net<- subsetCommunication(cellchat)
colnames(df.net)
write.csv(df.net,file ='SCENIC/cellchat.csv',quote=F)

# df.net <- read.csv(file = "SCENIC/cellchat.csv" )
unique(df.net$interaction_name)
unique(df.net$source)
unique(df.net$target)


#returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways

#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
#gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.

#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,1), xpd=TRUE)
pppp <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
p2 <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

dev.off()
pdf(file = 'SCENIC/03netVisual_circle_N.pdf', width = 5.6,height = 5)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

dev.off()
pdf(file = 'SCENIC/04netVisual_circle_Interaction weights.pdf', width = 5.6,height = 5)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

ggsave(plot = pppp,file ='SCENIC/netVisual_circle_N.pdf', device = "pdf",width = 5.6,height = 5)

ggsave(plot = p2,file ='SCENIC/netVisual_circle_Interaction weights.pdf',device= "pdf", width = 5.6,height = 5)




dev.off()
table(scRNA_chat$celltype)
table(cellchat@idents)



levels(cellchat@idents)
p <- netVisual_bubble(cellchat, sources.use = c(5:7), 
                      targets.use = c(1:3), remove.isolate = FALSE)+coord_flip()
p

dev.off()
pdf(file = 'SCENIC/02bubble_Bcells_subtype.pdf', width = 10,height = 5)
p
dev.off()






#必须把上一个图关掉
dev.off()

levels(cellchat@idents)
colnames(df.net)
table(df.net$pathway_name)
netVisual_aggregate(cellchat, signaling = 'VISFATIN',
                    targets.use = c('Ferroptosis_high VSMCs','Ferroptosis_low VSMCs'),
                    sources.use = c(' ECs',' Fibroblasts', ' Macrophages'),vertex.receiver = vertex.receiver)



p_bubble= netVisual_bubble(cellchat,
                           targets.use = c('Ferroptosis_high VSMCs','Ferroptosis_low VSMCs'),
                           sources.use = c(' ECs',' Fibroblasts', ' Macrophages'),
                           
                           remove.isolate = FALSE)+coord_flip()
p_bubble
ggsave(p_bubble,file ='SCENIC/bubble_Bcells_subtype.pdf',width = 5.6,height = 5)

#必须把上一个图关掉
dev.off()

levels(cellchat@idents)
colnames(df.net)
table(df.net$pathway_name)


netVisual_aggregate(cellchat, signaling = 'SEMA3',
                    targets.use = c('Ferroptosis_high VSMCs','Ferroptosis_low VSMCs'),
                    sources.use = c(' ECs',' Fibroblasts', ' Macrophages'),vertex.receiver = vertex.receiver)

bubble_mk <- netVisual_aggregate(cellchat, signaling = 'ANGPTL', targets.use = c('Ferroptosis_high VSMCs','Ferroptosis_low VSMCs'),
                    sources.use = c(' ECs',' Fibroblasts', ' Macrophages'), vertex.receiver = vertex.receiver,layout = 'chord')
### 有时候报错，是右下角不够大
bubble_mk
ggsave(bubble_mk,file ='SCENIC/bubble_mk.pdf',width = 8,height = 8)





