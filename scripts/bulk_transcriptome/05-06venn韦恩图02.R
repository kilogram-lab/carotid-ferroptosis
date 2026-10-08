

write.csv(diffLab,file = "./01input/diffLab_GSE164050.csv")

gse164050 <- read.csv(file = "./01input/diffLab_GSE164050.csv")
gse126627 <- read.csv(file = "./01input/diffLab_GSE126627.csv")

SVMRFE_gene <- read.csv(file = "./02output/svmlao/4feature_svm.csv")
SVMRFE_gene$FeatureName
RF_gene <- read.csv(file = "./02output/svmlao/3feature_RF.csv")
RF_gene$Gene

SCrna_VSMCs_deg <- read.csv(file = "./01input/scrna.dge.VSMCs.csv")
SCrna_VSMCs_deg$X

SCrna_allgene <- read.csv(file = "./01input/scrna_allgene.csv")


SCrna_VSMCs_deg$X

intersect(SVMRFE_gene$FeatureName[1:20],intersect(RF_gene$Gene[1:20],SCrna_allgene$x))  %in% gse126627$X

gene_used <- intersect(SVMRFE_gene$FeatureName[1:20],intersect(RF_gene$Gene[1:20],SCrna_allgene$x))[intersect(SVMRFE_gene$FeatureName[1:20],intersect(RF_gene$Gene[1:20],SCrna_allgene$x))  %in% gse126627$X]

paste0(gene_used, collapse = ",")
dput(gene_used)



genelist




intersec01 <- intersect(gse164050$X,gse126627$X)
interseco2 <- intersect(intersec01, genelist)
