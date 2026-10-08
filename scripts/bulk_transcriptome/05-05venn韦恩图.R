rm(list = ls())
library(tidyverse)

cup_list <- read.csv(file = "GeneCards-Cuproptosis-SearchResults.csv")
genelist <- cup_list$Gene.Symbol
genelist <- gsub("(\\w)(\\w+)", "\\U\\1\\L\\2", genelist, perl = TRUE)



( "Plin2", "Atg13",  "Cdo1")

dput(SCrna_VSMCs_deg$X)
genelist <- c(dput(SCrna_VSMCs_deg$X),"Plin2", "Atg13",  "Cdo1")
genelist <- genelist[genelist != "Ano6"]
# 计算最大长度
max_length <- max(length(genelist), length(importance_df_sorted$Gene[1:5]))

RF_gene <- c(importance_df_sorted$Gene[1:20]  ,"Plin2", "Atg13",  "Cdo1")
  
SVM_gene <- c(top.features$FeatureName[1:20]  ,"Plin2", "Atg13",  "Cdo1")

dat04 <- data.frame(ScRNADEG = c(genelist, rep("", max_length -length(genelist) )))
# 创建数据框 dat05
dat05 <- data.frame(RF = c(RF_gene, rep("", max_length - length(RF_gene))))

# 创建数据框 dat06
dat06 <- data.frame(SVM = c(SVM_gene, rep("", max_length - length(SVM_gene))))

merged_dat <- cbind(dat04, dat05, dat06)


# 打印合并后的数据框
print(merged_dat)

write.csv(merged_dat,  row.names = F, file = "./02output/05wg_venn_dat02.csv")

# 机器学习交集后 c("Pdhx",  "Hmgb1",  "Slc25a3","Fdx1")
# 验证后 ("Pdhx",  "Fdx1")在Neo组下调