

library(tidyverse)

rm(list = ls())
library(limma)
load(file = "exprSet_rmdup.Rdata")
load(file = "cli01_GSE164050.Rdata")

# 机器学习交集后 c("Pdhx",  "Hmgb1",  "Slc25a3","Fdx1")
# 验证后 ("Pdhx",  "Fdx1")在Neo组下调

lassoGene <- c("Pdhx",  "Fdx1")

lassoGene <-c( "Plin2", "Atg13",  "Cdo1")
table(lassoGene %in% rownames(exprSet))
exprSet[1:4,1:4]
rt <- exprSet
### 整合临床信息
if F {
  # rt=read.table("data.exp.txt",header=T,sep="\t",row.names=1)            #读取文件
  # rt$futime=rt$futime/365   #以年为单位；若以月为单位，则除以12
  # class(rt$futime)
  # class(rt$fustat)
  # 
  
  
  rt <- exprSet
  colnames(rt)[1:2]
  clin <- cli01_GSE37745
  colnames(clin)
  clin$futime <- 0
  table(is.na(clin$`days to determined death status`))
  # clin$futime <- ifelse(clin$vital_status== "Alive",clin$last_contact_days_to,clin$death_days_to)
  clin$futime <- clin$`days to determined death status`
  clin$futime <- clin$futime/365
  
  
  clin$fustat <- 0
  clin$fustat <- ifelse(clin$dead== "no",0,1)
  str(exprSet)
  
  # exprSet02 <- apply(exprSet, 2, as.numeric)
  # exprSet02 <- as.data.frame(exprSet02)
  # exprSet02[1:4,1:4]
  # rownames(exprSet02) <- rownames(exprSet)
  # exprSet02[1:4,1:4]
  
  # exprSet <- exprSet02
  clin$sample <- clin$accession
  exprSet <- exprSet
  identical(clin$sample, colnames(exprSet))
  exprSet <- as.data.frame(t(exprSet))
  identical(clin$sample, rownames(exprSet))
  
  # max(apply(exprSet, 2, max))
  # reverse(3)
  # max(apply(, 2, max))
  # 
  # 
  # max(abs(exprSet))
  
  rt <- cbind(clin[,c("futime", "fustat")],exprSet)
  rownames(rt) <- rownames(exprSet)
  class(rt$futime)
  class(rt$fustat)
  
  
  which(is.na(rt$futime))
  which(is.na(rt$fustat))
  rt <- na.omit(rt) ### 196 个样本
  rt <- rt[rt$futime !=0,]
  class(rt$futime)
  class(rt$fustat)
  table(rt$fustat)
  rt$fustat <- as.integer(rt$fustat)
  class(rt$futime)
  class(rt$fustat)
  
  
  
  
  
}

table(lassoGene %in% colnames(rt))
which(colnames(rt) %in% lassoGene)
colnames(rt)[which(colnames(rt) %in% lassoGene)]

dput(lassoGene)

rt <- as.data.frame(t(exprSet))
rt$fustat <- cli01$sample
rt$fustat <- as.numeric(factor(rt$fustat))
rt$fustat <- rt$fustat - 1
FinalGeneExp = rt %>%
  select(lassoGene)

# FinalGeneExp = rt[,lassoGene]


# myFun = function(x){crossprod(as.numeric(x),actCoef)}
# riskScore = apply(FinalGeneExp,1,myFun)

if F {
  ############  logistic建模
  train_dat <- cbind(rt$fustat, FinalGeneExp)
  colnames(train_dat)[1] <- "group"
  model_LR <- glm(group ~ .,
                  data = train_dat,
                  family = "binomial")
  
  
  # 用之前建好的模型
  riskScore_LR <- predict(model_LR, newdata = train_dat,type="response")
  
  if F {
  ############  随机森林建模
  model_RF <- randomForest(group ~ ., 
                           data = train_dat,
                           ntree = 1000, # 树的数目，例文为1000
                           nPerm = 50, # 扰动次数，一般为50
                           mtry = floor(sqrt(ncol(train_dat)-1)), 
                           proximity = T,
                           importance = T)
  riskScore_RF <- predict(model_RF, newdata = train_dat,type="response")
  }
}

riskScore <- riskScore_LR
# riskScore <- riskScore_RF
outCol = c( "fustat", lassoGene)
risk = as.vector(ifelse(riskScore > median(riskScore), "high", "low"))
dat = cbind(rt[,outCol], riskScore=as.vector(riskScore), risk)


# km_dat <- data.frame(status = rt$fustat, time = rt$futime , group = risk)
# write.csv(km_dat , file = "./lasso/km_dat_GEO.csv" , row.names = F, col.names = T)

roc_dat <- data.frame(status = rt$fustat,  group = riskScore)
write.csv(roc_dat , file = "./02output/06-03roc_LR_dat_GSE164.csv" , row.names = F, col.names = T)

roc_dat <- data.frame(status = rt$fustat,  group = dat$Pdhx)
write.csv(roc_dat , file = "./02output/06-01roc_Pdhx_dat_GSE164.csv" , row.names = F, col.names = T)

roc_dat <- data.frame(status = rt$fustat,  group = dat$Fdx1)
write.csv(roc_dat , file = "./02output/06-02roc_Fdx1_dat_GSE164.csv" , row.names = F, col.names = T)

exprSet["Pdhx",]










