### 前接03_06_cluster

load(file = "clusterOut.Rdata")
clusterOut

dat_cluster <- data.frame(sample = rownames(clusterOut)[-1], geneCluster = clusterOut$geneCluster[-1])


clin$futime <- ifelse(clin$vital_status== "Alive",clin$last_contact_days_to,clin$death_days_to)
clin$futime <- clin$futime/365
# clin$fustat <- 0
clin$fustat <- ifelse(clin$vital_status== "Alive",0,1)

clin$sample[1:3]
clin$fustat[1:3]
clin$futime[1:3]

rt <- cbind(clin[,c("sample","futime", "fustat")])



rt <- na.omit(rt) 
rt <- rt[rt$futime !=0,] ### 991 个样本
class(rt$futime)
class(rt$fustat)
table(rt$fustat)


clusterOut_now <- data.frame(sample = rownames(clusterOut)[which(rownames(clusterOut)  %in% rt$sample)], cluster_group = clusterOut[rt$sample,])
 



identical(clusterOut_now$sample, rt$sample)




riskScore <- clusterOut_now$cluster_group
# riskScore <- riskScore_RF

risk = as.vector(riskScore)


#     dir.create("cluster")
km_dat <- data.frame(status = rt$fustat, time = rt$futime , group = risk)
write.csv(km_dat , file = "./cluster/km_culster_dat04.csv" , row.names = F, col.names = T)




# timeroc_dat <- data.frame(status = rt$fustat, time = round(rt$futime*365) , group = riskScore)
# write.csv(timeroc_dat , file = "./cluster/timeroc_cluster_dat.csv" , row.names = F, col.names = T)
# 







