################################################
################################################
### 作者：果子
### 更新时间：2020-01-04
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人博客: httpLOCAL_PROJECT_PATH

###############
###############
##volcano火山图
###############
###############
rm(list = ls())
##用ggplot2
library(ggplot2)
library(ggrepel)
library(dplyr)
load(file = "diffLab.Rda")
## 加载表达数据
load(file = "exprSet_rmdup.Rdata")


vol_dat <- data.frame(gene = rownames(allDiff), logFC = allDiff$logFC, padj = allDiff$P.Value
)


write.csv(vol_dat, file = "./02output/02vol_dat.csv", row.names= F, col.names= T)









data <- allDiff
data$gene <- rownames(data)
## 仔细观察data数据
## 如果是你自己的数据，至少有三列
## logFC，P.Value，gene
ggplot(data=data, aes(x=logFC, y =-log10(P.Value))) +
  ## 三个部分分别画点
  geom_point(data=subset(data,abs(data$logFC) <= 1),aes(size=abs(logFC)),color="black",alpha=0.1) +
  geom_point(data=subset(data,data$P.Value<0.05 & data$logFC > 1),aes(size=abs(logFC)),color="red",alpha=0.2) +
  geom_point(data=subset(data,data$P.Value<0.05 & data$logFC < -1),aes(size=abs(logFC)),color="green",alpha=0.2) +
  ## 画线
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+
  ## 主题
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black"))+
  labs(x="log2 (fold change)",y="-log10 (q-value)")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position='none')+
  ## 标签
  geom_text_repel(data=subset(data, abs(logFC) > 3), aes(label=gene),col="black",alpha = 0.8)


##换一种风格
library(ggplot2)
library(ggrepel)
data <- allDiff
data$gene <- rownames(data)
data$significant <- as.factor(data$adj.P.Val<0.05 & abs(data$logFC) > 0.5)
data$gene <- rownames(data)
pvolvano <-  ggplot(data=data, aes(x=logFC, y =-log10(adj.P.Val),color=significant)) +
  geom_point(alpha=0.8, size=1.2,col="black")+
  geom_point(data=subset(data, logFC > 1),alpha=0.8, size=1.2,col="red")+
  geom_point(data=subset(data, logFC < -1),alpha=0.6, size=1.2,col="blue")+
  labs(x="log2 (fold change)",y="-log10 (adj.P.Val)")+
  theme(plot.title = element_text(hjust = 0.4))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(0.5,-0.5),lty=4,lwd=0.6,alpha=0.8)+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black")) +
  geom_point(data=subset(data, abs(logFC) >= 4),alpha=0.8, size=4,col="green")+
  geom_text_repel(data=subset(data, abs(logFC) > 4), 
                  aes(label=gene),col="black",alpha = 0.8)

png(filename = "volvano_carotid.png", height = 5000, width =  5000, res = 600)
pvolvano
dev.off()

ggsave("./02output/01pic/pic02-b-vol.pdf", pvolvano, width=8 ,height=8)

####自己修改，p小的text不显示
ggplot(data=data, aes(x=logFC, y =-log10(adj.P.Val),color=significant)) +
  geom_point(alpha=0.8, size=1.2,col="black")+
  geom_point(data=subset(data, logFC > 1),alpha=0.8, size=1.2,col="red")+
  geom_point(data=subset(data, logFC < -1),alpha=0.6, size=1.2,col="blue")+
  labs(x="log2 (fold change)",y="-log10 (adj.P.Val)")+
  theme(plot.title = element_text(hjust = 0.4))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(0.5,-0.5),lty=4,lwd=0.6,alpha=0.8)+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black")) +
  geom_point(data=subset(data,adj.P.Val<0.05 & abs(logFC) >= 3),alpha=0.8, size=4,col="green")+
  geom_text_repel(data=subset(data,adj.P.Val<0.05 & abs(logFC) > 3), 
                  aes(label=gene),col="black",alpha = 0.8, size=2.5)
## 加载R包
library(export)
## 导成PPT可编辑的格式
graph2ppt(file="volcano02.pptx")

## GEO教程长期更新的链接是这个:
## httpLOCAL_PROJECT_PATH