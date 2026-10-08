library(mlr3)
library(tidyverse)
# install.packages("mlr3")
# # 或者
# remotes::install_github("mlr-org/mlr3")

# 加载必要的库
library(mlr3)
library(mlr3learners)
library(mlr3filters)
library(mlr3viz)
library(mlr3fselect)
# 安装 mlr3 包（如果尚未安装）
# install.packages("mlr3")

# 导入 mlr3 包
library(mlr3)

## 读取铁死亡基因
cup_list <- read.csv(file = "./01input/db_ferdb.csv")


genelist <- cup_list$Symbol
genelist <- gsub("(\\w)(\\w+)", "\\U\\1\\L\\2", genelist, perl = TRUE)

diffLab <- allDiff %>% 
  rownames_to_column() %>% 
  filter(P.Value < 0.03) %>% 
  filter(abs(logFC) >0.3) %>% 
  column_to_rownames()

intersec01 <- intersect(genelist, rownames(diffLab))


heatdata <- exprSet[intersec01,]
data <- t(heatdata)

data <- as.data.frame(data)
str(data)


cli01$group <- cli01$sample
data <- cbind(data, target = cli01$group)
data$target <- factor(data$target)



# 创建任务
# 创建任务
task <- TaskClassif$new(id = "my_task", backend = data, target = "target")



# 创建随机森林学习器
learner <- lrn("classif.ranger", importance = "impurity")

# 训练模型
learner$train(task)

# 获取特征重要性
importance <- learner$importance()
importance_RF <- learner$importance()
# 打印特征重要性
print(importance)


autoplot(task,type = "pairs")
# 特征选择
selected_features <- names(importance)[importance > median(importance)]

# 创建新任务，只包含选定的特征
new_task <- task$select(selected_features)



###
# 创建一个数据框
importance_df <- data.frame(Gene = names(importance), Value = importance)

# 按照Value的降序排列
importance_df_sorted <- importance_df[order(importance_df$Value, decreasing = TRUE), ]

# 打印结果
print(importance_df_sorted)

importance_df_sorted$Gene[1:5]
write.csv(importance_df_sorted,"./02output/svmlao/3feature_RF.csv")
########### 画棒棒糖图

###调整配色
#### httpLOCAL_PROJECT_PATH
library(RColorBrewer) 
library(viridis)
library(wesanderson)
library(ggplot2)
library(ggsci)
library("scales")
pal= pal_npg("nrc")(10)
show_col(pal)
pal[1:2]

pal=pal_nejm("default")(10)
)+ scale_color_jco()





ggplot(data = importance_df_sorted[1:23,], aes(x = reorder(Gene, -Value), y = Value, fill = Value)) +
  geom_segment(aes(xend = Gene, yend = 0), color = "gray") +
  geom_point(size = 3, shape = 21) +
  scale_fill_gradient(low = pal[1], high = pal[2]) +
  labs(x = "Gene", y = "Importance", title = "Feature Importance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ scale_color_npg()
pbang <- ggplot(data = importance_df_sorted[1:20,], aes(x = reorder(Gene, -Value), y = Value, fill = Value)) +
  geom_segment(aes(xend = Gene, yend = 0), color = "gray") +
  geom_point(size = 3, shape = 21) +
  scale_fill_gradient(low = pal[1], high = pal[2]) +
  labs(x = "Gene", y = "Importance", title = "Feature Importance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ scale_color_npg()



ggsave("./02output/01pic/pic05-a-bangRF.pdf", pbang, width=8 ,height=8)





# GSE126627
# genelist[!(genelist %in% rownames(exprSet02))]
# [1] "Slc31a1" "Lipt1"   "Arid1a"  "Nlrp3"  

# GSE174098
# genelist <- readRDS(file = 'genelist.RDS')
# table(genelist %in% rownames(scRNA))
RF_gene <- importance_df$Gene

genelist
intersec01
# 不在另一个GSE的基因
c("Slc31a1", "Lipt1" ,  "Arid1a", "Nlrp3")
other_genes <- setdiff(intersec01, c("Slc31a1", "Lipt1", "Arid1a", "Nlrp3"))
print(other_genes)


Pdhx
Hmgb1
Slc25a3