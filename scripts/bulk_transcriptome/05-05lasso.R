# 安装 glmnet 包（如果尚未安装）
# install.packages("glmnet")

# 导入 glmnet 包
library(glmnet)

# # 读取数据
# data <- read.csv("data.csv")

# 将数据划分为特征和目标变量
features <- data[, -ncol(data)]
target <- data[, ncol(data)]

# 将特征数据转换为数值型矩阵
features <- as.matrix(features)



# 运行 LASSO 回归进行特征选择
lasso_model <- cv.glmnet(features, target, alpha = 1, family = "binomial")

# 获取特征选择结果
selected_features <- coef(lasso_model, s = "lambda.min")[-1, ]
print(selected_features)
plot(lasso_model, xvar = "lambda", label = TRUE)



