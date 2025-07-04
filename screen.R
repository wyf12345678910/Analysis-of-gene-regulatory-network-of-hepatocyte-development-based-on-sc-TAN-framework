# 计算非零表达比例
nonZeroExp <- rowSums(exp > 0) / ncol(exp)

# 计算MAD
geneMAD <- apply(exp, 1, mad)

# 根据数据特性动态设置阈值
nonZeroThreshold <- quantile(nonZeroExp, probs = 0.5) # 选择非零表达比例在前50%的基因
madThreshold <- median(geneMAD) # 选择MAD高于中位数的基因

# 应用过滤条件
selectedGenes <- which(nonZeroExp > nonZeroThreshold & geneMAD > madThreshold)
expFiltered <- exp[selectedGenes, ]
# 非零表达比例分布
hist(nonZeroExp, main = "Distribution of Non-zero Expression Proportion", xlab = "Non-zero Expression Proportion")

# MAD分布
hist(geneMAD, main = "Distribution of MAD", xlab = "MAD")

library(factoextra)
expFiltered_log <- log1p(expFiltered) # 对数转换以稳定方差
pcaResults <- prcomp(t(expFiltered_log))
fviz_pca_ind(pcaResults, title = "PCA of Filtered Expression Data")

library(randomForest)

# 运行随机森林
# 注意，这里不再需要转置exp矩阵
set.seed(123) # 确保可重复性
rf <- randomForest(t(exp), as.factor(type), importance = F, ntree = 500)

# 获取特征重要性并选择最重要的基因
importanceScores <- importance(rf)
topGenes <- names(sort(importanceScores, decreasing = TRUE)[1:1000])

# 根据选择的基因过滤表达矩阵
expFiltered_rf <- exp[topGenes, ]
length(intersect(unique(tf),rownames(expFiltered_rf)))

