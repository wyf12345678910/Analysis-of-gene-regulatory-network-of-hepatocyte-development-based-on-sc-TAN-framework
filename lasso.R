library(BiocManager)
library(sf)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(biomaRt)
library(stringr)
library(WGCNA)
library(limma)
library(scRNAtoolVis)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(VennDiagram)
library(ggvenn)
library(ggrepel)
library(pheatmap)
library(Seurat)
library(GENIE3)
library(monocle)
library(FactoMineR)#画主成分分析图需要加载这两个包
library(factoextra) 
library(sva)
library(glmnet)
rm(list = ls())
setwd("D:/R/R工作区/cell/对比")
liver = read.table(file="GSE96981_data.human.liver.csv",header=T,sep=",")
TFlist1<-read.table(file="Homo_sapiens_TF.txt",header=T,sep="\t")
TFlist2<-read.table(file="trrust_rawdata.human.txt",header=F,sep="\t")
type <- liver$category
table(type)
# fetal_10 fetal_17     hHep   hLivEC
# 83      123       82       15 
# sorted_indices <- order(type)
# sorted_type <- type[sorted_indices]
exp_liver = liver[,c(-1:-13)]
rownames(exp_liver)=liver[,1]
dat=as.data.frame(exp_liver)
# pca的统一操作走起
dat.pca <- PCA(dat, graph = FALSE)
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = type, # color by groups
             #palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
exp2 <- removeBatchEffect(t(exp_liver), type)
dat1=as.data.frame(t(exp2))
# pca的统一操作走起
dat.pca1 <- PCA(dat1, graph = FALSE)
fviz_pca_ind(dat.pca1,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = type, # color by groups
             #palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
TFlist=intersect(unique(TFlist1[,2]),rownames(exp2))
tf=intersect(unique(TFlist1[,2]),rownames(exp2))
length(intersect(unique(TFlist),rownames(exp2)))
# type <- sorted_type
# TFlist=intersect(unique(TFlist2[,1]),rownames(exp))
# tf=intersect(unique(TFlist2[,1]),rownames(exp))
# length(intersect(unique(TFlist),rownames(exp)))
# 假设原始数据是基因×细胞矩阵，需要转置
exp_matrix <- t(as.matrix(exp_liver))  # 转置为细胞×基因

# 检查矩阵结构
dim(exp_matrix)  # 应显示 [细胞数, 基因数]
head(rownames(exp_matrix))  # 确认是细胞名

# 定义批次变量（确保是因子类型）
batch <- factor(liver[,2])

# 定义协变量
mod <- model.matrix(~ as.factor(type))

# 运行ComBat（注意输入矩阵已转置）
y <- ComBat(
  dat = exp_matrix,
  batch = batch,
  mod = mod,
  par.prior = FALSE
)
exp2=y
# 后续分析（如t-SNE）
library(Rtsne)
tsne_result <- Rtsne(t(y), dims = 2)  # 注意：此时y已是细胞×基因矩阵，无需再次转置
#######################################################
del_gen1<-c()
for(i in 1:dim(exp2)[1]){
  if((length(which(exp2[i,]==0))>=(dim(exp2)[2]-3)))
    #如果某个基因在所有细胞中的表达值为零的细胞数量（列数量为dim(exp)[2]）大于等于（细胞总数-3），则执行以下代码块。
    del_gen1<-c(del_gen1,i)
  #将符合条件的基因索引添加到del_gen1向量中
}

length(del_gen1)#4640
exp2<-exp2[-del_gen1,]
#从exp数据框中移除del_gen1向量中指定的行（基因），得到一个经过初步过滤的基因表达矩阵

del_gen<-c()
for(i in 1:dim(exp2)[1]){
  if(var(as.numeric(exp2[i,]))<0.9)
    #如果某个基因在所有细胞中的表达值的方差（使用var()函数计算）小于0.9，则执行以下代码块。
    del_gen<-c(del_gen,i)#将符合条件的基因索引添加到del_gen向量中
}
length(del_gen)#4015
exp2<-exp2[-del_gen,]
dim(exp2)#10272   333

# 假设exp是基因表达矩阵，行为基因，列为样本
# 设定阈值为0.5 CPM，在至少5%的样本中有表达
threshold <- 0.5
min_samples <- ncol(exp2) * 0.05
#
# # 筛选出至少在10%的样本中表达量大于1 CPM的基因
filtered_genes <- rowSums(exp2 >= threshold) >= min_samples
exp2 <- exp2[filtered_genes, ]
dim(exp2)
#对细胞进行基因表达的标准化处理(是为了消除不同细胞之间的技术差异和批次效应，以便更准确地比较和分析基因表达水平。)
# FPKM to TPM,FPKM和TPM均为常见的两种标准化方法
# fpkmToTpm <- function(fpkm)
# {
#   exp2(log(fpkm) - log(sum(fpkm)) + log(1e6))
# }
# # 定义了一个名为fpkmToTpm的函数，用于将FPKM（Fragments Per Kilobase of transcript per Million mapped reads）转换为TPM(Transcripts Per Million)
# tpms <- apply(exp2,2,fpkmToTpm)
# #应用fpkmToTpm函数将基因表达矩阵exp2中的每列进行转换，得到TPM值。
# tpms<-exp2
m.mad <- apply(exp2,1,mad)
exp2<- exp2[which(m.mad >max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.005)),]
dim(exp2)
length(intersect(unique(tf),rownames(exp2)))
m.mad <- apply(exp2,1,mad)
dim(exp2)
#筛选时不设置mad最小值，直接使用前75%的基因或者探针。
exp2 <- exp2[which(m.mad >  max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
dim(exp2)
###########
select_list<-c()
overlap_list<-list()#用来存储每个分组条件下的差异表达基因的交集
deg_list <- list()
for(i in 1:length(table(type))){
  group_list<-rep("type0",dim(exp2)[2])
  group_list[which(type==names(table(type))[i])]<-"type1"
  #将符合当前分组条件的细胞（根据type变量的值）标记为"type1"
  design <- model.matrix(~group_list)
  #根据group_list创建一个模型矩阵design，用于线性模型（linear model）分析
  colnames(design) <- levels(group_list)
  #将design的列名设置为group_list的水平值（levels）
  rownames(design) <- colnames(exp2)
  #将design的行名设置为exp2的列名，即细胞的ID
  
  fit<-lmFit(exp2,design)
  fit<-eBayes(fit)
  #limma包差异分析
  
  deg<-topTable(fit,coef=2,adjust='BH',number = Inf)
  #根据经验贝叶斯结果，从头部（top）选择最显著的差异表达基因（默认选择第2个系数）
  deg_list[[i]] <-deg  
  head(deg) 
  logFC_1<-rbind(deg[deg$logFC>1,],deg[deg$logFC< -1,])
  # logFC_1<-deg[deg$logFC>1,]
  #选择对数折叠变化（log-fold change）大于1的差异表达基因。
  #pvalue_0.1<-deg[deg$adj.P.Val<0.00005,]
  pvalue_0.05<-deg[deg$P.Value<0.05,]#可以分两种情况讨论不同的p值
  #选择调整后的P值小于0.05的差异表达基因。
  #select_list<-c(select_list,c(rownames(logFC_1),rownames(pvalue_0.1)))
  select_list[[i]]<-c(intersect(rownames(logFC_1),rownames(pvalue_0.05)))
  #将差异表达基因的交集添加到select_list中
  #select_list<-c(select_list,rownames(logFC_1))
  #select_list<-c(select_list,rownames(pvalue_0.1))
  datExpr<-exp2[rownames(exp2)%in%select_list[[i]],]
}
datExpr<-unique(datExpr)#去除select_list中的重复元素
length(select_list)#10332


#从exp2中选取在select_list中出现的基因行，得到一个经过筛选和选择的基因
dim(datExpr)
##########lasso
# datExpr=t(datExpr)
group_name <- names(table(type))[1]  # 可改为 "fetal_10" 等
group_cells <- which(type == group_name)

# 提取该组的表达子矩阵
group_DEGs <- rownames(datExpr)
target_genes <- intersect(group_DEGs, rownames(datExpr))
dat_subset <- datExpr[target_genes, group_cells]  # DEG表达

# TF 表达子集
TFs <- intersect(tf, rownames(datExpr))
TF_expr <- datExpr[TFs, group_cells]

# Step 2: LASSO 模型构建 TF→Gene 网络
lasso_result <- list()

for (gene in rownames(dat_subset)) {
  Y <- as.numeric(dat_subset[gene, ])
  X <- t(TF_expr)

  fit <- try(cv.glmnet(X, Y, alpha = 1), silent = TRUE)
  if (!inherits(fit, "try-error")) {
    coefs <- coef(fit, s = "lambda.1se")
    selected_TFs <- rownames(coefs)[which(abs(coefs) > 0.05 & rownames(coefs) != "(Intercept)")]
    if (length(selected_TFs) > 0) {
      lasso_result[[gene]] <- selected_TFs
    }
  }
}

# Step 3: 输出调控网络信息
all_TFs <- unlist(lasso_result)
cat("共识别出", length(lasso_result), "个被调控基因\n")
cat("共识别出", length(unique(all_TFs)), "个 TF\n")
cat("共识别出", length(all_TFs), "条 TF→gene 调控边\n")

# Top TFs
tf_freq <- sort(table(all_TFs), decreasing = TRUE)
cat("Top 10 调控最广的 TF:\n")
print(head(tf_freq, 10))

tf_data <- data.frame(TF = unique(all_TFs))

# 将结果保存为CSV文件
write.csv(tf_data, "TFs_in_regulatory_network.csv", row.names = FALSE)
#############MRNET
library(parmigene)
library(igraph)

# 1. 确保数据格式正确
datExpr <- as.matrix(datExpr)
n_genes <- ncol(datExpr)  # 获取基因数量
n_samples <- nrow(datExpr)  # 获取样本数量

# 2. 转置数据：使行为基因，列为样本
datExpr_t <- t(datExpr)  # 转置后：行=基因，列=样本

# 3. 计算互信息矩阵
mi_matrix <- knnmi.all(datExpr_t)  # 现在正确计算基因间的互信息

# 4. 使用MRNET算法构建网络
mrnet_matrix <- mrnet(mi_matrix)

# 5. 筛选显著关联
cutoff <- quantile(mrnet_matrix, 0.99, na.rm = TRUE)
adj_matrix <- ifelse(mrnet_matrix > cutoff, 1, 0)
diag(adj_matrix) <- 0  # 移除自环

# 6. 构建基因调控网络
# 现在邻接矩阵维度应该与基因数量一致
if (nrow(adj_matrix) == n_genes) {
  gene_names <- colnames(datExpr)  # 原始数据列名是基因名
  rownames(adj_matrix) <- gene_names
  colnames(adj_matrix) <- gene_names
} else {
  stop("邻接矩阵维度与基因数量不匹配。原始基因数: ", n_genes, 
       "，邻接矩阵维度: ", nrow(adj_matrix), "x", ncol(adj_matrix))
}

# 7. 创建网络对象
net <- graph_from_adjacency_matrix(
  adj_matrix,
  mode = "undirected",  # MRNET生成对称矩阵
  diag = FALSE
)

# 8. 简化网络
net <- simplify(net)

# 9. 识别关键基因（高度基因）
node_degrees <- degree(net)
key_genes <- names(sort(node_degrees, decreasing = TRUE)[1:min(50, length(node_degrees))])

# 输出结果
cat("找到", length(key_genes), "个关键基因:\n")
print(key_genes)

# 10. 可视化网络（仅显示关键基因）
V(net)$size <- 3
V(net)$color <- "lightblue"
V(net)$label <- NA  # 默认不显示标签

# 标记关键基因
V(net)[key_genes]$size <- 8
V(net)[key_genes]$color <- "tomato"
V(net)[key_genes]$label <- key_genes

# 创建布局 - 对于大型网络使用更高效的布局算法
if (vcount(net) > 500) {
  layout <- layout_with_drl(net)  # 适合大型网络
} else {
  layout <- layout_with_fr(net)   # 适合小型网络
}

# 绘图
plot(net, 
     layout = layout,
     edge.arrow.size = 0.2,
     vertex.label.cex = 0.7,
     main = "基因调控网络 (MRNET方法)")
# 获取关键基因的度中心性信息
gene_degrees <- degree(net)[key_genes]

# 创建数据框
key_genes_df <- data.frame(
  Gene_Symbol = key_genes,
  Degree = gene_degrees,
  Rank = rank(-gene_degrees)  # 按度降序排名
)

# 按度降序排列
key_genes_df <- key_genes_df[order(-key_genes_df$Degree), ]

# 添加行号作为ID
key_genes_df$ID <- 1:nrow(key_genes_df)

# 重新排列列顺序
key_genes_df <- key_genes_df[, c("ID", "Gene_Symbol", "Degree", "Rank")]

# 保存为CSV文件
write.csv(key_genes_df, "Key_Genes_List.csv", row.names = FALSE)

# 输出结果
cat("已保存", nrow(key_genes_df), "个关键基因到 Key_Genes_List.csv 文件\n")
cat("文件路径:", getwd(), "\n")

# 查看前10个关键基因
cat("\n前10个关键基因:\n")
print(head(key_genes_df, 10))

############################拟时序分析##############################
p_data <- data.frame(cell_id=colnames(datExpr),cell_type=type)
rownames(p_data)<-p_data[,1]
f_data <- data.frame(gene_short_name = rownames(datExpr),row.names = rownames(datExpr))
#使用exp数据框的行索引创建一个数据框，其中包含基因的短名（即去除了转录本ID的外部基因名）

pd <- new("AnnotatedDataFrame", data = p_data)
fd <- new("AnnotatedDataFrame", data = f_data)
#cds <- newCellDataSet(as.matrix(2.^exp-1),
cds <- newCellDataSet(matrix(unlist(datExpr),ncol=dim(datExpr)[2],nrow=dim(datExpr)[1]),
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.1,
                      expressionFamily = tobit(Lower = 0.1))
#使用newCellDataSet()函数创建一个Monocle对象cds，其中包含基因表达矩阵、样本信息、基因信息和一些其他参数

# Next, use it to estimate RNA counts
rpc_matrix <- relative2abs(cds, method = "num_genes")

# Now, make a new CellDataSet using the RNA counts
cds <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
#使用newCellDataSet()函数创建一个Monocle对象cds，其中包含基因表达矩阵、样本信息、基因信息和一些其他参数

##step2: size factor and Dispersions
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds
##step3: filter low-quality genes
cds <- detectGenes(cds, min_expr =0.1) #add a col "num_cells_expressed" in fData(cds)
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 2)) #filter out genes expressed in less than 2 cells
length(expressed_genes)
cds <- cds[expressed_genes,]
cds
print(head(pData(cds)))
##使用monocle选择的高变基因
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1&dispersion_empirical >= 1*dispersion_fit)$gene_id
#选择平均表达大于等于0.1且经验离散度大于等于1倍离散度拟合值的基因，并提取它们的基因ID。这些基因可以被视为具有高变异性的基因。
length(disp.genes)
cds <- setOrderingFilter(cds, disp.genes)
#将基因过滤为在上一步中选择的高变异性基因集，以便对其进行后续分析。
#cds <- setOrderingFilter(cds, disp.genes[1:500])
ordering_genes_plot <- plot_ordering_genes(cds)

# 使用 theme() 函数调整图表的字体大小
ordering_genes_plot <- ordering_genes_plot + 
  theme(
    plot.title = element_text(size = 20),       # 标题字体大小
    axis.title = element_text(size = 26),       # 轴标签字体大小
    axis.text = element_text(size = 18),        # 轴刻度字体大小
    legend.title = element_text(size = 16),     # 图例标题字体大小
    legend.text = element_text(size = 14)       # 图例文字字体大小
  )

# 打印图表
print(ordering_genes_plot)
# ggsave("高突基因.pdf", ordering_genes_plot, width = 6.89, height = 5.67, units = "in",dpi=300)
# plot_pc_variance_explained(cds,return_all = F)


##reduce dimension
cds <- reduceDimension(cds, reduction_method = 'DDRTree')
cds <- orderCells(cds)
plot_cell_trajectory(cds,color_by="type",size=1,show_backbone=TRUE)
plot_cell_trajectory(cds,color_by="State",show_backbone=FALSE)
cds <- orderCells(cds,root_state=4)

#参数和图像有关系
plot_cell_trajectory(cds,color_by="Pseudotime",show_backbone=FALSE)

# gene1 <- c("SON","HMGB2","NFE2","KLF1","LYL1","TFDP2","SP1","TADA2A")
# gene2 <- c("IKZF1","PKMYT1","SLC14A1","SLC22A4","ARHGEF39")
# gene3 <- c("MATN2","MEIS2","MEN1","MTF1","PACS1")
# gene4 <- c("IGFBP3","ACTA2","ARHGAP24","BGN","CCDC80","COL1A1","COL1A2")
# plot_pseudotime_heatmap(cds[gene1,],show_rownames = TRUE,return_heatmap= FALSE)
# plot_pseudotime_heatmap(cds[gene2,],show_rownames = TRUE,return_heatmap= FALSE)
# 
# plot_pseudotime_heatmap(cds[gene3,],show_rownames = TRUE,return_heatmap= FALSE)
# plot_pseudotime_heatmap(cds[gene4,],show_rownames = TRUE,return_heatmap= FALSE)
# 
# 
# my_cds_subset=cds
# # 拟时序数据和细胞位置在pData 中
# head(pData(my_cds_subset))
# 
# # 这个differentialGeneTest会比较耗费时间，测试每个基因的拟时序表达
# my_pseudotime_de <- differentialGeneTest(my_cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 4 )#cores调用的核心数
# 
# head(my_pseudotime_de)
# phe=pData(cds)
# head(phe)
# table(phe$State,phe$cell_type)
# 
# library(dplyr)
# #%>% 管道函数 把左边的值发送给右边的表达式，并作为右件表达式函数的第一个参数
# # my_pseudotime_de %>% arrange(qval) %>% head()
# # 
# # # 保存前六的基因
# # my_pseudotime_de %>% arrange(qval) %>% head() %>% select(gene_short_name) -> my_pseudotime_gene
# # my_pseudotime_gene=my_pseudotime_gene[,1]
# # my_pseudotime_gene
# 
# # #绘制一个或多个基因的拟时序
# # plot_genes_in_pseudotime(my_cds_subset[my_pseudotime_gene,])
# # ggsave('monocle_top6_pseudotime_by_state.pdf')
# # 
#%>% 管道函数 把左边的值发送给右边的表达式，并作为右件表达式函数的第一个参数
# my_pseudotime_de %>% arrange(qval) %>% tail()
# 
# # 保存前六的基因
# my_pseudotime_de %>% arrange(qval) %>% tail() %>% select(gene_short_name) -> my_pseudotime_gene
# my_pseudotime_gene=my_pseudotime_gene[,1]
# my_pseudotime_gene

#绘制一个或多个基因的拟时序
# my_pseudotime_gene=c("MT1G","HBA1","CAT","ATP5I")
# plot_genes_in_pseudotime(my_cds_subset[my_pseudotime_gene,])
# # ggsave('monocle_down6_pseudotime_by_state.pdf')
# 
# 
# 
# ##cell ordering
cells_pseudotime<-data.frame(cell_id=as.character(cds@phenoData@data$cell_id),pseudotime=cds@phenoData@data$Pseudotime)
rownames(cells_pseudotime)<-cells_pseudotime[,1]

type_pseudotime<-type[order(cells_pseudotime[,"pseudotime"]) ]
#type_pseudotime<-gsub("\\d","",type_pseudotime)

cells_pseudotime<-cells_pseudotime[order(cells_pseudotime[,"pseudotime"]),] 
cells_pseudotime<-cbind(cells_pseudotime[,1],type_pseudotime,cells_pseudotime[,2])
colnames(cells_pseudotime)<-c("cell_id","cell_type","pseudotime")
rownames(cells_pseudotime)<-cells_pseudotime[,1]
cells_pseudotime<-as.data.frame(cells_pseudotime)
write.table(cells_pseudotime, "./order_cells-pseudotime_hepatoLineage_Monocle.csv", row.names=TRUE,col.names=TRUE,sep = "," )
cells_pseudotime<-read.csv("./order_cells-pseudotime_hepatoLineage_Monocle.csv",header=T)


###select proper window size
win_size<-3
max_num<-ceiling(dim(datExpr)[2]/win_size)
colnames(datExpr)
exp_order<-datExpr[,rownames(cells_pseudotime)]
exp_new<-c()
cells_pseudotime_new<-c()
type_new<-c()
#tables<-c()
for(i in 1:max_num){
  cell_in_win<-c((win_size*(i-1)+1):min((win_size*i),dim(exp_order)[2]))
  tables1<-length(table(cells_pseudotime[cell_in_win,2]))
  j=1
  if(tables1>1){
    if((cells_pseudotime[cell_in_win[j],2]==cells_pseudotime[cell_in_win[j+1],2])&(cells_pseudotime[cell_in_win[j+1],2]!=cells_pseudotime[cell_in_win[j+2],2]))
    {
      exp_new<-cbind(exp_new,rowMeans(exp_order[,cell_in_win[c(j,j+1)]]),exp_order[,cell_in_win[j+2]])
      cells_pseudotime_new<-c(cells_pseudotime_new,mean(as.numeric(cells_pseudotime[cell_in_win[c(j,j+1)],3])),as.numeric(cells_pseudotime[cell_in_win[j+2],3]))
      type_new<-c(type_new,cells_pseudotime[cell_in_win[c(j,j+2)],2])
    }
    if((cells_pseudotime[cell_in_win[j],2]!=cells_pseudotime[cell_in_win[j+1],2])&(cells_pseudotime[cell_in_win[j+1],2]==cells_pseudotime[cell_in_win[j+2],2])){
      exp_new<-cbind(exp_new,exp_order[,cell_in_win[j]],rowMeans(exp_order[,cell_in_win[c(j+1,j+2)]]))
      cells_pseudotime_new<-c(cells_pseudotime_new,as.numeric(cells_pseudotime[cell_in_win[j],3]),mean(as.numeric(cells_pseudotime[cell_in_win[c(j+1,j+2)],3])))
      
      type_new<-c(type_new,cells_pseudotime[cell_in_win[c(j,j+1)],2])
    }
    if((cells_pseudotime[cell_in_win[j],2]!=cells_pseudotime[cell_in_win[j+1],2])&(cells_pseudotime[cell_in_win[j+1],2]!=cells_pseudotime[cell_in_win[j+2],2])){
      exp_new<-cbind(exp_new,exp_order[,cell_in_win])
      cells_pseudotime_new<-c(cells_pseudotime_new,as.numeric(cells_pseudotime[cell_in_win,3]))
      type_new<-c(type_new,cells_pseudotime[cell_in_win,2])
    }
    
  }else{
    exp_new<-cbind(exp_new,rowMeans(exp_order[,cell_in_win]))
    cells_pseudotime_new<-c(cells_pseudotime_new,mean(as.numeric(cells_pseudotime[cell_in_win,3])))
    type_new<-c(type_new,cells_pseudotime[cell_in_win[j],2])
  }
  
  #tables<-c(tables,length(table(cells_pseudotime[(win_size*(i-1)+1):(win_size*i),2])))
  # 根据需要添加更多的打印语句
}
#table(tables)
cells_pseudotime_new<-data.frame(cell_id=paste("t",c(1:(length(cells_pseudotime_new))),sep=""),pseudotime=cells_pseudotime_new)
rownames(cells_pseudotime_new)<-cells_pseudotime_new[,1]
colnames(exp_new)<-paste("t",c(1:(dim(cells_pseudotime_new))[1]),sep="")
exp_new<-as.data.frame(exp_new)

write.table(cells_pseudotime_new, "./order_cells_new-pseudotime_hepatoLineage_Monocle.csv", row.names=TRUE,col.names=TRUE,sep = "," )
cells_pseudotime_new<-read.csv("./order_cells_new-pseudotime_hepatoLineage_Monocle.csv",header=T)
#按照每组中出现频率最大的细胞类型命名的错误率
error_rate<-c()
for(win_size in 2:15){
  max_num<-ceiling(dim(datExpr)[2]/win_size)
  
  exp_new1<-c()
  cells_pseudotime_new1<-c()
  type_new1<-c()
  wrong_num<-0
  for(i in 1:max_num){
    cell_in_win<-c((win_size*(i-1)+1):min((win_size*i),dim(exp_order)[2]))
    tables1<-sort(table(cells_pseudotime[cell_in_win,2]),decreasing=T)
    if(length(tables1)>1)
      wrong_num<-wrong_num+1
    if(length(cell_in_win)>1)
      exp_new1<-cbind(exp_new1,rowMeans(exp_order[,cell_in_win]))
    if(length(cell_in_win)==1)
      exp_new1<-cbind(exp_new1,exp_order[,cell_in_win])
    cells_pseudotime_new1<-c(cells_pseudotime_new1,mean(as.numeric(cells_pseudotime[cell_in_win,3])))
    type_new1<-c(type_new1,names(tables1)[1])
    
    #tables<-c(tables,length(table(cells_pseudotime[(win_size*(i-1)+1):(win_size*i),2])))
    
  }
  error_rate<-c(error_rate,wrong_num/length(type_new1))
}
error_rate_table<-data.frame(win_size=c(2:15),error_rate=error_rate)
plot(x=c(2:15),y=error_rate,xlab="window size",ylab="error rate",type="l")

# 初始化存储错误率的向量
# error_rate <- c()
# # 循环不同的窗口大小
# for(win_size in 2:15){
#   max_num <- ceiling(dim(datExpr)[2] / win_size)
#   # 初始化存储每个窗口伪时间标准差的向量
#   std_dev_pseudotime <- numeric(max_num)
#   # 遍历每个窗口
#   for(i in 1:max_num){
#     # 确定当前窗口内的细胞
#     cell_in_win <- seq((win_size * (i - 1) + 1), min((win_size * i), dim(datExpr)[2]))
#     # 计算并存储窗口内伪时间的标准差
#     std_dev_pseudotime[i] <- sd(cells_pseudotime[cell_in_win, "pseudotime"])
#   }
#   # 计算并存储当前窗口大小的平均伪时间标准差
#   error_rate[win_size - 1] <- mean(std_dev_pseudotime, na.rm = TRUE)
# }
# # 创建一个包含窗口大小和相应错误率的数据框
# error_rate_table <- data.frame(win_size = 2:15, error_rate = error_rate)
# # 绘制窗口大小与错误率的关系图
# plot(error_rate_table$win_size, error_rate_table$error_rate, type = "b", 
#      xlab = "Window Size", ylab = "Average Pseudo-time Std Dev",
#      main = "Error Rate by Window Size")

###########################WGCNA detect module##########################################
exp_new<-t(exp_new)
powers = c(c(1:30))
# Call the network topology analysis function
sft = pickSoftThreshold(exp_new, powerVector = powers, verbose = 5,
                        networkType = "signed",RsquaredCut = 0.85)
#pickSoftThreshold说明：Rows correspond to samples and columns to genes.
#根据基因表达矩阵datExpr以及指定的幂次数组，计算拓扑分析结果。pickSoftThreshold函数将进行拓扑分析，并返回一个包含结果的对象sft。这个函数还可以通过verbose参数来控制输出的详细程度，以及通过networkType参数指定网络类型，这里是"signed"代表有符号网络。RsquaredCut参数指定拟合的R平方阈值，用于确定Soft Threshold的值。

# Plot the results:
##sizeGrWindow(9, 5)0
# png(filename = "Soft Threshold.png",width=1000,height=500)
par(mfrow = c(1,2));
cex1 = 0.85;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
#绘制拓扑分析结果中的Scale-Free拟合指标图。在此图中，x轴表示Soft Threshold的值，y轴表示Scale-Free拟合指标的负值与符号之积。每个值都用红色的文字表示。在图中标记每个点的幂次值

# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")#在图中添加一条红色水平线，表示R平方的阈值为0.90
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold(power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#绘制拓扑分析结果中的平均连接性图。在此图中，x轴表示Soft Threshold的值，y轴表示平均连接性。每个值都用红色的文字表示。在图中标记每个点的幂次值。
# dev.off()

sft$powerEstimate#mef.beta=9
#这行代码获取先前拓扑分析的结果对象sft中估计的最佳幂次值。

cor <- WGCNA::cor
#将"WGCNA"包中的cor函数赋给cor变量。
net = blockwiseModules(
  exp_new,
  power = 25,
  maxBlockSize = 2000,
  TOMType = "signed",
  networkType = "signed",
  minModuleSize = 50,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = F,
  verbose = 3,
)

#使用blockwiseModules函数进行拓扑重建模块分析。这里提供了各种参数，如幂次值（根据先前估计的最佳幂次值），最大块大小，TOM的类型和网络类型等。分析的结果存储在net对象中。
cor<-stats::cor
#重新将"stats"包中的cor函数赋给cor变量，覆盖之前的"WGCNA"包中的cor函数。
#save(net,file="net_mes.Rdata")
#load("net_mef.Rdata")
table(net$colors)
#输出每个模块的基因数目的频数分布表。
# length(intersect(colnames(exp_new)[which(net$colors==0)],tf))
#输出模块0中与载体基因(tf)的交集长度。

net1<-net
for(i in 1:length(table(net$colors))){
  if(table(net$colors)[i]>1000){
    
    geneinmodule<-which(net$colors==names(table(net$colors)[i]))
    mdatExpr<-exp_new[,geneinmodule]
    #mdatExpr<-datExpr[geneinmodule,]原始代码但越界
    #adjacency1 = adjacency(mdatExpr, power =sft$powerEstimate,type = "signed")
    adjacency1 = adjacency(mdatExpr, power =sft$powerEstimate,type = "signed")
    TOM1 = TOMsimilarity(adjacency1);
    dissTOM1 = 1-TOM1
    geneTree = hclust(as.dist(dissTOM1));
    minModuleSize = 500;
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM1,pamStage = FALSE,
                                deepSplit = 2, pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize );
    if(length(table(dynamicMods))>1)
      for(j in 1:length(table(dynamicMods)))
        net1$colors[geneinmodule[which(dynamicMods==names(table(dynamicMods)[j]))]]<-max(as.numeric(names(table(net1$colors))))+1
  }
}
#嵌套在循环内的代码对于包含大量基因的模块进行进一步处理和分析。其中，geneinmodule获取属于当前模块的基因索引，mdatExpr是根据这些基因生成的新的基因表达矩阵，adjacency1计算这个新矩阵的邻接矩阵，TOM1计算邻接矩阵的TOM相似性，geneTree是基于TOM相似性计算的基因聚类树，minModuleSize是指定的最小模块大小，dynamicMods利用动态树切分算法进行模块划分。


table(net1$colors)
#输出更新后的模块中基因数目的频数分布表。
#save(net1,file="net_mef1.Rdata")
#load("net_mef1.Rdata")
for (k in 0:6)
  print(length(intersect(colnames(exp_new)[which(net1$colors==k)],tf)))
#输出模块14中与载体基因(tf)的交集长度。

nGenes <- ncol(exp_new)#将datExpr矩阵的列数赋给nGenes，表示基因的数量。
nSamples <- nrow(exp_new)#将datExpr矩阵的行数赋给nSamples，表示样本的数量。
design <- as.data.frame(model.matrix(~0+ type_new))
#根据type变量创建一个适用于线性模型的设计矩阵，其中将其转化为一个数据框。
design<-design[,paste("type_new",names(table(type)),sep="")]
#根据模式type在设计矩阵中选择特定的列。
#moduleColors<-labels2colors(dynamicMods)
moduleColors<-labels2colors(net1$colors)
# plotDendroAndColors(geneTree, moduleColors, "Dynamic Tree Cut",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05,
#                     main = "Gene dendrogram and module colors")
#利用net1$colors中的模块标签（颜色）创建一个用于模块颜色的向量，即根据模块标签给每个样本分配相应的颜色。
# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(exp_new, moduleColors)$eigengenes
#计算给定单个数据集中模块的特征基因(第一主成分)
#根据给定的单个数据集和模块颜色标签，计算模块的特征基因（即每个模块的第一主成分）
MEs <- orderMEs(MEs0);
##不同颜色的模块的ME值矩 (样本vs模块)，
#将模块特征基因（MEs0）按照模块的顺序（模块颜色标签排序）重新排序。
moduleTraitCor <- cor(MEs, design , use = "p");
#计算模块特征基因（MEs）与设计矩阵（design）之间的相关性，使用皮尔逊相关系数
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
#基于相关性计算模块特征基因与设计矩阵之间的P值。

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
#生成一个文本矩阵，包含相关性值和P值的格式化字符串。
par(mar = c(5, 10, 4, 2) + 0.1)  # 设置边距，增加左边距以显示y轴标签

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(table(type)),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 1.5,  # 增大方块内数字的大小
  cex.lab.x = 1,  # 缩小x轴标签的大小
  cex.lab.y = 1,  # 缩小y轴标签的大小
  zlim = c(-1, 1),  # 设置颜色映射范围
  yColorWidth = 0.01,
  xLabelsPosition = "bottom",
  xColorWidth = 0.5,
  main = paste("Module-trait relationships")
)
#绘制带有标签、颜色和相关性值的热图，显示模块特征基因与设计矩阵之间的相关性。
table(moduleColors)#输出模块颜色向量中每个模块的基因数量的频数分布表。
#################################GENIE3################################################
toptf_list<-list()
all_transcription_factors <- c()  # 创建一个向量来存储所有筛选出来的转录因子
allmodule<-c()
for(i in 1:length(table(type_new))){
  cell_type = as.data.frame(design[,i])
  #将design数据的第i列转换成数据框cell_type
  names(cell_type) = names(table(type_new))[i]
  #将cell_type数据框的列名设置为table(type_new)中对应列名的第i个
  sig_module<-signif(moduleTraitCor, 2)[,i]
  #从moduleTraitCor中提取与第i列相关的显著模块
  module <- substring(names(sig_module[which(sig_module>=0.5)]), 3)
  #根据模块与第i列相关的相关性阈值（大于等于0.5），提取显著模块的名称。
  if(length(module)==0)
    module <- substring(names(sig_module[which.max(sig_module)]), 3)
  #将与第i列相关性最大的模块作为显著模块。
  modNames = substring(names(MEs), 3)
  #从MEs中提取所有模块的名称。
  for(j in 1:length(module))#遍历显著模块的每个元素。
  {
    if(!any(allmodule==module[j])){
      column = match(module[j], modNames);
      #在modNames中找到与第j个显著模块名称匹配的列的索引。
      moduleGenes = moduleColors==module[j];#确定与当前模块匹配的基因。
      exprMatr<-t(exp_new[,moduleGenes])
      #从exp_new中提取所有行和与当前模块匹配的列，得到一个特定模块的转置矩阵。
      corMatr<-cor(t(exprMatr))
      #计算特定模块的转置矩阵的相关系数矩阵。
      set.seed(42)
      #regulator<-intersect(alltf,rownames(exprMatr))源代码
      regulator <- rownames(exprMatr)
      #找到alltf和exprMatr中细胞特征的共同转录因子。
      if(length(regulator)>1){
        weightMat<-GENIE3(exprMatr,regulators=regulator)
        #使用GENIE3算法计算基因表达与转录因子之间的因果关系权重。
        #save(weightMat,file=paste("./weightMat",i,"_",j,".Rdata",sep=""))
        linklist<-getLinkList(weightMat,threshold=0.0245)
        #根据阈值筛选出重要的连接。
        colnames(linklist)<-c("Source","Target","Weight")
        linklist_csv <- paste0("./linkList", i, "_", j, ".csv")
        write.csv(linklist, file = linklist_csv, row.names = FALSE)
        topn<-min(50,length(table(linklist[,1])))
        #确定要选择的前n个转录因子（最多50个）。
        toptf<-names(sort(table(linklist[,1]),decreasing = T)[1:topn])
        #从连接列表中选择出现次数最多的前n个转录因子，并将它们的名称存储在toptf中。
        all_transcription_factors <- unique(c(all_transcription_factors, toptf))
        # 收集转录因子
        toptf_list[[(i+j-1)]]<-toptf
        #将选定的转录因子列表toptf添加到toptf_list中的合适位置。
        #length(intersect(toptf,tf))
        if(length(toptf)>1){
          names(toptf)<-rep(names(cell_type),length(toptf))
          #将转录因子列表的名称设置为与当前单元类型相对应的名称。
          # pdfname<-paste("./venn",i,"-",j,".pdf",sep="")
          # intersect_TF<-venn.diagram(list(TF_CONTRAST=tf,TF_SELECT=na.omit(toptf)),
          #                            filename=NULL,col="transparent",margin=0.1,
          #                            main.cex = 2,fill=c('#668B8B','#00868B'),disable.logging=F)
          #根据转录因子列表和tf创建一个韦恩图。
          # pdf(pdfname)
          # grid.draw(intersect_TF)
          # dev.off()
          #heatmap
          exp_tf<-t(exp_new)[colnames(exp_new)%in%toptf,]
          #从exp_new中提取与所选转录因子名称相对应的列。
          group_sample<-data.frame(type_new)
          #创建一个名为group_sample的数据框，其中包含type_new的数据。
          rownames(group_sample)<-colnames(exp_tf)
          #将exp_tf的列名设置为group_sample的行名。
          colnames(group_sample)<-"type_new"
          #将group_sample的列名设置为"type_new"。
          
          
          group_sample$type_new <- as.factor(group_sample$type_new)
          sorted_index <- order(group_sample$type_new)
          
          # 按照排序后的索引重排exp_tf的列
          exp_tf_sorted <- exp_tf[, sorted_index]
          group_sample_sorted <- group_sample[sorted_index, , drop=FALSE]
          
          # 绘制热图，关闭列的聚类
          pdfname<-paste("./pheatmap",i,"-",j,".pdf",sep="")
          pdf(pdfname)
          pheatmap(exp_tf_sorted,
                   scale = "row",
                   color = colorRampPalette(c("blue","white","red"))(50),
                   display_numbers = F,
                   cluster_cols = F, # 关闭列的聚类功能，因为已经手动排序
                   cluster_rows = T,
                   show_colnames = F,
                   show_rownames = T,
                   annotation_col = group_sample_sorted) # 使用排序后的注释信息
          dev.off()
          cat("i=",i,"j=",j,"intersect=",intersect(tf,na.omit(toptf)),"\n")
        }
      }
    }
  }
  allmodule<-c(allmodule,module) 
}
dev.off()
exp_tf<-t(exp_new)[colnames(exp_new)%in%all_transcription_factors,]
#从exp_new中提取与所选转录因子名称相对应的列。
group_sample<-data.frame(type_new)
#创建一个名为group_sample的数据框，其中包含type_new的数据。
rownames(group_sample)<-colnames(exp_tf)
#将exp_tf的列名设置为group_sample的行名。
colnames(group_sample)<-"type_new"
#将group_sample的列名设置为"type_new"。


group_sample$type_new <- as.factor(group_sample$type_new)
sorted_index <- order(group_sample$type_new)

# 按照排序后的索引重排exp_tf的列
exp_tf_sorted <- exp_tf[, sorted_index]
group_sample_sorted <- group_sample[sorted_index, , drop=FALSE]

ann_colors <- list(
  type_new = c("fetal_10" = "#DC143C",  # 红色
               "fetal_17" = "#0000FF",  # 蓝色
               "hHep" = "#20B2AA",     # 绿色
               "hLivEC" = "#FFA500")    # 橙色
)

# 进一步调整颜色梯度，使红色和蓝色更加纯粹
color_palette <- colorRampPalette(c("#00008B", "#0000FF", "white", "#FF0000", "#8B0000"))(100)

# 绘制热图，关闭列的聚类
pheat<-pheatmap(exp_tf_sorted,
                scale = "row",
                color = color_palette,
                display_numbers = F,
                cluster_cols = F, # 关闭列的聚类功能，因为已经手动排序
                cluster_rows = T,
                show_colnames = F,
                show_rownames = F,
                annotation_col = group_sample_sorted,
                annotation_colors = ann_colors)
# ggsave("热图.pdf", pheat, width = 175, height = 144, units = "mm",dpi=300)
# save(toptf_list,file="toptf_list.Rdata")
# load("toptf_list.Rdata")
