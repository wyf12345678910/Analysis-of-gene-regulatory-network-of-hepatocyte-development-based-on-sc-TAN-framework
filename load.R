# BiocManager::install("WGCNA")
# BiocManager::install("limma")
# BiocManager::install("biomaRt")
# BiocManager::install("GENIE3")
# BiocManager::install("monocle")
# BiocManager::install("GO.db")
# BiocManager::install("impute")
# BiocManager::install("preprocessCore")
# BiocManager::install("qlcMatrix")
# BiocManager::install("sf")
# BiocManager::install("proxy")
# BiocManager::install("grid")
# BiocManager::install("futile.logger")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# #下载monocle依赖的其它包
# BiocManager::install("monocle",force=T)
# BiocManager::install("factoextra")
# BiocManager::install("randomForest")
# BiocManager::install("LEAP")
# BiocManager::install("Seurat")
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
# exp = t(exp_liver[sorted_indices, ])
exp2=t(exp_liver)
TFlist=intersect(unique(TFlist1[,2]),rownames(exp2))
tf=intersect(unique(TFlist1[,2]),rownames(exp2))
length(intersect(unique(TFlist),rownames(exp2)))
# type <- sorted_type
# TFlist=intersect(unique(TFlist2[,1]),rownames(exp))
# tf=intersect(unique(TFlist2[,1]),rownames(exp))
# length(intersect(unique(TFlist),rownames(exp)))
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
# 
# del_gen<-c()
# for(i in 1:dim(exp1)[1]){
#   if(var(as.numeric(exp1[i,]))<0.9)
#     #如果某个基因在所有细胞中的表达值的方差（使用var()函数计算）小于0.9，则执行以下代码块。
#     del_gen<-c(del_gen,i)#将符合条件的基因索引添加到del_gen向量中
# }
# length(del_gen)#4015
# exp2<-exp1[-del_gen,]
# dim(exp2)#10272   333

# 假设exp是基因表达矩阵，行为基因，列为样本
# 设定阈值为0.5 CPM，在至少5%的样本中有表达
# threshold <- 0.5
# min_samples <- ncol(exp) * 0.05
# 
# # 筛选出至少在10%的样本中表达量大于1 CPM的基因
# filtered_genes <- rowSums(exp >= threshold) >= min_samples
# exp2 <- exp[filtered_genes, ]
# dim(exp2)
#对细胞进行基因表达的标准化处理(是为了消除不同细胞之间的技术差异和批次效应，以便更准确地比较和分析基因表达水平。)
#FPKM to TPM,FPKM和TPM均为常见的两种标准化方法
# fpkmToTpm <- function(fpkm)
# {
#   exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
# }
# 定义了一个名为fpkmToTpm的函数，用于将FPKM（Fragments Per Kilobase of transcript per Million mapped reads）转换为TPM(Transcripts Per Million)
# tpms <- apply(exp2,2,fpkmToTpm)
#应用fpkmToTpm函数将基因表达矩阵exp2中的每列进行转换，得到TPM值。
#tpms<-exp2
# m.mad <- apply(exp,1,mad)
# exp2<- exp[which(m.mad >max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.005)),]
# dim(exp2)
# length(intersect(unique(tf),rownames(exp2)))
# m.mad <- apply(exp,1,mad)
# dim(exp)
# #筛选时不设置mad最小值，直接使用前75%的基因或者探针。
# exp2 <- exp[which(m.mad >  max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
dim(exp2)
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

# mexp<-cbind(rowMeans(exp2[,which(type==names(table(type))[1])],na.rm=T),
#             rowMeans(exp2[,which(type==names(table(type))[2])],na.rm=T),
#             rowMeans(exp2[,which(type==names(table(type))[3])],na.rm=T),
#             rowMeans(exp2[,which(type==names(table(type))[4])],na.rm=T))
# 
# mexp<-log(mexp+1)
# 
# colnames(mexp)<-names(table(type))
# 
# datExpr2 <- FindVariableFeatures(mexp, selection.method = "mvp")
# #使用"mvp"方法，根据基因表达矩阵mexp找出具有变异性的基因。结果保存在datExpr2中。
# datExpr2<-datExpr2[order(datExpr2[,3],decreasing=T),]
# #据第三列（变异性）降序对datExpr2进行排序
# #save(datExpr2,file="datExpr_mvp.Rdata")
# #load("datExpr_mvp.Rdata")
# 
# top<-datExpr2[c(1:10300),]
# 
# gene_totol<-union(rownames(datExpr1),rownames(top))
# 
# data_total<-exp[rownames(exp)%in%gene_totol,]
# dim(data_total)
# datExpr<-data_total
#save(datExpr,file="data_total.Rdata")
#load("./data_total.Rdata")
# length(intersect(unique(tf),rownames(datExpr)))

# sampleTree = hclust(dist(datExpr), method = "average")
# old_par <- par()  # 保存当前的图形参数设置
# par(mar = c(2, 4, 4, 2))  # 减少下方的空白
# # 计算树状图中高度的最小和最大值
# heights <- sampleTree$height
# ylim_range <- c(min(heights), max(heights))
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", ylab="", labels=F, ylim=ylim_range)
# 
# 
# 
# # p1<-ggplot(liver, aes(x = liver[,3], y = liver[,4])) +
# #   geom_point(aes(color = factor(type)), size = 1.5) +  # 假设有一个分类变量cluster
# #   theme_minimal() +
# #   labs(title = "t-SNE plot of cells",
# #        x = "t-SNE Dimension 1",
# #        y = "t-SNE Dimension 2",
# #        color = "type") +
# #   theme(legend.position = "right")
# # p1
# 
# p2 <- ggplot(liver, aes(x = liver[,3], y = liver[,4])) +
#   geom_point(aes(color = factor(type)), size = 1.5) +  
#   theme_minimal() +
#   labs(title = "t-SNE plot of cells",
#        x = "t-SNE Dimension 1",
#        y = "t-SNE Dimension 2",
#        color = "type") +
#   theme(legend.position = "right",
#         panel.grid.major = element_blank(),  # 删除主要网格线
#         panel.grid.minor = element_blank(),  # 删除次要网格线
#         panel.border = element_blank(),       # 删除面板边框
#         axis.line = element_line(colour = "black")) +  # 添加坐标轴线，增强可读性
#   scale_color_manual(values = c("#FA7F6F",  "#8ECFC9","#82B0D2" ,"#E7DAD2"))  # 设置颜色
# p2
# 
# # deg$type <- ifelse(deg$logFC > 2 & deg$adj.P.Val < 0.01, "up",
# #                        ifelse(deg$logFC < -2 & deg$adj.P.Val < 0.01, "down", "not-sig")
# # )
# # deg$gene <- rownames(deg) 
# # p2 <- ggplot(deg, aes(logFC, -log10(adj.P.Val)))+
# #   geom_point(aes(color=type))+
# #   scale_color_manual(values = c("blue","black","red"))+
# #   geom_hline(yintercept = -log10(0.01),linetype=2)+
# #   geom_vline(xintercept = c(-2,2), linetype=2)+
# #   geom_text_repel(data = subset(deg, abs(logFC) > 4), 
# #                   aes(label=gene),col="black",alpha = 0.8)+
# #   theme_bw()
# # p2
# 
# plots <- list()
# info<-list()
# titles <- c("Fetal_10", "Fetal_17", "hHep", "hLivEC")  # 不同的标题名字
# for (i in 1:4) {
#   deg <- deg_list[[i]]
#   deg$type <- ifelse(deg$logFC > 1 & deg$adj.P.Val < 0.05, "up",
#                      ifelse(deg$logFC < -1 & deg$adj.P.Val < 0.05, "down", "not-sig"))
#   deg$gene <- rownames(deg)
#   
#   # 为up和down类型的基因分别选出最大的5个
#   top_up <- deg[deg$type == "up",][order(deg[deg$type == "up", "logFC"], decreasing = TRUE)[1:5], ]
#   top_down <- deg[deg$type == "down",][order(deg[deg$type == "down", "logFC"], decreasing = FALSE)[1:5], ]
#   info[[i]] <- rbind(deg[deg$type == "up",], deg[deg$type == "down",])
#   # 合并top_up和top_down
#   top_genes <- rbind(top_up, top_down)
#   
#   # 创建火山图
#   plots[[i]] <- ggplot(deg, aes(logFC, -log10(adj.P.Val))) +
#     geom_point(aes(color = type)) +
#     scale_color_manual(values = c("#8983BF", "#E7DAD2", "#FA7F6F")) +
#     geom_hline(yintercept = -log10(0.05), linetype = 2) +
#     geom_vline(xintercept = c(-1, 1), linetype = 2) +
#     geom_text_repel(data = top_genes,
#                     aes(label = gene), color = "black", alpha = 0.8,
#                     max.overlaps = 10)+  # 增加允许的重叠数量
#     theme_bw() +
#     xlab(expression(paste("log"[2], "FC"))) +
#     ylab(expression(paste("log"[10], "adj.P.Val")))+
#     labs(title = titles[i])+  # 使用不同的标题
#     theme(plot.title = element_text(size = 12),  # 缩小标题字号
#           legend.text = element_text(size = 8))  # 缩小标识字号
#   
#   # 打印图形
# }
# grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], nrow = 2, ncol = 2)
# 
# part1_names <- select_list[[1]]
# part2_names <- select_list[[2]]
# part3_names <- select_list[[3]]
# part4_names <- select_list[[4]]
# vennlist_names <- list(Fetal_10= part1_names, Fetal_17 = part2_names, hHep = part3_names, hLivEC = part4_names)
# # 使用 ggvenn 绘制韦恩图
# venn.plot <- venn.diagram(
#   x = vennlist_names,
#   category.names = names(vennlist_names),
#   filename = NULL,
#   output = TRUE,
#   height = 480,
#   width = 480,
#   resolution = 300,
#   col = "transparent",
#   lty = "blank",
#   fill = c(brewer.pal(8, 'Set2')[2:5]),
#   cat.col = c(brewer.pal(8, 'Set2')[2:5]),
#   cat.cex = 1.5,
#   cex = 1.5,
#   fontface = "bold",
#   cat.fontface = "bold"
# )
# # 绘制图形
# grid.draw(venn.plot)
# 
