# sampleTree = hclust(dist(datExpr), method = "average")
# old_par <- par()  # 保存当前的图形参数设置
# par(mar = c(2, 4, 4, 2))  # 减少下方的空白
# # 计算树状图中高度的最小和最大值
# heights <- sampleTree$height
# ylim_range <- c(min(heights), max(heights))
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", ylab="", labels=F, ylim=ylim_range)


setwd("D:/R/R工作区/cell/对比/新建文件夹")

# p1<-ggplot(liver, aes(x = liver[,3], y = liver[,4])) +
#   geom_point(aes(color = factor(type)), size = 1.5) +  # 假设有一个分类变量cluster
#   theme_minimal() +
#   labs(title = "t-SNE plot of cells",
#        x = "t-SNE Dimension 1",
#        y = "t-SNE Dimension 2",
#        color = "type") +
#   theme(legend.position = "right")
# p1

###tSNE
library(ggplot2)

# 创建 ggplot 图
p2 <- ggplot(liver, aes(x = liver[,3], y = liver[,4])) +
  geom_point(aes(color = factor(type)), size = 1.5) +  
  theme_minimal(base_size = 18) +  # 设置基础字体大小
  labs(title = "t-SNE plot of cells",
       x = "t-SNE Dimension 1",
       y = "t-SNE Dimension 2",
       color = "type") +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),  # 删除主要网格线
    panel.grid.minor = element_blank(),  # 删除次要网格线
    panel.border = element_blank(),       # 删除面板边框
    axis.line = element_line(colour = "black"),
    plot.title = element_text(size = 22, face = "bold"),  # 增大标题字体大小
    axis.title = element_text(size = 18),  # 增大轴标题字体大小
    axis.text = element_text(size = 16),  # 增大轴文本字体大小
    legend.title = element_text(size = 18),  # 增大图例标题字体大小
    legend.text = element_text(size = 16)  # 增大图例文本字体大小
  ) +
  scale_color_manual(values = c("#DC143C","#0000FF","#20B2AA","#FFA500"))  # 设置颜色

# 显示图形
print(p2)
# ggsave("tSNE.pdf", p2, width = 6.89, height = 5.67, units = "in",dpi=300)
# deg$type <- ifelse(deg$logFC > 2 & deg$adj.P.Val < 0.01, "up",
#                        ifelse(deg$logFC < -2 & deg$adj.P.Val < 0.01, "down", "not-sig")
# )
# deg$gene <- rownames(deg) 
# p2 <- ggplot(deg, aes(logFC, -log10(adj.P.Val)))+
#   geom_point(aes(color=type))+
#   scale_color_manual(values = c("blue","black","red"))+
#   geom_hline(yintercept = -log10(0.01),linetype=2)+
#   geom_vline(xintercept = c(-2,2), linetype=2)+
#   geom_text_repel(data = subset(deg, abs(logFC) > 4), 
#                   aes(label=gene),col="black",alpha = 0.8)+
#   theme_bw()
# p2

###火山图
# plots <- list()
# info <- list()
# titles <- c("Fetal_10", "Fetal_17", "hHep", "hLivEC")  # 不同的标题名字
# 
# # 找到所有子图中 y 轴的最大值
# max_y <- max(sapply(deg_list, function(deg) {
#   max(-log10(deg$adj.P.Val), na.rm = TRUE)
# }))
# 
# # 存储每组上调和下调基因的数量
# up_counts <- numeric(length(titles))
# down_counts <- numeric(length(titles))
# 
# for (i in 1:4) {
#   deg <- deg_list[[i]]
#   deg$type <- ifelse(deg$logFC > 1 & deg$adj.P.Val < 0.05, "up",
#                      ifelse(deg$logFC < -1 & deg$adj.P.Val < 0.05, "down", "not-sig"))
#   deg$gene <- rownames(deg)
#   
#   # 计算上调和下调基因数量
#   up_counts[i] <- sum(deg$type == "up")
#   down_counts[i] <- sum(deg$type == "down")
#   
#   # 为up和down类型的基因分别选出最大的5个
#   top_up <- deg[deg$type == "up",][order(deg[deg$type == "up", "logFC"], decreasing = TRUE)[1:5], ]
#   top_down <- deg[deg$type == "down",][order(deg[deg$type == "down", "logFC"], decreasing = FALSE)[1:5], ]
#   info[[i]] <- rbind(deg[deg$type == "up",], deg[deg$type == "down",])
#   
#   # 合并top_up和top_down
#   top_genes <- rbind(top_up, top_down)
#   
#   # 创建火山图
#   p <- ggplot(deg, aes(logFC, -log10(adj.P.Val))) +
#     geom_point(aes(color = type)) +
#     scale_color_manual(values = c("up" = "#FA7F6F", "down" = "#8983BF", "not-sig" = "#E7DAD2")) +
#     geom_hline(yintercept = -log10(0.05), linetype = 2) +
#     geom_vline(xintercept = c(-1, 1), linetype = 2) +
#     geom_text_repel(data = top_genes,
#                     aes(label = gene), color = "black", alpha = 0.8,
#                     max.overlaps = 10) +
#     coord_cartesian(ylim = c(0, max_y)) +  # 统一纵坐标范围
#     theme_bw() +
#     xlab(expression(paste("log"[2], "FC"))) +
#     ylab(expression(paste("log"[10], "adj.P.Val"))) +
#     labs(title = titles[i]) +  # 使用不同的标题
#     theme(plot.title = element_text(size = 12),  # 缩小标题字号
#           legend.text = element_text(size = 8))  # 缩小标识字号
#   
#   plots[[i]] <- p
# }
# 
# # 提取第一个图的图例
# legend <- cowplot::get_legend(plots[[1]])
# 
# # 移除其他图的图例
# for (i in 1:4) {
#   plots[[i]] <- plots[[i]] + theme(legend.position = "none")
# }
# 
# # 将四个图并排显示，并在最后添加图例
# p3 <- plot_grid(plot_grid(plotlist = plots, nrow = 2, ncol = 2), 
#                 legend, 
#                 ncol = 2, 
#                 rel_heights = c(1, 0.2))
# 
# print(p3)
library(cowplot)
library(ggplot2)
library(ggrepel)

plots <- list()
info <- list()
titles <- c("Fetal_10", "Fetal_17", "hHep", "hLivEC")  # 不同的标题名字

# 存储每组上调和下调基因的数量
up_counts <- numeric(length(titles))
down_counts <- numeric(length(titles))

for (i in 1:4) {
  deg <- deg_list[[i]]
  deg$type <- ifelse(deg$logFC > 1 & deg$adj.P.Val < 0.05, "up",
                     ifelse(deg$logFC < -1 & deg$adj.P.Val < 0.05, "down", "not-sig"))
  deg$gene <- rownames(deg)
  
  # 计算上调和下调基因数量
  up_counts[i] <- sum(deg$type == "up")
  down_counts[i] <- sum(deg$type == "down")
  
  # 为 up 和 down 类型的基因分别选出最大的 5 个
  top_up <- deg[deg$type == "up",][order(deg[deg$type == "up", "logFC"], decreasing = TRUE)[1:5], ]
  top_down <- deg[deg$type == "down",][order(deg[deg$type == "down", "logFC"], decreasing = FALSE)[1:5], ]
  info[[i]] <- rbind(deg[deg$type == "up",], deg[deg$type == "down",])
  
  # 合并 top_up 和 top_down
  top_genes <- rbind(top_up, top_down)
  
  # 设置不同的y轴上限
  y_axis_limit <- if (i == 3) 150 else if (i == 4) 40 else 100
  
  # 创建火山图
  p <- ggplot(deg, aes(logFC, -log10(adj.P.Val))) +
    geom_point(aes(color = type)) +
    scale_color_manual(values = c("up" = "#FA7F6F", "down" = "#8983BF", "not-sig" = "#E7DAD2")) +
    geom_hline(yintercept = -log10(0.05), linetype = 2) +
    geom_vline(xintercept = c(-1, 1), linetype = 2) +
    geom_label_repel(data = top_genes,
                     aes(label = gene),
                     nudge_x = 0.1,  # 调整水平偏移量
                     nudge_y = 0.1,  # 调整垂直偏移量
                     color = "black", alpha = 0.8,
                     max.overlaps = 20,  # 增加最大重叠数量
                     segment.color = "grey",  # 细线颜色
                     segment.size = 0.6,  # 细线尺寸，增加粗细
                     force = 5,  # 增加排斥力
                     box.padding = 0.5,  # 增加标签框的填充
                     point.padding = 0.5) +  # 增加点的填充
    coord_cartesian(ylim = c(0, y_axis_limit)) +  # 调整y轴范围
    scale_y_continuous(breaks = seq(0, y_axis_limit, by = if (y_axis_limit == 50) 10 else 50)) +  # 设置 y 轴范围和刻度
    theme_bw() +
    xlab(expression(paste("log"[2], "FC"))) +
    ylab(expression(paste("-log"[10], "adj.P.Val"))) +  # 放大 y 轴标签
    labs(title = titles[i]) +  # 使用不同的标题
    theme(plot.title = element_text(size = 15),  # 缩小标题字号
          legend.text = element_text(size = 12),  # 调整标识字号
          legend.title = element_text(size = 12),  # 调整图例标题
          axis.title.y = element_text(size = 16),  # 放大 y 轴标题
          axis.title.x = element_text(size = 16),  # 放大 x 轴标题
          legend.key.size = unit(1, "lines"))  # 调整图例颜色圈的大小
  
  plots[[i]] <- p
}

# 提取第一个图的图例
legend <- cowplot::get_legend(plots[[1]])

# 移除其他图的图例
for (i in 1:4) {
  plots[[i]] <- plots[[i]] + theme(legend.position = "none")
}

# 将四个图并排显示，并在右侧添加图例
p3 <- plot_grid(plotlist = plots, nrow = 1, ncol = 4, align = 'hv', axis = 'tblr')  # 将四个图并排显示

# 将图例添加到组合图中
p3 <- plot_grid(p3, legend, ncol = 2, rel_widths = c(1, 0.1))  # 调整图例的相对宽度

# 打印组合图
print(p3)
# ggsave("火山图.pdf", p3, width = 12.89, height = 7.67, units = "in",dpi=300)

####柱状图
data <- as.matrix(data.frame(Fetal_10 = c(up_counts[1], down_counts[1]),         
                             Fetal_17 = c(up_counts[2], down_counts[2]),
                             hHep = c(up_counts[3], down_counts[3]),
                             hLivEC = c(up_counts[4], down_counts[4])))
row.names(data) <- c("up","down") 

# 设置更专业的颜色方案和图表主题
colors <- c("red", "blue")  # 选用清晰区分的颜色
theme_set(theme_minimal(base_size = 14))  # 使用更大的基础字体大小，清晰的主题

# 绘制柱状图
# 增加右侧边距
par(mar = c(5, 5, 5, 15))  # 增加右侧边距以放置图例

# 重新绘制条形图
bar_heights <- barplot(data,
                       col = colors,
                       beside = TRUE,
                       ylim = c(0, max(data) * 1.2),  # 留出更多空间防止数字超出
                       las = 1,  # 坐标轴标签水平显示
                       cex.names = 1.2)  # x轴标签大小

# 添加数字
for (j in 1:ncol(data)) {
  for (i in 1:nrow(data)) {
    text(x = bar_heights[i, j], y = data[i, j] + 0.01 * max(data), labels = data[i, j], pos = 3, cex = 1.2)
  }
}

# 添加图例并调整其大小
legend("topright",
       inset = c(0, 0),  # 将图例向右移动
       legend = c("Up", "Down"),
       fill = colors,
       box.lty = 1,  # 图例框的线型
       box.lwd = 0.2,  # 缩小图例框的线条宽度
       cex = 1.2,  # 缩小图例字体
       y.intersp = 0.6,  # 缩小图例项之间的垂直间距
       x.intersp = 0.5)  # 缩小图例项之间的水平间距

# 添加描述性标签并调整字体大小
title(main = "Differential Gene Expression",
      cex.main = 2,  # 调整标题大小
      xlab = "Type",
      cex.lab = 1.5,  # 调整x轴和y轴标签大小
      ylab = "Number of DEGs",
      cex.axis = 1.4)  # 调整坐标轴刻度标签大小

# 添加轻微的网格线以增强可读性
grid(nx = NA, ny = NULL)


#ggsave2("柱状图.pdf", p4, width = 175, height = 144, units = "mm",dpi=300)


###Venn图
part1_names <- select_list[[1]]
part2_names <- select_list[[2]]
part3_names <- select_list[[3]]
part4_names <- select_list[[4]]
vennlist_names <- list(Fetal_10 = part1_names, Fetal_17 = part2_names, hHep = part3_names, hLivEC = part4_names)

# 使用 VennDiagram 包绘制韦恩图
venn.plot <- venn.diagram(
  x = vennlist_names,
  category.names = names(vennlist_names),
  filename = NULL,
  output = TRUE,
  height = 480,
  width = 480,
  resolution = 300,
  col = "transparent",
  lty = "blank",
  fill = c(brewer.pal(8, 'Set2')[2:5]),
  cat.col = c(brewer.pal(8, 'Set2')[2:5]),
  cat.cex = 1.5,
  cex = 1.5,
  fontface = "bold",
  cat.fontface = "bold"
)

# 绘制韦恩图
grid.newpage()
p5 <- grid.draw(venn.plot)
print(p5)
# ggsave("Venn.pdf", p5, width = 6.89, height = 5.67, units = "in",dpi=300)

# # 循环遍历 info 列表中的每个子列表，输出为单独的 CSV 文件
# for (i in 1:length(info)) {
#   # 动态生成文件名，例如 info_part_1.csv, info_part_2.csv, etc.
#   file_name <- paste0("info_part_", i, ".csv")
#   
#   # 将每个子列表输出为 CSV 文件，且不保存行名
#   write.csv(info[[i]], file = file_name, row.names = FALSE)
# }
# output_dir <- "D:/R/R工作区/cell/对比/新建文件夹/"
# 
# # 输出每个子列表为单独的 CSV 文件到指定路径
# for (i in 1:length(info)) {
#   # 动态生成文件名，并包含完整路径
#   file_name <- paste0(output_dir, "info_part_", i, ".csv")
#   
#   # 保存每个子列表为 CSV 文件，且不保存行名
#   write.csv(info[[i]], file = file_name, row.names = FALSE)
# }


#GO富集分析
library(openxlsx)#读取.xlsx文件
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
dev.off()
# toptf_main_list<-list()
# toptf_main_list[[1]]<-c("AHSP","ALAS2","ANK1","BSG","DNAJC9","GLRX5","GYPA","GYPB","HBA1","HBA2","HBG1",
#                         "HBG2","HBQ1","HMBS","HMGB2","KLF1","MYL4","PRDX2","RFESD","RHAG","SLC25A39","UROD")
# toptf_main_list[[2]]<-c("AIMP2","AK2","ALG3","C1QBP","CCT5","EIF6","ENO1","FAM136A"
#                         ,"MIF","NPM1","PHAX","PITPNB","PPP1R14B","PSMB2","PSMG1","RPS7","EIF4EBP1")
# toptf_main_list[[3]]<-c("EIF3K","EIF5AL1","LSM7","NOL7","NUTF2","POLR2G","PPIA","RPS2","SNRPC","TUBA1B","TUBA1C","XRCC6")
# toptf_main_list[[4]]<-c("ADH1B","ADH4","CES1","CYP2C8","CYP2C9","HRG","ITIH4","SULT2A1")
# GO_database <- 'org.Hs.eg.db' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
# KEGG_database <- 'hsa' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html
# gene<-list()
# info1<-c()
# for (i in 1:4){
# info1<-info[[i]][toptf_main_list[[i]],c(8,1,4,5)]
# gene[[i]]=bitr(info1$gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
# GO<-enrichGO( gene[[i]]$ENTREZID,#GO富集分析
#               OrgDb = GO_database,
#               keyType = "ENTREZID",#设定读取的gene ID类型
#               ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
#               pvalueCutoff = 0.05,#设定p值阈值
#               qvalueCutoff = 0.05,#设定q值阈值
#               readable = T)
# KEGG<-enrichKEGG(gene[[i]]$ENTREZID,#KEGG富集分析
#                  organism = KEGG_database,
#                  pvalueCutoff = 0.05,
#                  qvalueCutoff = 0.05)
# print(barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free"))#柱状图
# print(barplot(KEGG,showCategory = 40,title = 'KEGG Pathway'))
# print(dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free"))#点状图
# print(dotplot(KEGG))# print(GO)
# # 提取GO分析结果
# genedata<-data.frame(ID=info1$gene,logFC=info1$logFC)
# GO_results <- as.data.frame(GO)
# GOplotIn_BP<-GO[1:10,c(2,3,7,9)] #提取GO富集BP的前10行,提取ID,Description,p.adjust,GeneID四列
# GOplotIn_CC<-GO[39:43,c(2,3,7,9)]#提取GO富集CC的前10行,提取ID,Description,p.adjust,GeneID四列
# GOplotIn_MF<-GO[44:63,c(2,3,7,9)]#提取GO富集MF的前10行,提取ID,Description,p.adjust,GeneID四列
# GOplotIn_BP$geneID <-str_replace_all(GOplotIn_BP$geneID,'/',',') #把GeneID列中的’/’替换成‘,’
# GOplotIn_CC$geneID <-str_replace_all(GOplotIn_CC$geneID,'/',',')
# GOplotIn_MF$geneID <-str_replace_all(GOplotIn_MF$geneID,'/',',')
# names(GOplotIn_BP)<-c('ID','Term','adj_pval','Genes')#修改列名,后面弦图绘制的时候需要这样的格式
# names(GOplotIn_CC)<-c('ID','Term','adj_pval','Genes')
# names(GOplotIn_MF)<-c('ID','Term','adj_pval','Genes')
# GOplotIn_BP$Category = "BP"#分类信息
# GOplotIn_CC$Category = "CC"
# GOplotIn_MF$Category = "MF"
# circ_BP<-GOplot::circle_dat(GOplotIn_BP,genedata) #GOplot导入数据格式整理
# circ_CC<-GOplot::circle_dat(GOplotIn_CC,genedata)
# circ_MF<-GOplot::circle_dat(GOplotIn_MF,genedata)
# chord_BP<-chord_dat(data = circ_BP,genes = genedata) #生成含有选定基因的数据框
# chord_CC<-chord_dat(data = circ_CC,genes = genedata)
# chord_MF<-chord_dat(data = circ_MF,genes = genedata)
# GOChord(data = chord_BP,#弦图
#         title = 'GO-Biological Process',space = 0.01,#GO Term间距
#         limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
#         lfc.col = c('red','white','blue'), #上下调基因颜色
#         process.label = 10) #GO Term字体大小
# GOChord(data = chord_CC,title = 'GO-Cellular Component',space = 0.01,
#         limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
#         lfc.col = c('red','white','blue'),
#         process.label = 10)
# GOChord(data = chord_MF,title = 'GO-Mollecular Function',space = 0.01,
#         limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
#         lfc.col = c('red','white','blue'),
#         process.label = 10)
# }


toptf_main_list<-list()
toptf_main_list[[1]]<-c("AHSP","ALAS2","ANK1","BSG","DNAJC9","GLRX5","GYPA","GYPB","HBA1","HBA2","HBG1",
                        "HBG2","HBQ1","HMBS","HMGB2","KLF1","MYL4","PRDX2","RFESD","RHAG","SLC25A39","UROD")
toptf_main_list[[2]]<-c("AIMP2","AK2","ALG3","C1QBP","CCT5","EIF6","ENO1","FAM136A"
,"MIF","NPM1","PHAX","PITPNB","PPP1R14B","PSMB2","PSMG1","RPS7","EIF4EBP1")
toptf_main_list[[3]]<-c("EIF3K","EIF5AL1","LSM7","NOL7","NUTF2","POLR2G","PPIA","RPS2","SNRPC","TUBA1B","TUBA1C","XRCC6")
toptf_main_list[[4]]<-c("ADH1B","ADH4","CES1","CYP2C8","CYP2C9","HRG","ITIH4","SULT2A1")
GO_database <- 'org.Hs.eg.db' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
KEGG_database <- 'hsa' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html
gene<-list()
info1<-c()
for (i in 1:3) {
  info1<-info[[i]][toptf_main_list[[i]],c(8,1,4,5)]
  gene[[i]]=bitr(info1$gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  GO<-enrichGO( gene[[i]]$ENTREZID,#GO富集分析
                OrgDb = GO_database,
                keyType = "ENTREZID",#设定读取的gene ID类型
                ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                pvalueCutoff = 0.05,#设定p值阈值
                qvalueCutoff = 0.05,#设定q值阈值
                readable = T)
  print(GO$result)
  KEGG<-enrichKEGG(gene[[i]]$ENTREZID,#KEGG富集分析
                   organism = KEGG_database,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)
  print(barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free"))#柱状图
  print(barplot(KEGG,showCategory = 40,title = 'KEGG Pathway'))
  print(dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free"))#点状图
  print(dotplot(KEGG))# print(GO)
  # 提取GO分析结果
  GO_results <- as.data.frame(GO)

  # 创建一个term数据框，包含GO ID、描述、p值等信息
  terms <- data.frame(Category = GO_results$ONTOLOGY,
                      ID = GO_results$ID,
                      Term = GO_results$Description,
                      Genes = str_replace_all(GO_results$geneID, "/", ","),
                      adj_pval = GO_results$p.adjust)

  # 创建一个基因数据框，包含基因ID和logFC等
  genes <- data.frame(ID = unlist(str_split(terms$Genes, ",")),
                      logFC = rnorm(length(unlist(str_split(terms$Genes, ",")))))  # 示例数据
  # 绘制弦图
  circ <- circle_dat(terms, genes)

  # 使用GOCircle函数绘制弦图
  GOCircle(circ)
  # 自定义颜色和基因数目限制
  GOCircle(circ, nsub = 10)  # 仅显示前10个GO term
}

#####################################################################
library(dplyr)

# Read the CSV file
data <- read.csv("linkList1_1.csv")

# Count the occurrences in the Source column
source_counts <- data %>%
  group_by(Source) %>%
  summarize(Source_Count = n()) %>%
  arrange(desc(Source_Count))

# Count the occurrences in the Target column
target_counts <- data %>%
  group_by(Target) %>%
  summarize(Target_Count = n())

# Merge the source and target counts, ensuring all Source genes are included
result <- left_join(source_counts, target_counts, by = c("Source" = "Target"))

# Replace NA values in Target_Count with 0 (for genes that don't appear as Target)
result <- result %>%
  mutate(Target_Count = ifelse(is.na(Target_Count), 0, Target_Count))

# View the final result
print(result)
# Load necessary libraries
library(ggplot2)

# Assuming 'result' contains the Source (outgoing) and Target (incoming) counts for all genes,
# Replace 'result' with your actual data
data <- data.frame(
  Gene = 1:nrow(result),  # X-axis: Sequential gene positions
  Outgoing = result$Source_Count,  # Y-axis (top): Outgoing connections
  Incoming = result$Target_Count   # Y-axis (bottom): Incoming connections
)

# Create the plot
ggplot(data, aes(x = Gene)) +
  # Outgoing points and lines (top section)
  geom_point(aes(y = Outgoing), color = "#F7776E", size = 2.5) +  # Increased point size
  geom_segment(aes(x = Gene, xend = Gene, y = 0, yend = Outgoing), color = "gray") +
  # Incoming points and lines (bottom section, inverted to keep everything positive)
  geom_point(aes(y = -Incoming), color = "#6C8FC6", size = 2.5) +  # Increased point size
  geom_segment(aes(x = Gene, xend = Gene, y = 0, yend = -Incoming), color = "gray") +
  # Customizing the axis to show positive numbers on both sides
  scale_y_continuous(labels = abs, breaks = seq(0, max(c(data$Outgoing, data$Incoming)), by = 25)) +
  # Adjust x-axis to remove the vertical line at x=0 by removing margin (expand=c(0,0))
  scale_x_continuous(breaks = seq(0, nrow(result), by = 50), expand = c(0, 0), limits = c(1, nrow(result) + 1)) +  # Adjust x-axis limits for padding
  # Customizing axes and plot appearance
  theme_minimal(base_size = 15) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 2),  # Black frame around plot
    axis.line = element_line(size = 1.2, color = "gray30"),  # Lighter black axis lines
    panel.grid.major = element_line(color = "lightgray", size = 0.8),  # Thicker grid lines
    panel.grid.minor = element_blank(),  # No minor grid lines
    axis.text = element_text(size = 12),  # Adjust text size for clarity
    plot.margin = margin(5, 15, 5, 5)  # Adjust right margin to ensure full display
  ) +
  labs(x = "Modules", y = "Counts")
###################################################################################

library(dplyr)

# Read the CSV file
data <- read.csv("linkList2_1.csv")

# Count the occurrences in the Source column
source_counts <- data %>%
  group_by(Source) %>%
  summarize(Source_Count = n()) %>%
  arrange(desc(Source_Count))

# Count the occurrences in the Target column
target_counts <- data %>%
  group_by(Target) %>%
  summarize(Target_Count = n())

# Merge the source and target counts, ensuring all Source genes are included
result <- left_join(source_counts, target_counts, by = c("Source" = "Target"))

# Replace NA values in Target_Count with 0 (for genes that don't appear as Target)
result <- result %>%
  mutate(Target_Count = ifelse(is.na(Target_Count), 0, Target_Count))

# View the final result
print(result)
# Load necessary libraries
library(ggplot2)

# Assuming 'result' contains the Source (outgoing) and Target (incoming) counts for all genes,
# Replace 'result' with your actual data
data <- data.frame(
  Gene = 1:nrow(result),  # X-axis: Sequential gene positions
  Outgoing = result$Source_Count,  # Y-axis (top): Outgoing connections
  Incoming = result$Target_Count   # Y-axis (bottom): Incoming connections
)

# Create the plot
ggplot(data, aes(x = Gene)) +
  # Outgoing points and lines (top section)
  geom_point(aes(y = Outgoing), color = "#1DBDC6", size = 2.5) +  # Increased point size
  geom_segment(aes(x = Gene, xend = Gene, y = 0, yend = Outgoing), color = "gray") +
  # Incoming points and lines (bottom section, inverted to keep everything positive)
  geom_point(aes(y = -Incoming), color = "#F6C490", size = 2.5) +  # Increased point size
  geom_segment(aes(x = Gene, xend = Gene, y = 0, yend = -Incoming), color = "gray") +
  # Customizing the axis to show positive numbers on both sides
  scale_y_continuous(labels = abs, breaks = seq(0, max(c(data$Outgoing, data$Incoming)), by = 25)) +
  # Adjust x-axis to remove the vertical line at x=0 by removing margin (expand=c(0,0))
  scale_x_continuous(breaks = seq(0, nrow(result), by = 50), expand = c(0, 0), limits = c(1, nrow(result) + 1)) +  # Adjust x-axis limits for padding
  # Customizing axes and plot appearance
  theme_minimal(base_size = 15) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 2),  # Black frame around plot
    axis.line = element_line(size = 1.2, color = "gray30"),  # Lighter black axis lines
    panel.grid.major = element_line(color = "lightgray", size = 0.8),  # Thicker grid lines
    panel.grid.minor = element_blank(),  # No minor grid lines
    axis.text = element_text(size = 12),  # Adjust text size for clarity
    plot.margin = margin(5, 15, 5, 5)  # Adjust right margin to ensure full display
  ) +
  labs(x = "Modules", y = "Counts")
##############################################################################

library(dplyr)

# Read the CSV file
data <- read.csv("linkList2_3.csv")

# Count the occurrences in the Source column
source_counts <- data %>%
  group_by(Source) %>%
  summarize(Source_Count = n()) %>%
  arrange(desc(Source_Count))

# Count the occurrences in the Target column
target_counts <- data %>%
  group_by(Target) %>%
  summarize(Target_Count = n())

# Merge the source and target counts, ensuring all Source genes are included
result <- left_join(source_counts, target_counts, by = c("Source" = "Target"))

# Replace NA values in Target_Count with 0 (for genes that don't appear as Target)
result <- result %>%
  mutate(Target_Count = ifelse(is.na(Target_Count), 0, Target_Count))

# View the final result
print(result)
# Load necessary libraries
library(ggplot2)

# Assuming 'result' contains the Source (outgoing) and Target (incoming) counts for all genes,
# Replace 'result' with your actual data
data <- data.frame(
  Gene = 1:nrow(result),  # X-axis: Sequential gene positions
  Outgoing = result$Source_Count,  # Y-axis (top): Outgoing connections
  Incoming = result$Target_Count   # Y-axis (bottom): Incoming connections
)

# Create the plot
ggplot(data, aes(x = Gene)) +
  # Outgoing points and lines (top section)
  geom_point(aes(y = Outgoing), color = "#F6C6CB", size = 2.5) +  # Increased point size
  geom_segment(aes(x = Gene, xend = Gene, y = 0, yend = Outgoing), color = "gray") +
  # Incoming points and lines (bottom section, inverted to keep everything positive)
  geom_point(aes(y = -Incoming), color = "#6CB4AE", size = 2.5) +  # Increased point size
  geom_segment(aes(x = Gene, xend = Gene, y = 0, yend = -Incoming), color = "gray") +
  # Customizing the axis to show positive numbers on both sides
  scale_y_continuous(labels = abs, breaks = seq(0, max(c(data$Outgoing, data$Incoming)), by = 25)) +
  # Adjust x-axis to remove the vertical line at x=0 by removing margin (expand=c(0,0))
  scale_x_continuous(breaks = seq(0, nrow(result), by = 50), expand = c(0, 0), limits = c(1, nrow(result) + 1)) +  # Adjust x-axis limits for padding
  # Customizing axes and plot appearance
  theme_minimal(base_size = 15) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 2),  # Black frame around plot
    axis.line = element_line(size = 1.2, color = "gray30"),  # Lighter black axis lines
    panel.grid.major = element_line(color = "lightgray", size = 0.8),  # Thicker grid lines
    panel.grid.minor = element_blank(),  # No minor grid lines
    axis.text = element_text(size = 12),  # Adjust text size for clarity
    plot.margin = margin(5, 15, 5, 5)  # Adjust right margin to ensure full display
  ) +
  labs(x = "Modules", y = "Counts")
################################################################################
library(dplyr)

# Read the CSV file
data <- read.csv("linkList3_1.csv")

# Count the occurrences in the Source column
source_counts <- data %>%
  group_by(Source) %>%
  summarize(Source_Count = n()) %>%
  arrange(desc(Source_Count))

# Count the occurrences in the Target column
target_counts <- data %>%
  group_by(Target) %>%
  summarize(Target_Count = n())

# Merge the source and target counts, ensuring all Source genes are included
result <- left_join(source_counts, target_counts, by = c("Source" = "Target"))

# Replace NA values in Target_Count with 0 (for genes that don't appear as Target)
result <- result %>%
  mutate(Target_Count = ifelse(is.na(Target_Count), 0, Target_Count))

# View the final result
print(result)
# Load necessary libraries
library(ggplot2)

# Assuming 'result' contains the Source (outgoing) and Target (incoming) counts for all genes,
# Replace 'result' with your actual data
data <- data.frame(
  Gene = 1:nrow(result),  # X-axis: Sequential gene positions
  Outgoing = result$Source_Count,  # Y-axis (top): Outgoing connections
  Incoming = result$Target_Count   # Y-axis (bottom): Incoming connections
)

# Create the plot
ggplot(data, aes(x = Gene)) +
  # Outgoing points and lines (top section)
  geom_point(aes(y = Outgoing), color = "#FDD400", size = 2.5) +  # Increased point size
  geom_segment(aes(x = Gene, xend = Gene, y = 0, yend = Outgoing), color = "gray") +
  # Incoming points and lines (bottom section, inverted to keep everything positive)
  geom_point(aes(y = -Incoming), color = "#8A3684", size = 2.5) +  # Increased point size
  geom_segment(aes(x = Gene, xend = Gene, y = 0, yend = -Incoming), color = "gray") +
  # Customizing the axis to show positive numbers on both sides
  scale_y_continuous(labels = abs, breaks = seq(0, max(c(data$Outgoing, data$Incoming)), by = 25)) +
  # Adjust x-axis to remove the vertical line at x=0 by removing margin (expand=c(0,0))
  scale_x_continuous(breaks = seq(0, nrow(result), by = 50), expand = c(0, 0), limits = c(1, nrow(result) + 1)) +  # Adjust x-axis limits for padding
  # Customizing axes and plot appearance
  theme_minimal(base_size = 15) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 2),  # Black frame around plot
    axis.line = element_line(size = 1.2, color = "gray30"),  # Lighter black axis lines
    panel.grid.major = element_line(color = "lightgray", size = 0.8),  # Thicker grid lines
    panel.grid.minor = element_blank(),  # No minor grid lines
    axis.text = element_text(size = 12),  # Adjust text size for clarity
    plot.margin = margin(5, 15, 5, 5)  # Adjust right margin to ensure full display
  ) +
  labs(x = "Modules", y = "Counts")
#######################################################
# 加载所需包
library(ggplot2)
library(dplyr)

# 创建数据框
gene_data <- data.frame(
  Group = c("Fetal_10", "Fetal_17", "hHep", "hLivEC"),
  Total = c(163, 182, 146, 76),
  Annotated = c(148, 167, 136, 75)
)

# 计算未注释基因数
gene_data <- gene_data %>%
  mutate(Unannotated = Total - Annotated)

# 将数据转换为长格式以便ggplot使用
gene_data_long <- gene_data %>%
  select(Group, Annotated, Unannotated) %>%
  tidyr::pivot_longer(cols = c("Annotated", "Unannotated"),
                      names_to = "Status", values_to = "Count")

# 定义颜色方案
fill_colors <- c("Fetal_10.Annotated" = "#E64825", "Fetal_10.Unannotated" = "#F6C0CC",
                 "Fetal_17.Annotated" = "#715EA9", "Fetal_17.Unannotated" = "#CAC0E1",
                 "hHep.Annotated" = "#6A9ACE", "hHep.Unannotated" = "#ABDAEC",
                 "hLivEC.Annotated" = "#1E803D", "hLivEC.Unannotated" = "#97D1A0")

# 绘制堆叠条形图
ggplot(gene_data_long, aes(x = Group, y = Count, fill = interaction(Group, Status))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = fill_colors) +
  labs(title = "Proportion of Annotated and Unannotated Genes",
       x = "Gene Groups", y = "Gene Count",) +
  theme_minimal() +
  theme(plot.caption = element_text(hjust = 0.5))  # 居中显示题注
