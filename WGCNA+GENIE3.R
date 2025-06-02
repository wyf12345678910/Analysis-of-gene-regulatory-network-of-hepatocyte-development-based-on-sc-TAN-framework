###########################WGCNA detect module##########################################
exp_new<-t(exp_new)
powers = c(c(1:30))
# Call the network topology analysis function
sft = pickSoftThreshold(exp_new, powerVector = powers, verbose = 5,
                        networkType = "signed",RsquaredCut = 0.8)
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
  power = 27,
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
ggsave("热图.pdf", pheat, width = 175, height = 144, units = "mm",dpi=300)
save(toptf_list,file="toptf_list.Rdata")
load("toptf_list.Rdata")
