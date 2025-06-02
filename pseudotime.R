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

