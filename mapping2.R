library(ggsci)
library(ggpubr)
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080")

q1=plot_cell_trajectory(cds, color_by = "cell_type") +  scale_color_manual(values=colour) + 
  theme(
    plot.title = element_text(size = 20),       # 标题字体大小
    axis.title = element_text(size = 23),       # 轴标签字体大小
    axis.text = element_text(size = 20),        # 轴刻度字体大小
    legend.title = element_text(size = 16),     # 图例标题字体大小
    legend.text = element_text(size = 23)       # 图例文字字体大小
  )
q2=plot_cell_trajectory(cds, color_by = "State") + 
  scale_color_manual(values=colour) + 
  theme(
    plot.title = element_text(size = 20),       # 标题字体大小
    axis.title = element_text(size = 23),       # 轴标签字体大小
    axis.text = element_text(size = 20),        # 轴刻度字体大小
    legend.title = element_text(size = 16),     # 图例标题字体大小
    legend.text = element_text(size = 23)       # 图例文字字体大小
  )
arrows <- data.frame(
  x = c(8, 2, 5),   # 箭头起点X坐标
  y = c(4, -1, -2),     # 箭头起点Y坐标
  xend = c(3, -8, 6), # 箭头终点X坐标
  yend = c(0, 1, -7)   # 箭头终点Y坐标
)
q1=q1 + geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend),
                     arrow = arrow(type = "open", length = unit(0.2, "inches")),
                     color = "black",  # 设置为浅灰色
                     alpha = 5,       # 增加透明度
                     size = 1)        # 设置线的粗细
q2=q2 + geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend),
                     arrow = arrow(type = "open", length = unit(0.2, "inches")),
                     color = "black",  # 设置为浅灰色
                     alpha = 5,       # 增加透明度
                     size = 1)        # 设置线的粗细
q1|q2
# ggsave("拟时序图.pdf", q1|q2, width = 246, height = 175, units = "mm",dpi=300)

plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)
q3<-plot_complex_cell_trajectory(cds,x=1,y=2,
                                 color_by="cell_type")+
  scale_color_manual(values=colour)+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 23))       # 图例文字字体大小
print(q3)
# ggsave("树图.pdf", q3, width = 174, height = 144, units = "mm",dpi=300)

q4<-ggplot(pData(cds),aes(Pseudotime,colour=State,fill=State))+geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()
print(q4)
qq4<-ggplot(pData(cds),aes(Pseudotime,colour=cell_type,fill=cell_type))+geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()
print(qq4)

s.genes <- c("HBM","ORM1","AFP","CAT","ADH1B","TUBA1B","MT1G","FCN3")
q5 <- plot_genes_jitter(cds[s.genes,], grouping = "State", color_by = "State")+ 
  theme(
    plot.title = element_text(size = 24),
    axis.text = element_text(size = 23),        # 轴刻度字体大小
    legend.title = element_text(size = 16),     # 图例标题字体大小
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    legend.text = element_text(size = 23),       # 图例文字字体大小,       
    strip.text = element_text(size = 24)        # 小图标题字体大小
  )
q6 <- plot_genes_violin(cds[s.genes,], grouping = "State", color_by = "State")+ 
  theme(
    plot.title = element_text(size = 24),
    axis.text = element_text(size = 23),        # 轴刻度字体大小
    legend.title = element_text(size = 16),     # 图例标题字体大小
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    legend.text = element_text(size = 23),       # 图例文字字体大小
    strip.text = element_text(size = 24)        # 小图标题字体大小
  )
q7 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "State")+ 
  theme(
    plot.title = element_text(size = 24),
    axis.text = element_text(size = 23),        # 轴刻度字体大小
    legend.title = element_text(size = 16),     # 图例标题字体大小
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    legend.text = element_text(size = 23),       # 图例文字字体大小
    strip.text = element_text(size = 24)        # 小图标题字体大小
  )
plotc <- q5|q6|q7
print(plotc)
# ggsave("小提琴图.pdf", plotc, width = 381, height = 481, units = "mm",dpi=300)

# 加载所需的库
library(ggplot2)

# 创建示例数据
data <- data.frame(
  Category = c("1-1", "1-1", "2-3", "2-3"),
  Direction = c("Out", "In", "Out", "In"),
  Count = c(51, -12, 53, -8)
)

# 绘制双向柱状图
q8 <- ggplot(data, aes(x = Category, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_manual(values = c("Out" = "lightcoral", "In" = "lightblue"), name = "Direction") +
  labs(title = "Counts of 'Out' and 'In' by Gene regulation Networks", x = "Category", y = "Counts") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  geom_text(aes(label = abs(Count)), position = position_stack(vjust = 0.5), color = "black")

# 显示图形
print(q8)

