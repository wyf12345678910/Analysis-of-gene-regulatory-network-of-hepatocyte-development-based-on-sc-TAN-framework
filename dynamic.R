
#############################################################################
LHrmP <- expression(
  1,
  2*tt,
  (4 * tt^2 - 2),
  (8 * tt^3 - 12 * tt),
  (16 * tt^4 - 48 * tt^2 + 12),
  (32 * tt^5 - 160 * tt^3 + 120 * tt)
  # 继续定义更高次的赫米特多项式...
)

# 修改 Legendre() 函数以使用赫米特多项式
Legendre <- function(par, r, times) {
  tnum <- length(times)
  f <- c()
  for (t in 1:tnum) {
    tt <- -1 + 2 * (times[t] - times[1]) / (times[tnum] - times[1])
    ff <- 0
    for (i in 1:r) {
      ff <- ff + par[i] * eval(LHrmP[i])
    }
    f <- c(f, ff)
  }
  return(f)
}

s.mle<-function(par,data,r,times){
  y <- Legendre(par,r,times)
  yi <- data
  res <- sum((yi-y)^2)
  return(res)
}#最小二乘

get_Legendre_par<-function(initial_f_par,data,r,times){
  
  a <- optim(initial_f_par,s.mle,data=data,r=r,times=times,method = "Nelder-Mead",control=list(maxit=2000,trace=FALSE))
  #optim函数来执行最小化操作，以找到Legendre多项式拟合的最佳参数
  curve_par_i<-a$par
  return(list(par=a$par,res=a$value))
}#对给定的输入数据进行Legendre多项式拟合，并返回拟合结果

detect.dynamic<-function(x){
  
  peak<-c()
  trough<-c()
  diffx<-diff(x, lag=1)
  signx<-sign(diffx)
  for(ii in 1:(length(signx)-1)){
    if(signx[ii]>signx[ii+1])
      peak<-c(peak,ii)
    if(signx[ii]<signx[ii+1]){
      trough<-c(trough,ii)
    }
  }
  if(length(signx>0)==length(signx)){
    result1<-c(-1,length(x))
  }
  if(length(signx<0)==length(signx)){
    result1<-c(1,-length(x))
  }
  if((length(peak)==1)&(length(trough)==0)){
    result1<-c(-1,peak,-length(x))
  }
  if((length(peak)==0)&(length(trough)==1)){
    result1<-c(1,-trough,length(x))
  }
  if((length(peak)==1)&(length(trough)==1)){
    if(peak>trough)
      result1<-c(1,-trough,peak,-length(x))
    if(peak<trough)
      result1<-c(-1,peak,-trough,length(x))
  }
  if((length(peak)==1)&(length(trough)==2)){
    result1<-c(1,-trough[1],peak,-trough[2],length(x))
  }
  if((length(peak)==2)&(length(trough)==1)){
    result1<-c(-1,peak[1],-trough,peak[2],-length(x))
  }
  
  return(result=result1)
}


detect.patterns<-function(dy,fit_result2,type,type.name,r,times,plot1=TRUE,plot2=TRUE){
  
  dytype<-c()
  basic.patterns<-c()
  for(i in 1:length(dy)){
    type1<-type[abs(dy[[i]])]
    type1[which(dy[[i]]<0)]<-paste("low",type1[which(dy[[i]]<0)],sep="")
    dytype<-c(dytype,paste(type1,collapse="-"))
    
    flu1<-c("up","down","up","down","up","down")
    flu2<-c("down","up","down","up","down","up")
    if(dy[[i]][1]<0)
      basic.patterns1<-paste(flu1[1:length(dy[[i]])-1],collapse ="-")
    if(dy[[i]][1]>0)
      basic.patterns1<-paste(flu2[1:length(dy[[i]])-1],collapse ="-")
    basic.patterns<-c(basic.patterns,basic.patterns1)
  }
  
  #table(dytype)
  if(plot1==TRUE){
    for(i in 1:length(table(dytype))){
      tf_type<-which(dytype==names(table(dytype))[i])
      g1<-same.patterns.plot(fit_result2[,tf_type],r,times,type.name,type,genename=colnames(fit_result2)[tf_type])
      
      png(filename=paste(names(table(dytype))[i],".png",sep=""),width=600,height=500)
      multiplot (g1,cols=1)
      dev.off()
    }
  }
  
  pvalue<-c()
  dynamic.pattern<-c()
  for(i in 1:length(dy)){
    x<-fit_result2[,i]
    peaks<-dy[[i]][which(dy[[i]]>0)]
    maxvalue<-peaks[which.max(x[peaks])]
    if(type[maxvalue]==type.name[1]&maxvalue<=length(which(type==type[maxvalue]))/2){
      pvalue<-c(pvalue,NA)
      dynamic.pattern<-c(dynamic.pattern,type[maxvalue])
    }
    if(type[maxvalue]==type.name[1]&maxvalue>length(which(type==type[maxvalue]))/2){
      e1<-x[c(1:length(which(type==type[maxvalue])))]
      e2<-x[c((length(which(type==type[maxvalue]))+1):length(c(which(type==type.name[1]),which(type==type.name[2]))))]
      wilcox_test<-wilcox.test(e1, e2)
      pvalue1<-wilcox_test$p.value
      pvalue<-c(pvalue,pvalue1)
      if(pvalue1>0.01)
        dynamic.pattern<-c(dynamic.pattern,paste(type[maxvalue],type.name[2],sep="-"))
      if(pvalue1<0.01)
        dynamic.pattern<-c(dynamic.pattern,type[maxvalue])
    }
    if(type[maxvalue]==rev(type.name)[1]&maxvalue<length(which(type==type[maxvalue]))/2){
      e1<-x[c((length(which(type!=type[maxvalue]))+1):length(type))]
      e2<-x[c(length(c(which(type==type.name[length(table(type))]),which(type==type.name[length(table(type))-1]))):length(which(type!=type[maxvalue])))]
      wilcox_test<-wilcox.test(e1, e2) 
      pvalue1<-wilcox_test$p.value
      pvalue<-c(pvalue,pvalue1)
      if(pvalue1>0.01)
        dynamic.pattern<-c(dynamic.pattern,paste(type.name[length(table(type))-1],type[maxvalue],sep="-"))
      if(pvalue1<0.01)
        dynamic.pattern<-c(dynamic.pattern,type[maxvalue])
    }
    if(type[maxvalue]==rev(type.name)[1]&maxvalue>=length(which(type==type[maxvalue]))/2){
      pvalue<-c(pvalue,NA)
      dynamic.pattern<-c(dynamic.pattern,type[maxvalue])
    }
    if(type[maxvalue]!=rev(type.name)[1]&type[maxvalue]!=type.name[1]&maxvalue<length(which(type==type[maxvalue]))/2){
      e1<-x[which(type==type[maxvalue])]
      e2<-x[which(type==type.name[which(type.name==type[maxvalue])-1])]
      wilcox_test<-wilcox.test(e1, e2) 
      pvalue1<-wilcox_test$p.value
      pvalue<-c(pvalue,pvalue1)
      if(pvalue1>0.01)
        dynamic.pattern<-c(dynamic.pattern,paste(type.name[which(type.name==type[maxvalue])-1],type[maxvalue],sep="-"))
      if(pvalue1<0.01)
        dynamic.pattern<-c(dynamic.pattern,type[maxvalue])
    }
    if(type[maxvalue]!=rev(type.name)[1]&type[maxvalue]!=type.name[1]&maxvalue>=length(which(type==type[maxvalue]))/2){
      e1<-x[which(type==type[maxvalue])]
      e2<-x[which(type==type.name[which(type.name==type[maxvalue])+1])]
      wilcox_test<-wilcox.test(e1, e2) 
      pvalue1<-wilcox_test$p.value
      pvalue<-c(pvalue,pvalue1)
      if(pvalue1>0.01)
        dynamic.pattern<-c(dynamic.pattern,paste(type[maxvalue],type.name[which(type.name==type[maxvalue])+1],sep="-"))
      if(pvalue1<0.01)
        dynamic.pattern<-c(dynamic.pattern,type[maxvalue])
    }
  }
  
  if(plot2==TRUE){
    for(i in 1:length(table(dynamic.pattern))){
      tf_type<-which(dynamic.pattern==names(table(dynamic.pattern))[i])
      g1<-same.patterns.plot(fit_result2[,tf_type],r,times,type.name,type,genename=colnames(fit_result2[,tf_type]))
      
      png(filename=paste(names(table(dynamic.pattern))[i],".2",".png",sep=""),width=600,height=500)
      multiplot (g1,cols=1)
      dev.off()
    }
  }
  result_list<-data.frame(basic.patterns,pattern1=dytype,patterns2=dynamic.pattern,p.value=pvalue)
  rownames(result_list)<-colnames(fit_result2)
  return(result_list) 
}



same.patterns.plot<-function(fit,r,times,type.name,type,genename){
  fit_result<-as.data.frame(cbind(times,fit))
  timepoint1<-0
  for(ii in 1:length(type.name)){
    timepoint1<-c(timepoint1,length(which(type==type.name[ii])))
  }
  timepoint<-cumsum(timepoint1)
  colorlist<-brewer.pal(length(type.name),"Set3")
  g1 <- ggplot(fit_result)
  for(ii in 1:length(type.name)){
    ndata<-cbind(fit_result,tt=c(timepoint[ii],timepoint[ii+1],rep(0,length(type)-2)))
    g1 <- g1 + geom_rect(data=ndata,aes(xmin=tt[1], xmax=tt[2], ymin=-2, ymax=16), fill=colorlist[ii], color=NA, alpha=0.1)
    g1 <- g1 + annotate("text", x = (ndata$tt[1]+ndata$tt[2])/2, y =14, label =type.name[ii],family="serif",size=6)
  }
  if(dim(fit_result)[2]==2){
    p.data<-data.frame(time=times,fit=fit_result[,2])
    g1 <- g1 + geom_line(data=p.data,aes(x=time,y=fit),size=0.6,alpha=0.8)
    g1 <- g1 + annotate("text", x = length(times), y =rev(p.data$fit)[1], label =genename,family="serif",size=5)
  }
  if(dim(fit_result)[2]>2){
    for(j in 1:dim(fit)[2]){
      p.data<-data.frame(time=times,fit=fit_result[,j+1])
      #g1 <- g1 + geom_point(data=p.data,aes(x=time,y=true),colour=color_list[j],size=1,pch=1)
      g1 <- g1 + geom_line(data=p.data,aes(x=time,y=fit),size=0.6,alpha=0.8)
      g1 <- g1 + annotate("text", x = 0, y =p.data$fit[1], label =genename[j],family="serif",size=5)
    }
  }
  
  
  g1 <- g1 + scale_y_continuous(limits=c(-2,16),breaks=seq(-2,16,2),labels=seq(-2,16,2))
  g1 <- g1 + scale_x_continuous(limits=c(0,max(times)+2))
  g1 <- g1 + xlab("pseudotime")+ylab("log2(FPKM)") 
  
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}





#############################################TF拟合###############################
exp_new<-t(exp_new)
head(TFlist)
dim(datExpr)
table(type)
# type_news<-type[c(which(type=="fetal_10"),which(type=="fetal_17"),which(type=="hHep"),which(type=="hLivEC"))]
tfMatr<-exp_new[intersect(TFlist,rownames(exp_new)),]
#从datExpr中的基因表达矩阵中选择与TFlist中第二列相交的基因，并将结果保存在tfMatr中。
tfMatr<-tfMatr[,c(which(type_new=="fetal_10"),which(type_new=="fetal_17"),which(type_new=="hHep"),which(type_new=="hLivEC"))]
#为tfMatr选择与指定分组条件对应的列，并只保留这些列中的数据。
dim(tfMatr)
# type<-type_news

length(toptf_list)
r=5#将r的值设置为5，用作后续拟合过程中的参数
times=seq(1,dim(tfMatr)[2],1)
fit_result1<-c()#存储每次拟合的结果
for(j in 1:dim(tfMatr)[1]){
  data<-tfMatr[j,]
  result<-get_Legendre_par(initial_f_par=rep(0.01,r),data,r,times)#该函数为自编在下方
  fit_result1<-cbind(fit_result1,Legendre(result$par,r,times))
}
colnames(fit_result1)<-rownames(tfMatr)
dim(fit_result1)
type.name=c("fetal_10","fetal_17","hHep","hLivEC")
set.seed(1234)
cr<-kmeans(t(fit_result1),centers=5)
#对矩阵fit_result1进行k均值聚类分析，并将结果保存在变量cr中。分为5个聚类中心。
clus<-cr$cluster
#将聚类结果中的类别标签保存在变量clus中
library(NbClust)#用于评估聚类的最佳数量
nc<-NbClust(data=t(fit_result1),method ="kmeans",min.nc=2,max.nc=10,index="ch")
#使用NbClust函数对矩阵fit_result1进行聚类分析，尝试不同的聚类数目，使用"ch"指数评估聚类效果
nc$Best.nc
#输出评估结果中的最佳聚类数量
clus<-nc$Best.partition
#将最佳聚类数量的分区结果保存在变量clus中
library(mclust)#用于进行混合高斯模型聚类分析
cr<-Mclust(t(fit_result1))
#对矩阵fit_result1进行混合高斯模型聚类分析，并将结果保存在变量cr中
clus<-cr$classification
#将混合高斯模型聚类结果中的类别标签保存在变量clus中
g<-list()
#for(m in 1:length(table(clus))){
for(m in 7:9){
  times <- c(1:dim(tfMatr)[2])
  fit_result<-fit_result1[,which(clus==m)]
  true_value<- t(tfMatr[which(clus==m),])
  g[[m]]<-same.patterns.plot(fit_result1,r,times,type.name,type_new,genename=colnames(fit_result1))
}
cowplot::plot_grid(g[[1]],g[[2]],g[[3]],g[[4]],g[[5]],nrow=2)
cowplot::plot_grid(g[[7]],g[[8]],g[[9]],nrow=1)
g[[1]]


png("dynamic.png",width=2000,height=1500)
g[[9]]
dev.off()


x<-fit_result1[,clus2][,1]
plot(x=times,y=x,type="l")

g<-same.patterns.plot(fit_result1,r,times,type.name,type_new,genename=colnames(fit_result1))
g

cutval<-0.3
varval<-apply(tfMatr,MARGIN=1,var)
clus1<-which(varval<cutval)
clus2<-which(varval>=cutval)
#intersect(rownames(tfMatr)[clus1],c("Tcf7l2","Nfib","Id3","Sox11","Ascl1","Myt1","Tub","Insm1"))
g1<-same.patterns.plot(fit_result1[,clus1],r,times,type.name,type_new,genename=colnames(fit_result1[,clus1]))
g1
g2<-same.patterns.plot(fit_result1[,clus2],r,times,type.name,type_new,genename=colnames(fit_result1[,clus2]))
g2

#detect peak and trough
fit_result2<-fit_result1
#fit_result2<-fit_result1[,clus2]
tfMatr2<-t(tfMatr)
#tfMatr2<-t(tfMatr)[,clus2]
dy<-apply(fit_result2,MARGIN=2,detect.dynamic)

table(type)
plot(x=times,y=x)

pattern<-detect.patterns(dy,fit_result2,type_new,type.name,r,times,plot1=TRUE,plot2=TRUE)
