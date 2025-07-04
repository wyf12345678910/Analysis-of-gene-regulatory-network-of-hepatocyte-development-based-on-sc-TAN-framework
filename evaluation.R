library(ROCR)
library(infotheo)
library(ggplot2)


###our model
#predition network
linklist1<- readRDS("./hliver_result/SCENIC/allstep/int/type1/1.4_GENIE3_linkList.Rds")
linklist2<- readRDS("./hliver_result/SCENIC/allstep/int/type2/1.4_GENIE3_linkList.Rds")
linklist3<- readRDS("./hliver_result/SCENIC/allstep/int/type3/1.4_GENIE3_linkList.Rds")
linklist4<- readRDS("./hliver_result/SCENIC/allstep/int/type4/1.4_GENIE3_linkList.Rds")
linklist5<- readRDS("./hliver_result/SCENIC/allstep/int/type5-1/1.4_GENIE3_linkList.Rds")
linklist6<- readRDS("./hliver_result/SCENIC/allstep/int/type5-2/1.4_GENIE3_linkList.Rds")
linklist<-rbind(linklist1,linklist2,linklist3,linklist4,linklist5,linklist6)
dim(linklist)



###mef
tf_network1<-read.csv("D:/bjfu/paper/scRNA/Our Model/MEF/tf_network.csv")
dim(tf_network1)
tf_network<-tf_network1[tf_network1[,2]%in%rownames(exp_mef),]




##Our model
load("./Our model/GENIE3_on_net1_AnimalDB/allinklist.Rdata")
load("./Our model/GENIE3_on_net1_2DB/allinklist.Rdata")
load("./Our model/allinklist.Rdata")
dim(allinklist)
linklist<-data.frame(regulatoryGene=str_to_title(allinklist[,1]),
                     targetGene=str_to_title(allinklist[,2]),
                     weight=allinklist[,3])

##SCENIC
#link_table<-readRDS("./hliver_result/SCENIC/onlySCENIC/int/2.5_regulonTargetsInfo.Rds")
link_table<-readRDS("./SCENIC/int/2.5_regulonTargetsInfo.Rds")
linklist<-data.frame(TF=unlist(link_table[,1]),Target=unlist(link_table[,2]),weight=as.numeric(unlist(link_table[,10])))
dim(linklist)


##SCODE
linklist<-read.csv("./SCODE/cyto_edge.txt",sep="\t")
dim(linklist)


##LEAP
load("./LEAP//MAC_results.Rdata")
dim(MAC_results)
head(MAC_results)
load("./LEAP/data_leap.Rdata")
linklist<-data.frame(TF=rownames(data)[MAC_results[,3]],
                     Target=rownames(data)[MAC_results[,4]],
                     weight=as.numeric(MAC_results[,1]))
dim(linklist)


AUC(linklist,tf_network)
ARI(linklist,tf_network)
NMI(linklist,tf_network)



AUC<-function(linklist,tf_network){
  coname<-paste(linklist[,1],linklist[,2],sep="-")
  coname1<-paste(tf_network[,1],tf_network[,2],sep="-")
  link<-rep(0,dim(linklist)[1])
  link[which(coname%in%coname1)]<-1
  Net<-data.frame(TF=linklist[,1],Target=linklist[,2],weight=linklist[,3],link=as.factor(link))
  if(length(which(is.na(Net[,3])))!=0)
    Net<-Net[-which(is.na(Net[,3])),]
  pred<-prediction(Net$weight,Net$link)
  
  
  #AUC
  perf<-performance(pred,measure="tpr",x.measure="fpr")
  #plot(perf,colorize=T)
  
  f<-approxfun(x=unlist(perf@x.values),y=unlist(perf@y.values))
  auc<-integrate(f,lower = 0 , upper = 1 , subdivisions = 1000)         
  return(auc)
}


PR<-function(linklist,tf_network){
  coname<-paste(linklist[,1],linklist[,2],sep="-")
  coname1<-paste(tf_network[,1],tf_network[,2],sep="-")
  link<-rep(0,dim(linklist)[1])
  link[which(coname%in%coname1)]<-1
  Net<-data.frame(TF=linklist[,1],Target=linklist[,2],weight=linklist[,3],link=as.factor(link))
  if(length(which(is.na(Net[,3])))!=0)
    Net<-Net[-which(is.na(Net[,3])),]
  pred<-prediction(Net$weight,Net$link)
  
  perf<-performance(pred,measure="rec",x.measure="prec")
  #plot(perf,colorize=T)
  
  f<-approxfun(x=unlist(perf@x.values)[-1],y=unlist(perf@y.values)[-1])
  pr<-integrate(f,lower = min(unlist(perf@x.values)[-1]) , upper = max(unlist(perf@x.values)[-1]) , subdivisions = 1000)         
  return(pr)
}





ARI<-function(linklist,tf_network){
  coname<-paste(linklist[,1],linklist[,2],sep="-")
  coname1<-paste(tf_network[,1],tf_network[,2],sep="-")
  n11<-length(which(coname%in%coname1))
  n01<-dim(linklist)[1]-length(which(coname%in%coname1))
  n10<-length(coname1)-length(which(coname1%in%coname))
  n00<-0
  N<-matrix(c(n11,n01,n10,n00),ncol=2)
  n<-sum(N)
  c1<-n11+n10
  c0<-n01+n00
  c1t<-n11+n01
  c0t<-n10+n00
  t1<-c1*(c1-1)/2+c0*(c0-1)/2
  t2<-c1t*(c1t-1)/2+c0t*(c0t-1)/2
  t3<-2*t1*t2/(n*(n-1))
  ari<-((n11*(n11-1)/2+n01*(n01-1)/2+n10*(n10-1)/2)-t3)/((t1+t2)/2-t3)
  return(ari)
}


NMI<-function(linklist,tf_network){
  coname<-paste(linklist[,1],linklist[,2],sep="-")
  coname1<-paste(tf_network[,1],tf_network[,2],sep="-")
  link<-rep(0,dim(linklist)[1])
  link[which(coname%in%coname1)]<-1
  Net1<-data.frame(TF=linklist[,1],Target=linklist[,2],cp=rep(1,dim(linklist)[1]),ct=as.factor(link))
  remain<-dim(tf_network)[1]-length(which(coname1%in%coname))
  Net2<-data.frame(TF=tf_network[-c(which(coname1%in%coname)),1],
                   Target=tf_network[-c(which(coname1%in%coname)),2],
                   cp=rep(0,remain),
                   ct=as.factor(rep(1,remain)))
  Net<-rbind(Net1,Net2)
  cp<-Net[,3]
  ct<-Net[,4]
  I<-mutinformation(ct,cp)
  nmi<-I/sqrt(entropy(ct)*entropy(cp)) 
  return(nmi)
}


#AUC1待参考
AUC1<-function(linklist,tf_network){
  coname<-paste(linklist[,1],linklist[,2],sep="-")
  coname1<-paste(tf_network[,1],tf_network[,2],sep="-")
  link<-rep(0,dim(linklist)[1])
  link[which(coname%in%coname1)]<-1
  Net1<-data.frame(TF=linklist[,1],Target=linklist[,2],weight=linklist[,3],link=as.factor(link))
  if(length(which(is.na(Net1[,3])))!=0)
    Net1<-Net1[-which(is.na(Net1[,3])),]
  remain<-dim(tf_network)[1]-length(which(coname1%in%coname))
  Net2<-data.frame(TF=tf_network[-c(which(coname1%in%coname)),1],
                   Target=tf_network[-c(which(coname1%in%coname)),2],
                   weight=rep(0,remain),
                   link=as.factor(rep(1,remain)))
  Net<-rbind(Net1,Net2)
  
  pred<-prediction(Net$weight,Net$link)
  perf<-performance(pred,measure="tpr",x.measure="fpr")
  #plot(perf,colorize=T)
  
  f<-approxfun(x=unlist(perf@x.values),y=unlist(perf@y.values))
  auc<-integrate(f,lower = 0 , upper = 1 , subdivisions = 1000)         
  return(auc)
}


ARI1<-function(linklist,tf_network){
  coname<-paste(linklist[,1],linklist[,2],sep="-")
  coname1<-paste(tf_network[,1],tf_network[,2],sep="-")
  n<-dim(linklist)[1]
  n11<-length(which(coname%in%coname1))
  n01<-dim(linklist)[1]-length(which(coname%in%coname1))
  c1<-n11
  c0<-n01
  c1t<-n11+n01
  t1<-c1*(c1-1)/2+c0*(c0-1)/2
  t2<-c1t*(c1t-1)/2
  t3<-2*t1*t2/(n*(n-1))
  ari<-((n11*(n11-1)/2+n01*(n01-1)/2)-t3)/((t1+t2)/2-t3)
  return(ari)
}

ggplot(perf, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = 'red') +
  geom_segment(
    aes(x = 0, y = 0, xend = 1, yend = 1),
    linetype = "dotted",
    color = "grey50"
  ) +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  ggtitle("ROC Curve") +
  annotate("text", 
           x = 0.7, y = 0.2,
           label = paste0('AUC = ', round(auc, 3))) +
  theme_bw() 
