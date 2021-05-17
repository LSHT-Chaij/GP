library(fdrtool)
library(qvalue)

adazprocedure <- function(zvalues,alpha,poz){
  #zvalue 有点问题，但是p值的问题不大
  m=length(zvalues)
  z=rbind(zvalues,order(zvalues),rep(0,m))
  if (poz=="normal"){
    rownames(z)=c("z-values","rank","flag")
    q<-fdrtool(sort(zvalues),"normal",plot=FALSE)
  } else if (poz=="pvalue"){
    rownames(z)=c("p-values","rank","flag")
    q<-fdrtool(sort(zvalues),"pvalue",plot=FALSE)
  }
  lfdr=sort(q$lfdr)
  i=m
  while (i<=m) {
    Q=sum(lfdr[1:i])/i
    if (Q>alpha){
      i=i-1
    } else{
      break
    }
    if(i==0){
      break
    }
  }
  z[3,z[2,]<=i]=1
  return(z)
}

adazprocedure_q <- function(pvalues,alpha){
  m=length(pvalues)
  z=rbind(pvalues,order(pvalues),rep(0,m))
  rownames(z)=c("p-values","rank","flag")
  
  q<-qvalue(sort(pvalues))
  lfdr=sort(q$lfdr)
  i=m
  while (i<=m) {
    Q=sum(lfdr[1:i])/i
    if (Q>alpha){
      i=i-1
    } else{
      break
    }
  }
  z[3,z[2,]<=i]=1
  return(z)
}

#qvalue 里的 参数是用Storey方法估计的 只能对p值用，因为不知道分布信息，所以不能采用p值
#fdrtool 里面用的是censorfit估计的 可对z值使用
#用两个方法看看结果差异大不大

#p值转z值
n1=1000#样本容量
n2=1000#样本容量
m=2000#变量数/检验数
rep=500#重复试验次数
m1=floor(sqrt(m))#nonnull的个数
m1=44
alpha=seq(0.05,1,length.out = 20)#预设FDR水平
mu1=c(rep(5*sqrt(log(m)/n1),m1),rep(0,m-m1))
mu2=c(rep(2*sqrt(log(m)/n2),m1),rep(0,m-m1))

#采用zvalues方法
true=c(rep(1,m1),rep(0,m-m1))
fdr1=0
fdr2=0
alpha=0.1
for (i in 1:rep){
  #XY比较同方差情形
  X=matrix(mu1,m,n1)+matrix(rnorm(m*n1,0,1),m,n1)#生成X样本
  Y=matrix(mu2,m,n2)+matrix(rnorm(m*n2,0,1),m,n2)#生成Y样本
  
  #计算检验统计量——t统计量，要算各个变量之间的均值和标准差，注意var函数分母是n-1
  Xbar=apply(X,1,mean)
  Xvar=apply(X,1,var)
  Ybar=apply(Y,1,mean)
  Yvar=apply(Y,1,var)
  varpool=(Xvar*(n1-1)+Yvar*(n2-1))/(n1+n2-2)#计算方差整合样本方差
  T1=sqrt(n1*n2/(n1+n2)/varpool)*(Xbar-Ybar)#同方差模型XY检验统计量
  p1=(1-pt(abs(T1),n1+n2-2))*2
  q1=adazprocedure(p1,alpha,"pvalue")
  #q=adazprocedure_q(p1,alpha)
  cdf=ecdf(T1)
  T2=cdf(T1)
  #z=qnorm(T2)
  z=(T1-mean(T1))/sd(T1)
  q2=adazprocedure(T1,alpha,"normal")

  q1=rbind(q1,true)
  q2=rbind(q2,true)
  
  falsediscovery=sum(q1[3,]==1 & q1[4,]==0)
  rejections=sum(q1[3,]==1)
  fdp1=falsediscovery/max(rejections,1)
  print(fdp1)
  fdr1=fdr1+fdp1
  
  falsediscovery=sum(q2[3,]==1 & q2[4,]==0)
  rejections=sum(q2[3,]==1)
  fdp2=falsediscovery/max(rejections,1)
  print(fdp2)
  fdr2=fdr2+fdp2
}
fdr1=fdr1/rep
fdr1
fdr2=fdr2/rep
fdr2


