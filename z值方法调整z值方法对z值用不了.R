library(fdrtool)
library(qvalue)
n1=100#样本容量
n2=100#样本容量
m=2000#变量数/检验数
rep=500#重复试验次数
m1=floor(sqrt(m))#nonnull的个数
m1=440
alpha=seq(0.05,1,length.out = 20)#预设FDR水平
mu1=c(rep(4*sqrt(log(m)/n1),m1),rep(0,m-m1))
mu2=c(rep(2*sqrt(log(m)/n2),m1),rep(0,m-m1))

#采用zvalues方法
true=c(rep(1,m1),rep(0,m-m1))
fdr2=0
alpha=0.1
  #XY比较同方差情形
  X=matrix(mu1,m,n1)+matrix(rnorm(m*n1,0,1),m,n1)#生成X样本
  Y=matrix(mu2,m,n2)+matrix(rnorm(m*n2,0,1),m,n2)#生成Y样本
  
  #计算检验统计量——t统计量，要算各个变量之间的均值和标准差，注意var函数分母是n-1
  Xbar=apply(X,1,mean)
  Xvar=apply(X,1,var)
  Ybar=apply(Y,1,mean)
  Yvar=apply(Y,1,var)
  varpool=(Xvar*(n1-1)+Yvar*(n2-1))/(n1+n2-2)#计算方差整合样本方差
  T1=sqrt(n1*n2/(n1+n2)/varpool)*abs(Xbar-Ybar)#同方差模型XY检验统计量
  cdf=ecdf(T1)
  T2=cdf(T1)
  z=qnorm(T2)
  q2=adazprocedure(z,alpha,"normal")
  q2=rbind(q2,true)
  
  p1=(1-pt(abs(T1),n1+n2-2))*2
  q1=adazprocedure(p1,alpha,"pvalue")
  q1=rbind(q1,true)
  
  falsediscovery=sum(q2[3,]==1 & q2[4,]==0)
  falsediscovery
  
  rejections=sum(q2[3,]==1)
  rejections
  
  fdp2=falsediscovery/max(rejections,1)
  fdp2
  
  falsediscovery=sum(q1[3,]==1 & q1[4,]==0)
  falsediscovery
  
  rejections=sum(q1[3,]==1)
  rejections
  
  fdp1=falsediscovery/max(rejections,1)
  fdp1
  
  