#Adatptive BH procedure
##m0用Storey的方法估计
abhprocedure_s  <- function(pvalues,alpha,tau) {
  #输入p值向量，预设水平
  #输出拒绝原假设的下标
  n=length(pvalues)#统计检验数
  p=sort(pvalues)#对p值排序
  pp=rbind(pvalues,order(pvalues),rep(0,n))
  rownames(pp)=c("p-values","rank","flag")
  m0=min(length(p[p>tau])/(1-tau),n)
  print(m0)
  #做调整后的BHprocedure 
  pp=bhprocedure(pvalues,alpha*n/m0)
  return(pp)
}


#try
##Model 1 exactly  sparse
n1=100#样本容量
n2=100#样本容量
m=2000#变量数/检验数
rep=500#重复试验次数
m1=floor(sqrt(m))#nonnull的个数
alpha=seq(0.05,1,length.out = 20)#预设FDR水平
mu1=c(rep(10*sqrt(log(m)/n1),m1),rep(0,m-m1))
mu2=c(rep(2*sqrt(log(m)/n2),m1),rep(0,m-m1))

#真实标签
true=c(rep(1,m1),rep(0,m-m1))
fdr=0
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
  T1=sqrt(n1*n2/(n1+n2)/varpool)*abs(Xbar-Ybar)#同方差模型XY检验统计量
  
  #计算p值
  p1=(1-pt(T1,n1+n2-2))*2
  p=abhprocedure_s(p1,0.1,0.9)
  p=rbind(p,true)
  
  falsediscovery=sum(p[3,]==1 & p[4,]==0)
  rejections=sum(p[3,]==1)
  fdp=falsediscovery/max(rejections,1)
  print(fdp)
  fdr=fdr+fdp
}
fdr=fdr/rep
fdr
