#GFC Procedure
gfcprocedure  <- function(test_statistic,threshold,alpha) {
  #test_statistic是渐进正态统计量
  #threshold是考虑拒绝域的限制范围,2sqrtlogp或者sqrt(2logp-2loglogp)或者sqrt(2logp)
  #alpha是控制的FDR水平
  #输出拒绝原假设的下标
  test_statistic=abs(test_statistic)
  n=length(test_statistic)#统计检验数
  p=sort(abs(test_statistic))#对p值排序
  pp=rbind(test_statistic,order(test_statistic),rep(0,n))
  rownames(pp)=c("test_statistic","rank","flag")

  #计算that
  grid=seq(0,threshold,length.out = 10000)
  fdp_est=grid*0
  rejections=length(p[p<=threshold])
  bp=(2-2*pnorm(threshold))*n/max(rejections,1)
  if (bp<=alpha){
    for(i in 1:10000){
    rejections=length(p[p<=grid[i]])
    Gtailprop=2-2*pnorm(grid[i])
    fdp_est[i]=Gtailprop*n/max(rejections,1)
  }
  that=grid[which.min(fdp_est[fdp_est<=alpha])]
  } else{
    that=bp
  }
  #开始统计拒绝的假设
  pp[3,which(abs(pp[1,])>that)]=1
  return(pp)
}


#try
##Model 1 exactly  sparse
n1=100#样本容量
n2=100#样本容量
m=2000#变量数/检验数
rep=500#重复试验次数
m1=200#nonnull的个数
#m1=500,看看不稀疏的时候他和abh的差别
alpha=seq(0.05,1,length.out = 20)#预设FDR水平
mu1=c(rep(4*sqrt(log(m)/n1),m1),rep(0,m-m1))
mu2=c(rep(2*sqrt(log(m)/n2),m1),rep(0,m-m1))

#真实标签
true=c(rep(1,m1),rep(0,m-m1))
fdr=0
for (i in 1:50000){
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
  
  #计算p值
  test=gfcprocedure(T1,2*sqrt(log(m)),0.3)
  test=rbind(test,true)
  
  falsediscovery=sum(test[3,]==1 & test[4,]==0)
  rejections=sum(test[3,]==1)
  fdp=falsediscovery/max(rejections,1)
  print(rejections)
  print(fdp)
  print((rejections-falsediscovery)/m1)
  fdr=fdr+fdp
}
fdr=fdr/rep
fdr
