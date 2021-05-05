# Oracle and Plug-in p-values Procedure
#在稀疏情形不可用，因为本来是用G去估计G1
#稀疏时候基本上都是G0，G0个G1又有明显不同，所以会产生很大的偏差。
pluginprocedure <- function(pvalues,alpha,lambda){
 
  #估计epsilon
  m=length(pvalues)
  W=m-length(pvalues[pvalues<=lambda])#W(λ)
  epsilon_n=1-min(W/(m*(1-lambda)),0.99)#ε_n，为了稳定，限制在0~1
  #print(epsilon_n)
  #epsilon_n=0.022
  
  #估计G
  emprical_cdf=ecdf(pvalues) 
  
  p=sort(pvalues)
  pp=rbind(pvalues,order(pvalues),rep(0,m))
  rownames(pp)=c("p-values","rank","flag")
  #求阈值t的估计，求函数的上确界
  #首先先判断t在不在0~1，t>1的情形是线性函数，直接求解即可，0~1遍历
  #也就是判断1-epsilon和alpha的大小关系
  if (1-epsilon_n>alpha){
    #说明t在0~1
    t=seq(0,1,length.out = 10001)
    find=pluginQ(epsilon_n,emprical_cdf,t)<=alpha
    that=t[10001+1-which.max(rev(find))]
  } else{
    that=alpha/(1-epsilon_n)
  }
  #print(that)
  pp[3,pp[1,]<=that]=1
  return(pp)
}

#emprical Q function 
pluginQ <- function(epsilon_n,emprical_cdf,t){
  #返回Q
  return((1-epsilon_n)*t/max(emprical_cdf(t),1e-5))#保持算法稳定性
}

n1=100#样本容量
n2=100#样本容量
m=2000#变量数/检验数
rep=50#重复试验次数
m1=floor(sqrt(m))#nonnull的个数
M1=seq(50,1500,50)
fdr=0*M1
for (j in 1:length(M1)){
  m1=M1[j]
  print(paste("iter=",j))
  mu1=c(rep(7*sqrt(log(m)/n1),m1),rep(0,m-m1))
  mu2=c(rep(2*sqrt(log(m)/n2),m1),rep(0,m-m1))
  true=c(rep(1,m1),rep(0,m-m1))
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
    p=pluginprocedure(p1,0.1,0.5)
    p=rbind(p,true)
    
    falsediscovery=sum(p[3,]==1 & p[4,]==0)
    rejections=sum(p[3,]==1)
    fdp=falsediscovery/max(rejections,1)
    fdr[j]=fdr[j]+fdp
  }
  fdr[j]=fdr[j]/rep
}
m1=44
alpha=seq(0.05,1,length.out = 20)#预设FDR水平


#真实标签


