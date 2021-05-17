#GAP procedure
gapprocedure <- function(Ti,Si,alpha=0.1,K=3,N=10,tau=0.5,epsilon=1e-5){
  #Input:
  ##Ti test statistic
  ##Si screen statistic
  ##alpha FDR level default=0.5
  ##K group number default=3
  ##N grid dense level default=10
  ##tau estimate nonnull proportion parameter default=0.5
  ##epsilon Algorithm stability parameter default=1e-5
  
  #Output:
  ##decision
  
  #计算p值
  p=2*(1-pnorm(abs(Ti)))
  #创建grid set
  m=length(Ti)
  
  lambda=c(-Inf,seq(-4*sqrt(log(m)),4*sqrt(log(m)),1/N*sqrt(log(m))),Inf)
  l=length(lambda)-2
  #分组,分成K组,则要从grid set中抽取K-1个元素
  groupgridset=t(combn(1:l,K-1))#每行就是抽取的格点
  #接下来就是逐行遍历
  row=dim(groupgridset)[1]
  col=dim(groupgridset)[2]
  rejection=0*row#记录每个分组的拒绝数
  g=rbind(Ti,Si,rep(0,m),rep(0,m))#检验统计量
  rownames(g)=c("Ti","Si","group","group_num")
  #print("iteration=")
  for(i in 1:row){
    #print(i)
    #对每个分组进行计算rejection
    grididv=groupgridset[i,]+1#分点对应的下标，因为有-Inf在前面
    gridval=lambda[c(1,grididv,length(lambda))]#分组格点
    #进行分组,分配组下标
    for(j in 1:(length(gridval)-1)){
      g[3,g[2,]>=gridval[j]&g[2,]<gridval[j+1]]=j#分配组下标
    }
    #计算组内元素个数
    groupnum=rep(0,K)
    for(j in 1:(length(gridval)-1)){
      g[4,g[2,]>=gridval[j]&g[2,]<gridval[j+1]]=length(g[2,g[2,]>=gridval[j]&g[2,]<gridval[j+1]])
      groupnum[j]=length(g[2,g[2,]>=gridval[j]&g[2,]<gridval[j+1]])
    }
    #估计每个组内nonnull比例
    pi=0*groupnum
    for (j in 1:K){
      pi[j]=1-length(p[p>tau & g[2,]>=gridval[j]&g[2,]<gridval[j+1]])/(1-tau)/m
      #print(pi[j])
      pi[j]=min(max(epsilon,pi[j]),1-epsilon)
    }
    #计算每个组的权重
    weight=m*pi/(1-pi)/sum(groupnum*pi/(1-pi))
    #调整p值
    pw=p
    for (j in 1:K){
      pw[g[3,]==j]=pmin(p[g[3,]==j]/weight[j],1)#这个最小值出问题了
    }
    pw=sort(pw)
    #进行BH方法
    j=m
    while (j > 0){
      if (pw[j]>(j)/m*alpha){
        j=j-1
      } else{
        break
      }
    }
    
    #记录rejection
    rejection[i]=j
  
    #重置g
    g=rbind(Ti,Si,rep(0,m),rep(0,m))
  }
  #选取rejection最大的格点，进而算权重计算rejection，输出decision
  opt_grid=groupgridset[which.max(rejection),]
  gridval=lambda[c(1,opt_grid,length(lambda))]
  for(j in 1:(length(gridval)-1)){
    g[3,g[2,]>=gridval[j]&g[2,]<gridval[j+1]]=j
  }
  groupnum=rep(0,K)
  for(j in 1:(length(gridval)-1)){
    g[4,g[2,]>=gridval[j]&g[2,]<gridval[j+1]]=length(g[2,g[2,]>=gridval[j]&g[2,]<gridval[j+1]])
    groupnum[j]=length(g[2,g[2,]>=gridval[j]&g[2,]<gridval[j+1]])
  }
  pi=0*groupnum
  for (j in 1:K){
    pi[j]=1-length(p[p>tau & g[2,]>=gridval[j]&g[2,]<gridval[j+1]])/(1-tau)/m
    pi[j]=min(max(epsilon,pi[j]),1-epsilon)
  }
  weight=m*pi/(1-pi)/sum(groupnum*pi/(1-pi))
  pw=p
  for (j in 1:K){
    pw[g[3,]==j]=pmin(p[g[3,]==j]/weight[j],1)
  }
  pw=sort(pw)
  j=m
  while (j > 0){
    if (pw[j]>(j)/m*alpha){
      j=j-1
    } else{
      break
    }
  }
  p=rbind(p,order(p),rep(0,m))
  p[3,p[2,]<=j]=1
  return(p[3,])
}



#try

n1=100#样本容量
n2=100#样本容量
m=20000#变量数/检验数
rep=500#重复试验次数
m1=floor(sqrt(m))+100#nonnull的个数
alpha=seq(0.05,1,length.out = 20)#预设FDR水平
mu1=c(rep(5*sqrt(log(m)/n1),m1),rep(0,m-m1))
mu2=c(rep(2*sqrt(log(m)/n2),m1),rep(0,m-m1))

#真实标签
true=c(rep(1,m1),rep(0,m-m1))

fdr=0
for (i in 1:rep){
  #XY比较同方差情形
  X=matrix(mu1,m,n1)+matrix(rnorm(m*n1,0,1),m,n1)#生成X样本
  Y=matrix(mu2,m,n2)+matrix(rnorm(m*n2,0,1),m,n2)#生成Y样本
  #XZ比较异方差情形
  Z=matrix(mu2,m,n2)+matrix(rnorm(m*n2,0,0.5),m,n2)#生成Z样本
  
  #计算检验统计量——t统计量，要算各个变量之间的均值和标准差，注意var函数分母是n-1
  Xbar=apply(X,1,mean)
  Xvar=apply(X,1,var)
  Ybar=apply(Y,1,mean)
  Yvar=apply(Y,1,var)
  Zbar=apply(Z,1,mean)
  Zvar=apply(Z,1,var)
  varpool=(Xvar*(n1-1)+Yvar*(n2-1))/(n1+n2-2)#计算方差整合样本方差
  T1=sqrt(n1*n2/(n1+n2)/varpool)*(Xbar-Ybar)#同方差模型XY检验统计量
  T2=1/sqrt(Xvar/n1+Zvar/n2)*(Xbar-Zbar)#异方差模型XZ检验统计量
  S1=sqrt(n1^2/(n1+n2)/varpool)*(Xbar+n2/n1*Ybar)
  S2=sqrt(n1/(Xvar*(1+n2/n1*Xvar/Yvar)))*(Xbar+n2/n1*Xvar/Yvar*Ybar)
  #计算p值
  p1=(1-pt(T1,n1+n2-2))*2
  
  p2=(1-pnorm(T1,0,1))*2#异方差是渐进正态统计量，
  #关于p值和G，z值的关系可以看GAP supp第九页或者oracle comupoud decision的正文
  
  flag=gapprocedure(T1,S1)
  flag=rbind(flag,true)
  falsediscovery=sum(flag[1,]==1 & flag[2,]==0)
  rejections=sum(flag[1,]==1)
  fdp=falsediscovery/max(rejections,1)
  print(rejections)
  print(fdp)
  fdr=fdr+fdp
}
fdr=fdr/i
fdr
