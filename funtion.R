library(CARS)
library(qvalue)
library(fdrtool)

#BH Procedure
bhprocedure  <- function(pvalues,alpha) {
  #输入p值向量，预设水平
  #输出拒绝原假设的下标
  n=length(pvalues)#统计检验数
  p=sort(pvalues)#对p值排序
  pp=rbind(pvalues,order(pvalues),rep(0,n))
  rownames(pp)=c("p-values","rank","flag")
  i=n
  
  #不是连贯的 是要最大的成立即可
  while (i > 0){
    if (p[i]>(i)/n*alpha){
      i=i-1
    } else{
      break
    }
  }
  pp[3,which(pp[2,]<=i)]=1
  return(pp)
}


pluginprocedure <- function(pvalues,alpha,lambda){
  
  #估计epsilon
  m=length(pvalues)
  W=m-length(pvalues[pvalues<=lambda])#W(λ)
  epsilon_n=1-min(W/(m*(1-lambda)),0.99)#ε_n，为了稳定，限制在0~1
  
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

abhprocedure  <- function(pvalues,alpha) {
  #输入p值向量，预设水平
  #输出拒绝原假设的下标
  n=length(pvalues)#统计检验数
  p=sort(pvalues)#对p值排序
  pp=rbind(pvalues,order(pvalues),rep(0,n))
  rownames(pp)=c("p-values","rank","flag")
  i=0
  
  while (i <= n-1){
    if (p[i+1]>(i+1)/n*alpha){
      i=i+1
    }
    else{
      break
    }
  }
  if (i==n){
    return(pp)
  }
  else {
    S=0#创建向量S
    for (i in 1:n){
      S[i]=(1-p[i])/(n+1-i)
    }
    #接下来估计m0
    for (i in 2:n){
      if (S[i]<S[i-1]){
        break
      }
    }
    m0=min(m,floor(1/S[i]+1))
    #做调整后的BHprocedure 
    ppp=bhprocedure(pvalues,alpha*m/m0)
    return(ppp)
  }
}

  
  #GFC Procedure
  gfcprocedure  <- function(test_statistic,threshold,alpha) {
    #test_statistic是渐进正态统计量
    #threshold是考虑拒绝域的限制范围,2sqrtlogp或者sqrt(2logp-2loglogp)或者sqrt(2logp)
    #alpha是控制的FDR水平
    #输出拒绝原假设的下标
    n=length(test_statistic)#统计检验数
    p=sort(test_statistic)#对p值排序
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
  
  
  usprocedure <- function(Ti,Si,alpha,N){
    #Ti test statistic
    #Si screen statistic
    #alpha FDR level
    #grid number,N=10
    
    ##Step 1:计算检验统计量，输入Ti，Si
    m=length(Ti)
    p=2*(1-pnorm(abs(Ti)))
    p=rbind(p,order(p),order(p),rep(0,m))
    p0=sort(p[1,])
    ##Step 2:调参，λ∈(0,4sqrt(log(m)))
    lambda=seq(0,4*sqrt(log(m)),1/N*sqrt(log(m)))#参数范围
    rejection=0*lambda#记录各个参数下的拒绝数
    K=length(lambda)
    j=1
    k=m
    while (k > 0){
      if (p0[k]>(k)/m*alpha){
        k=k-1
      } else{
        break
      }
    }
    rejection[1]=k
    p[4,]=rep(0,m)
    for (j in 2:K){
      #分组数
      parameter=lambda[j]
      m1=length(Si[abs(Si)<=parameter])
      m2=m-m1
      #处理合并p1p2的flag，多两行order去做bh，就可以不影响
      p[2,abs(Si)<=parameter]=order(p[1,abs(Si)<=parameter])
      p[2,abs(Si)>parameter]=m+1
      
      p[3,abs(Si)>parameter]=order(p[1,abs(Si)>parameter])
      p[3,abs(Si)<=parameter]=m+1
      
      #由GAP supplementary lemma3 得出求that的过程等价于BH
      i=m1
      
      p1=p0[abs(Si)<=parameter]
      #print(length(p1))
      while (i > 0){
        if (p1[i]>(i)/m1*alpha){
          i=i-1
        }
        else{
          break
        }
      }
      p[4,which(p[2,]<=i)]=1
      
      i=m2
      p2=p0[abs(Si)>parameter]
      #print(length(p2))
      while (i > 0){
        if (p2[i]>(i)/m2*alpha){
          i=i-1
        }
        else{
          break
        }
      }
      p[4,which(p[3,]<=i)]=1
      rejection[j]=sum(p[4,])
      p[4,]=rep(0,m)
    }
    ##Step 3:选取最优参数
    p[4,]=rep(0,m)
    opt_lambda=lambda[which.max(rejection)]
    m1=length(Si[abs(Si)<=opt_lambda])
    m2=m-m1
    p[2,abs(Si)<=opt_lambda]=order(p[1,abs(Si)<=opt_lambda])
    p[2,abs(Si)>opt_lambda]=m+1
    
    p[3,abs(Si)>opt_lambda]=order(p[1,abs(Si)>opt_lambda])
    p[3,abs(Si)<=opt_lambda]=m+1
    
    i=m1
    p1=p0[abs(Si)<=opt_lambda]
    while (i > 0){
      if (p1[i]>(i)/m1*alpha){
        i=i-1
      } else{
        break
      }
    }
    p[4,which(p[2,]<=i)]=1
    
    i=m2
    p2=p0[abs(Si)>opt_lambda]
    while (i > 0){
      if (p2[i]>(i)/m2*alpha){
        i=i-1
      } else{
        break
      }
    }
    p[4,which(p[3,]<=i)]=1
    ##Step 4:输出结果
    return(p[4,])
  }
  
  
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
  