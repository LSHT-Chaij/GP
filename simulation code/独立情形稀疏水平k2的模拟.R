#模型设定成功：比较cars us bh az qv还行 gap较差
library(fdrtool)
library(qvalue)
library(CARS)
#k2 稀疏水平的变化

n1=50
n2=60
varx=1
vary=4
m=2000
alpha=0.1
rep=150

k1=1000#k1表示非零元的下标终点
k2=200#k2=seq(100,1000,100) 表示差别元的下标终点
mu1=0
mu2=0
for (i in 1:m){
  if (i %in% 1:k2){
    mu1[i]=7/sqrt(30)
    mu2[i]=2/sqrt(30)
  }
  else if (i %in% (k2+1):k1){
    mu1[i]=4/sqrt(30)
    mu2[i]=4/sqrt(30)
  }
  else {
    mu1[i]=0
    mu2[i]=0
  }
}

K2=seq(100,1000,100)
k1=1000#k1表示非零元的下标终点
N=length(K2)
alpha=0.1
true=c(rep(1,k2),rep(0,m-k2))
#fdrlevel=seq(0.05,1,0.05)#预设FDR水平
rejections=matrix(0,9,rep)
fdp=matrix(0,9,rep)
falsediscovery=matrix(0,9,rep)
rownames(fdp)=c("BH","Adaptive BH","Plugin","qvalue","Adaptive zvalue","GFC","US","CARS","GAP")
rownames(rejections)=c("BH","Adaptive BH","Plugin","qvalue","Adaptive zvalue","GFC","US","CARS","GAP")
rownames(falsediscovery)=c("BH","Adaptive BH","Plugin","qvalue","Adaptive zvalue","GFC","US","CARS","GAP")
N01=fdp
S=fdp
e=fdp
fdr=matrix(0,9,N)
rownames(fdr)=c("BH","Adaptive BH","Plugin","qvalue","Adaptive zvalue","GFC","US","CARS","GAP")

etp=fdr
fnr=fdr


#for(alpha in fdrlevel){
#  print(alpha)
#}
#对alpha循环





for (j in 1:length(K2)){
  k2=K2[j]#k2=seq(100,1000,100) 表示差别元的下标终点
  mu1=0
  mu2=0
  for (i in 1:m){
    if (i %in% 1:k2){
      mu1[i]=7/sqrt(30)
      mu2[i]=2/sqrt(30)
    }
    else if (i %in% (k2+1):k1){
      mu1[i]=4/sqrt(30)
      mu2[i]=4/sqrt(30)
    }
    else {
      mu1[i]=0
      mu2[i]=0
    }
  }
  
  true=c(rep(1,k2),rep(0,m-k2))
  #fdrlevel=seq(0.05,1,0.05)#预设FDR水平
  rejections=matrix(0,9,rep)
  fdp=matrix(0,9,rep)
  falsediscovery=matrix(0,9,rep)
  rownames(fdp)=c("BH","Adaptive BH","Plugin","qvalue","Adaptive zvalue","GFC","US","CARS","GAP")
  rownames(rejections)=c("BH","Adaptive BH","Plugin","qvalue","Adaptive zvalue","GFC","US","CARS","GAP")
  rownames(falsediscovery)=c("BH","Adaptive BH","Plugin","qvalue","Adaptive zvalue","GFC","US","CARS","GAP")
  N01=fdp
  S=fdp
  e=fdp
  
  
  
  #for(alpha in fdrlevel){
  #  print(alpha)
  #}
  #对alpha循环
  
  
  m1=k2
  
  print(paste("the slevel is",k2/k1))
  #内循环
  for (i in 1:rep){
    ##数据生成
    print(paste("the iteration is",i))
    #XY比较同方差情形
    X=matrix(mu1,m,n1)+matrix(rnorm(m*n1,0,sdx),m,n1)#生成X样本
    Y=matrix(mu2,m,n2)+matrix(rnorm(m*n2,0,sdy),m,n2)#生成Y样本
    
    Xbar=apply(X,1,mean)
    Xvar=apply(X,1,var)
    Ybar=apply(Y,1,mean)
    Yvar=apply(Y,1,var)
    T1=1/sqrt(Xvar/n1+Yvar/n2)*(Xbar-Ybar)#异方差模型XZ检验统计量
    S1=sqrt(n1/(Xvar*(1+n2/n1*Xvar/Yvar)))*(Xbar+n2/n1*Xvar/Yvar*Ybar)
    p1=(1-pnorm(abs(T1),0,1))*2#异方差是渐进正态统计量，
    #关于p值和G，z值的关系可以看GAP supp第九页或者oracle comupoud decision的正文
    
    ##计算fdp等指标
    #BH
    p=bhprocedure(p1,alpha)
    p=rbind(p,true)
    
    falsediscovery[1,i]=sum(p[3,]==1 & p[4,]==0)
    rejections[1,i]=sum(p[3,]==1)
    N01[1,i]=sum(p[3,]==0 & p[4,]==1)
    S[1,i]=sum(p[3,]==0)
    e[1,i]=sum(p[3,]==1 & p[4,]==1)/m1
    print("BH is over!")
    #ABH
    p=abhprocedure(p1,alpha)
    p=rbind(p,true)
    
    #QV
    q=qvalue(p1)
    q=rbind(q$qvalues,order(q$qvalues),rep(0,length(q$qvalues)))
    q[3,q[1,]<=alpha]=1
    q=rbind(q,true)
    falsediscovery[4,i]=sum(q[3,]==1 & q[4,]==0)
    rejections[4,i]=sum(q[3,]==1)
    N01[4,i]=sum(q[3,]==0 & q[4,]==1)
    S[4,i]=sum(q[3,]==0)
    e[4,i]=sum(q[3,]==1 & q[4,]==1)/m1
    print("QV is over!")
    #AZ
    q1=adazprocedure(p1,alpha,"pvalue")
    q1=rbind(q1,true)
    
    falsediscovery[5,i]=sum(q1[3,]==1 & q1[4,]==0)
    rejections[5,i]=sum(q1[3,]==1)
    N01[5,i]=sum(q1[3,]==0 & q1[4,]==1)
    S[5,i]=sum(q1[3,]==0)
    e[5,i]=sum(q1[3,]==1 & q1[4,]==1)/m1
    print("AZ is over!")
    #US
    flag=usprocedure(T1,S1,alpha,10)
    flag=rbind(flag,true)
    falsediscovery[7,i]=sum(flag[1,]==1 & flag[2,]==0)
    rejections[7,i]=sum(flag[1,]==1)
    N01[7,i]=sum(flag[1,]==0 & flag[2,]==1)
    S[7,i]=sum(flag[1,]==0)
    e[7,i]=sum(flag[1,]==1 & flag[2,]==1)/m1
    print("US is over!")
    
    
    
    #GAP
    flag=gapprocedure(T1,S1,alpha)
    flag=rbind(flag,true)
    falsediscovery[9,i]=sum(flag[1,]==1 & flag[2,]==0)
    rejections[9,i]=sum(flag[1,]==1)
    N01[9,i]=sum(flag[1,]==0 & flag[2,]==1)
    S[9,i]=sum(flag[1,]==0)
    e[9,i]=sum(flag[1,]==1 & flag[2,]==1)/m1
    print("GAP is over!")
    
    #CARS 也有可能连续两次出错，可以写一个递归函数，直到成功！
    fit<-try(CARS(X,Y,alpha,tau=0.2,option='regular'),silent=TRUE)
    if('try-error' %in% class(fit)){
      next
    }else{
      flag=CARS(X,Y,alpha,tau=0.2,option='regular')
      flag=flag$de
      flag=rbind(flag,true)
      falsediscovery[8,i]=sum(flag[1,]==1 & flag[2,]==0)
      rejections[8,i]=sum(flag[1,]==1)
      N01[8,i]=sum(flag[1,]==0 & flag[2,]==1)
      S[8,i]=sum(flag[1,]==0)
      e[8,i]=sum(flag[1,]==1 & flag[2,]==1)/m1
    }
    print("CARS is over!")
  }
  fdp=falsediscovery/pmax(rejections,1)
  fnp=N01/pmax(S,1)
  fdr[,j]=apply(fdp,1,mean)#行列要考虑一下
  fnr[,j]=apply(fnp,1,mean)Z
  etp[,j]=apply(e,1,mean)
}
write.csv(fnr,"C:/Users/Lenovo/Desktop/结果/fnr2.csv")
write.csv(etp,"C:/Users/Lenovo/Desktop/结果/etp2.csv")
write.csv(fdr,"C:/Users/Lenovo/Desktop/结果/fdr2.csv")
