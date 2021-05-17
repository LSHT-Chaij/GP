library(CARS)
library(qvalue)
library(fdrtool)
#Model 1 simulation：equal variance or unequal variance
n1=100#样本容量
n2=100#样本容量
m=2000#变量数/检验数
rep=50#重复试验次数
m1=200#nonnull的个数
fdrlevel=seq(0.05,1,0.05)#预设FDR水平
mu1=c(rep(4.5*sqrt(log(m)/n1),m1),rep(0,m-m1))
mu2=c(rep(2*sqrt(log(m)/n2),m1),rep(0,m-m1))
#真实标签
true=c(rep(1,m1),rep(0,m-m1))

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

N=length(fdrlevel)
fdr=matrix(0,9,N)
rownames(fdr)=c("BH","Adaptive BH","Plugin","qvalue","Adaptive zvalue","GFC","US","CARS","GAP")

etp=fdr
fnr=fdr
for (FDR in 1:N){
  alpha=fdrlevel[FDR]
  print(paste("the fdrlevel is",alpha))
  #内循环
  for (i in 1:rep){
    ##数据生成
    print(paste("the iteration is",i))
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
    p1=(1-pt(abs(T1),n1+n2-2))*2
    
    p2=(1-pnorm(T1,0,1))*2#异方差是渐进正态统计量，
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
    
    falsediscovery[2,i]=sum(p[3,]==1 & p[4,]==0)
    rejections[2,i]=sum(p[3,]==1)
    N01[2,i]=sum(p[3,]==0 & p[4,]==1)
    S[2,i]=sum(p[3,]==0)
    e[2,i]=sum(p[3,]==1 & p[4,]==1)/m1
    print("ABH is over!")
    #PI
    p=pluginprocedure(p1,alpha,0.5)#0.5是调参
    p=rbind(p,true)
    
    falsediscovery[3,i]=sum(p[3,]==1 & p[4,]==0)
    rejections[3,i]=sum(p[3,]==1)
    N01[3,i]=sum(p[3,]==0 & p[4,]==1)
    S[3,i]=sum(p[3,]==0)
    e[3,i]=sum(p[3,]==1 & p[4,]==1)/m1
    print("PI is over!")
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
    #GFC
    test=gfcprocedure(T1,2*sqrt(log(m)),alpha)
    test=rbind(test,true)
    falsediscovery[6,i]=sum(test[3,]==1 & test[4,]==0)
    rejections[6,i]=sum(test[3,]==1)
    N01[6,i]=sum(test[3,]==0 & test[4,]==1)
    S[6,i]=sum(test[3,]==0)
    e[6,i]=sum(test[3,]==1 & test[4,]==1)/m1
    print("GFC is over!")
    #US
    flag=usprocedure(T1,S1,alpha,10)
    flag=rbind(flag,true)
    falsediscovery[7,i]=sum(flag[1,]==1 & flag[2,]==0)
    rejections[7,i]=sum(flag[1,]==1)
    N01[7,i]=sum(flag[1,]==0 & flag[2,]==1)
    S[7,i]=sum(flag[1,]==0)
    e[7,i]=sum(flag[1,]==1 & flag[2,]==1)/m1
    print("US is over!")
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
    
    
    #GAP
    flag=gapprocedure(T1,S1,alpha)
    flag=rbind(flag,true)
    falsediscovery[9,i]=sum(flag[1,]==1 & flag[2,]==0)
    rejections[9,i]=sum(flag[1,]==1)
    N01[9,i]=sum(flag[1,]==0 & flag[2,]==1)
    S[9,i]=sum(flag[1,]==0)
    e[9,i]=sum(flag[1,]==1 & flag[2,]==1)/m1
    print("GAP is over!")
  }
  fdp=falsediscovery/pmax(rejections,1)
  fnp=N01/pmax(S,1)
  fdr[,FDR]=apply(fdp,1,mean)#行列要考虑一下
  fnr[,FDR]=apply(fnp,1,mean)
  etp[,FDR]=apply(e,1,mean)
}

library(xlsx)
write.xlsx(fnr,"D:/蔡健/毕设安排/数值模拟代码/数据结果/fnr.xlsx")
write.xlsx(etp,"D:/蔡健/毕设安排/数值模拟代码/数据结果/etp.xlsx")
write.xlsx(fdr,"D:/蔡健/毕设安排/数值模拟代码/数据结果/fdr.xlsx")


#Model 2 simulation：equal variance or unequal variance
n1=100#样本容量
n2=100#样本容量
m=2000#变量数/检验数
rep=10#重复试验次数
m1=200#nonnull的个数
fdrlevel=seq(0.05,1,0.05)#预设FDR水平
mu1=c(rep(2.5*sqrt(log(m)/n1),m1),rep(0,m-m1))
mu2=c(rep(1*sqrt(log(m)/n2),floor(m1/2)),rep(-0.5*sqrt(log(m)/n2),m1-floor(m1/2)),rep(0,m-m1))
#真实标签
true=c(rep(1,m1),rep(0,m-m1))

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

N=length(fdrlevel)
fdr=matrix(0,9,N)
rownames(fdr)=c("BH","Adaptive BH","Plugin","qvalue","Adaptive zvalue","GFC","US","CARS","GAP")

etp=fdr
fnr=fdr
for (FDR in 1:N){
  alpha=fdrlevel[FDR]
  print(paste("the fdrlevel is",alpha))
  #内循环
  for (i in 1:rep){
    ##数据生成
    print(paste("the iteration is",i))
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
    p1=(1-pt(abs(T1),n1+n2-2))*2
    
    p2=(1-pnorm(T1,0,1))*2#异方差是渐进正态统计量，
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
    
    falsediscovery[2,i]=sum(p[3,]==1 & p[4,]==0)
    rejections[2,i]=sum(p[3,]==1)
    N01[2,i]=sum(p[3,]==0 & p[4,]==1)
    S[2,i]=sum(p[3,]==0)
    e[2,i]=sum(p[3,]==1 & p[4,]==1)/m1
    print("ABH is over!")
    #PI
    p=pluginprocedure(p1,alpha,0.5)#0.5是调参
    p=rbind(p,true)
    
    falsediscovery[3,i]=sum(p[3,]==1 & p[4,]==0)
    rejections[3,i]=sum(p[3,]==1)
    N01[3,i]=sum(p[3,]==0 & p[4,]==1)
    S[3,i]=sum(p[3,]==0)
    e[3,i]=sum(p[3,]==1 & p[4,]==1)/m1
    print("PI is over!")
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
    #GFC
    test=gfcprocedure(T1,2*sqrt(log(m)),alpha)
    test=rbind(test,true)
    falsediscovery[6,i]=sum(test[3,]==1 & test[4,]==0)
    rejections[6,i]=sum(test[3,]==1)
    N01[6,i]=sum(test[3,]==0 & test[4,]==1)
    S[6,i]=sum(test[3,]==0)
    e[6,i]=sum(test[3,]==1 & test[4,]==1)/m1
    print("GFC is over!")
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
  fdr[,FDR]=apply(fdp,1,mean)#行列要考虑一下
  fnr[,FDR]=apply(fnp,1,mean)
  etp[,FDR]=apply(e,1,mean)
}

library(xlsx)
write.xlsx(fnr,"D:/蔡健/毕设安排/数值模拟代码/数据结果/fnr2.xlsx")
write.xlsx(etp,"D:/蔡健/毕设安排/数值模拟代码/数据结果/etp2.xlsx")
write.xlsx(fdr,"D:/蔡健/毕设安排/数值模拟代码/数据结果/fdr2.xlsx")




##Model 3 asymptotically sparse
n1=100#样本容量
n2=100#样本容量
m=2000#变量数/检验数
rep=20#重复试验次数
m1=200#nonnull的个数
alpha=seq(0.05,1,length.out = 20)#预设FDR水平

#设置均值向量mu1和mu2
mu1=0#初始化
mu2=0#初始化
for (i in 1:m){
  if (i <= m1){
    mu1[i]=4*sqrt(log(m)/n1)
    mu2[i]=2*sqrt(log(m)/n2)
  }
  else {
    mu1[i]=(i/m)*sqrt(log(m)/n1)
    mu2[i]=(i/m)*sqrt(log(m)/n1)#要求了两者是一样的
  }
}

#XY比较同方差情形
X=matrix(mu1,m,n1)+matrix(rnorm(m*n1,0,1),m,n1)#生成X样本
Y=matrix(mu2,m,n2)+matrix(rnorm(m*n2,0,1),m,n2)#生成Y样本
#XZ比较异方差情形
X=matrix(mu1,m,n1)+matrix(rnorm(m*n1,0,1),m,n1)#生成X样本
Z=matrix(mu2,m,n2)+matrix(rnorm(m*n2,0,0.5),m,n2)#生成Z样本

true=c(rep(1,m1),rep(0,m-m1))

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

N=length(fdrlevel)
fdr=matrix(0,9,N)
rownames(fdr)=c("BH","Adaptive BH","Plugin","qvalue","Adaptive zvalue","GFC","US","CARS","GAP")

etp=fdr
fnr=fdr
for (FDR in 17:N){
  alpha=fdrlevel[FDR]
  print(paste("the fdrlevel is",alpha))
  #内循环
  for (i in 1:rep){
    ##数据生成
    print(paste("the iteration is",i))
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
    p1=(1-pt(abs(T1),n1+n2-2))*2
    
    p2=(1-pnorm(T1,0,1))*2#异方差是渐进正态统计量，
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
    
    falsediscovery[2,i]=sum(p[3,]==1 & p[4,]==0)
    rejections[2,i]=sum(p[3,]==1)
    N01[2,i]=sum(p[3,]==0 & p[4,]==1)
    S[2,i]=sum(p[3,]==0)
    e[2,i]=sum(p[3,]==1 & p[4,]==1)/m1
    print("ABH is over!")
    #PI
    p=pluginprocedure(p1,alpha,0.5)#0.5是调参
    p=rbind(p,true)
    
    falsediscovery[3,i]=sum(p[3,]==1 & p[4,]==0)
    rejections[3,i]=sum(p[3,]==1)
    N01[3,i]=sum(p[3,]==0 & p[4,]==1)
    S[3,i]=sum(p[3,]==0)
    e[3,i]=sum(p[3,]==1 & p[4,]==1)/m1
    print("PI is over!")
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
    #GFC
    test=gfcprocedure(T1,2*sqrt(log(m)),alpha)
    test=rbind(test,true)
    falsediscovery[6,i]=sum(test[3,]==1 & test[4,]==0)
    rejections[6,i]=sum(test[3,]==1)
    N01[6,i]=sum(test[3,]==0 & test[4,]==1)
    S[6,i]=sum(test[3,]==0)
    e[6,i]=sum(test[3,]==1 & test[4,]==1)/m1
    print("GFC is over!")
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
  fdr[,FDR]=apply(fdp,1,mean)#行列要考虑一下
  fnr[,FDR]=apply(fnp,1,mean)
  etp[,FDR]=apply(e,1,mean)
}

library(xlsx)
write.xlsx(fnr,"D:/蔡健/毕设安排/数值模拟代码/数据结果/fnr3.xlsx")
write.xlsx(etp,"D:/蔡健/毕设安排/数值模拟代码/数据结果/etp3.xlsx")
write.xlsx(fdr,"D:/蔡健/毕设安排/数值模拟代码/数据结果/fdr3.xlsx")

##出错，试验cars的数据


#Model 5 CARS setting 1
n1=50
n2=60
sdx=1
sdy=2
m=5000
alpha=0.05
rep=20
k=1000#k=seq(100,1000,100) 表示稀疏程度
mu1=0
mu2=0

for (i in 1:m){
  if (i %in% 1:k){
    mu1[i]=9/sqrt(30)
    mu2[i]=2/sqrt(30)
  }
  else if (i %in% (k+1):(2*k)){
    mu1[i]=4/sqrt(30)
    mu2[i]=4/sqrt(30)
  }
  else {
    mu1[i]=0
    mu2[i]=0
  }
}#差别的比例为k/m，非零比例为2k/m
true=c(rep(1,k),rep(0,m-k))
fdrlevel=seq(0.05,1,0.05)#预设FDR水平
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

N=length(fdrlevel)
fdr=matrix(0,9,N)
rownames(fdr)=c("BH","Adaptive BH","Plugin","qvalue","Adaptive zvalue","GFC","US","CARS","GAP")

etp=fdr
fnr=fdr



alpha=0.1
i=1
for (FDR in 1:N){
  alpha=fdrlevel[FDR]
  print(paste("the fdrlevel is",alpha))
  #内循环
  for (i in 1:20){
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
    
    falsediscovery[2,i]=sum(p[3,]==1 & p[4,]==0)
    rejections[2,i]=sum(p[3,]==1)
    N01[2,i]=sum(p[3,]==0 & p[4,]==1)
    S[2,i]=sum(p[3,]==0)
    e[2,i]=sum(p[3,]==1 & p[4,]==1)/m1
    print("ABH is over!")
    #PI
    p=pluginprocedure(p1,alpha,0.5)#0.5是调参
    p=rbind(p,true)
    
    falsediscovery[3,i]=sum(p[3,]==1 & p[4,]==0)
    rejections[3,i]=sum(p[3,]==1)
    N01[3,i]=sum(p[3,]==0 & p[4,]==1)
    S[3,i]=sum(p[3,]==0)
    e[3,i]=sum(p[3,]==1 & p[4,]==1)/m1
    print("PI is over!")
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
    #GFC
    test=gfcprocedure(T1,2*sqrt(log(m)),alpha)
    test=rbind(test,true)
    falsediscovery[6,i]=sum(test[3,]==1 & test[4,]==0)
    rejections[6,i]=sum(test[3,]==1)
    N01[6,i]=sum(test[3,]==0 & test[4,]==1)
    S[6,i]=sum(test[3,]==0)
    e[6,i]=sum(test[3,]==1 & test[4,]==1)/m1
    print("GFC is over!")
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
  fdr[,FDR]=apply(fdp,1,mean)#行列要考虑一下
  fnr[,FDR]=apply(fnp,1,mean)
  etp[,FDR]=apply(e,1,mean)
}

library(xlsx)
write.xlsx(fnr,"D:/蔡健/毕设安排/数值模拟代码/数据结果/fnrcars1.xlsx")
write.xlsx(etp,"D:/蔡健/毕设安排/数值模拟代码/数据结果/etpcars1.xlsx")
write.xlsx(fdr,"D:/蔡健/毕设安排/数值模拟代码/数据结果/fdrcars1.xlsx")


#Model 6
n1=50
n2=60
varx=1
vary=4
m=2000
alpha=0.1
rep=5

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



true=c(rep(1,k2),rep(0,m-k2))
fdrlevel=seq(0.05,1,0.05)#预设FDR水平
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

N=length(fdrlevel)
fdr=matrix(0,9,N)
rownames(fdr)=c("BH","Adaptive BH","Plugin","qvalue","Adaptive zvalue","GFC","US","CARS","GAP")

etp=fdr
fnr=fdr
m1=k2

alpha=0.1

for (FDR in 1:N){
  alpha=fdrlevel[FDR]
  print(paste("the fdrlevel is",alpha))
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
  fdr[,FDR]=apply(fdp,1,mean)#行列要考虑一下
  fnr[,FDR]=apply(fnp,1,mean)
  etp[,FDR]=apply(e,1,mean)
}

library(xlsx)
write.xlsx(fnr,"D:/蔡健/毕设安排/数值模拟代码/数据结果/fnrcars2.xlsx")
write.xlsx(etp,"D:/蔡健/毕设安排/数值模拟代码/数据结果/etpcars2.xlsx")
write.xlsx(fdr,"D:/蔡健/毕设安排/数值模拟代码/数据结果/fdrcars2.xlsx")

