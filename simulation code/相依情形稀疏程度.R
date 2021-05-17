#Model 8 GAP setting 1 
#以上模型可能是独立的检验，所以可以任意换下标，以下考虑相关情形
#这里样本量不明，要多试验


#GAP的setting 2 按下面的数据去跑可以说明GAP优于BH等方法，但是US更好，这还是个迷


#fdr控制效果还好，要说明比US，BH好会更妙


n1=100
n2=120
m=2000
rep=200
k0=floor(p/3)#生成多元正态分布随机向量，每行一个观测
library(MASS)
library(qvalue)
library(fdrtool)
library(CARS)
Sigma=matrix(NA,m,m)
for (i in 1:m){
  for (j in 1:m){
    Sigma[i,j]=0.8^(abs(i-j))
  }
}
beta=4
p=m
library(MASS)
Sigma=matrix(0,p,p)#之后赋值的时候就不用另外赋值0了，因为Sigma0不一定能到p阶方阵！

#形成分块对角阵
list2 <- NULL

for (i in 1:k0){
  list2[[i]] <- matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),3,3)
}

library(Matrix)

Sigma0 <- as.matrix(bdiag(list2))
Sigma[1:(3*k0),1:(3*k0)] =Sigma0#不能有零，要把剩下的补齐，否则会有缺失值
changdu=3*k0+1
Sigma[changdu:m,changdu:m]=diag(0.8,m-3*k0)
alpha=0.1
rep=100
K=seq(50,150,10)
B=seq(4,6,0.5)
N=length(K)
fdr=matrix(0,9,N)
rownames(fdr)=c("BH","Adaptive BH","Plugin","qvalue","Adaptive zvalue","GFC","US","CARS","GAP")
etp=fdr
fnr=fdr
ppp=50
beta=3.6
for (j in 1:N){
  ppp=K[j]#beta=6.5:9
  #猜测后面还有1900个相同的，因为H1要稀疏，文中p=2000,GAP要求s1和s2稀疏
  mu1=c(rep(3,ppp),rep(-3,ppp),rep(0,m-2*ppp))
  mu2=c(rep(beta,ppp),rep(-beta,ppp),rep(0,m-2*ppp))
  print(paste("iter_diff_level=",beta))
  m1=2*ppp
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
  for (i in 1:rep){
    ##数据生成
    print(paste("the iteration is",i))
    #XY比较同方差情形
    X=t(mvrnorm(n=n1, mu1, Sigma))#生成X样本
    Y=t(mvrnorm(n=n2, mu2, Sigma))#生成Y样本
    
    Xbar=apply(X,1,mean)
    Xvar=apply(X,1,var)
    Ybar=apply(Y,1,mean)
    Yvar=apply(Y,1,var)
    T1=1/sqrt(Xvar/n1+Yvar/n2)*(Xbar-Ybar)#异方差模型XZ检验统计量
    S1=sqrt(n1/(Xvar*(1+n2/n1*Xvar/Yvar)))*(Xbar+n2/n1*Xvar/Yvar*Ybar)
    p1=(1-pnorm(abs(T1),0,1))*2#异方差是渐进正态统计量，
    #关于p值和G，z值的关系可以看GAP supp第九页或者oracle comupoud decision的正文
    
    ##计算fdp等指标
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
  fdr[,j]=apply(fdp,1,mean)#行列要考虑一下
  fnr[,j]=apply(fnp,1,mean)
  etp[,j]=apply(e,1,mean)
}
write.csv(fnr,"C:/Users/Lenovo/Desktop/结果/fnr4.csv")
write.csv(etp,"C:/Users/Lenovo/Desktop/结果/etp4.csv")
write.csv(fdr,"C:/Users/Lenovo/Desktop/结果/fdr4.csv")


