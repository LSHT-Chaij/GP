#US procedrue
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

n1=100#样本容量
n2=100#样本容量
m=2000#变量数/检验数
rep=500#重复试验次数
m1=500#nonnull的个数
alpha=seq(0.05,1,length.out = 20)#预设FDR水平
mu1=c(rep(4*sqrt(log(m)/n1),m1),rep(0,m-m1))
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

flag=usprocedure(T1,S1,0.1,10)
flag=rbind(flag,true)
falsediscovery=sum(flag[1,]==1 & flag[2,]==0)
rejections=sum(flag[1,]==1)
fdp=falsediscovery/max(rejections,1)
print(rejections)
print(fdp)
fdr=fdr+fdp
}
fdr=fdr/rep
fdr

