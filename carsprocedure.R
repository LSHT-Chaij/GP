#CARS procedure
install.packages("CARS")
library(CARS)
CARS(X,Y,0.1,tau=0.1,option='regular')
CARS(X,Y,0.1,tau=0.5,option='sparse')#损失power

#XY是数据矩阵
n1=50
n2=60
sdx=1
sdy=2
m=5000
alpha=0.05
rep=500

k=100#k=seq(100,1000,100) 表示稀疏程度
mu1=0
mu2=0
for (i in 1:m){
  if (i %in% 1:k){
    mu1[i]=5/sqrt(30)
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

X=matrix(mu1,m,n1)+matrix(rnorm(m*n1,0,sdx),m,n1)#生成X样本
Y=matrix(mu2,m,n2)+matrix(rnorm(m*n2,0,sdy),m,n2)#生成Y样本
#真实标签
flag=CARS(X,Y,0.05,tau=0.1,option='regular')
true=c(rep(1,k),rep(0,m-k))
flag=flag$de
flag=rbind(flag,true)
falsediscovery=sum(flag[1,]==1 & flag[2,]==0)
rejections=sum(flag[1,]==1)
falsediscovery/rejections
