#模拟数据的生成
#H0 和 H1的区别要弄大一点，mu1 和mu2 差别大就对，差别小随机因素太大了


##Model 1 exactly  sparse
n1=100#样本容量
n2=100#样本容量
m=2000#变量数/检验数
rep=500#重复试验次数
m1=floor(sqrt(m))#nonnull的个数
alpha=seq(0.05,1,length.out = 20)#预设FDR水平
mu1=c(rep(3*sqrt(log(m)/n1),m1),rep(0,m-m1))
mu2=c(rep(2*sqrt(log(m)/n2),m1),rep(0,m-m1))
#真实标签
true=c(rep(1,m1),rep(0,m-m1))


#XY比较同方差情形
X=matrix(mu1,m,n1)+matrix(rnorm(m*n1,0,1),m,n1)#生成X样本
Y=matrix(mu2,m,n2)+matrix(rnorm(m*n2,0,1),m,n2)#生成Y样本
#XZ比较异方差情形
X=matrix(mu1,m,n1)+matrix(rnorm(m*n1,0,1),m,n1)#生成X样本
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

#计算p值
p1=(1-pt(T1,n1+n2-2))*2
p2=(1-pnorm(T1,0,1))*2#异方差是渐进正态统计量，
#关于p值和G，z值的关系可以看GAP supp第九页或者oracle comupoud decision的正文


##Model 2 exactly  sparse
n1=100#样本容量
n2=100#样本容量
m=2000#变量数/检验数
rep=500#重复试验次数
m1=floor(sqrt(m))#nonnull的个数
alpha=seq(0.05,1,length.out = 20)#预设FDR水平
mu1=c(rep(2*sqrt(log(m)/n1),m1),rep(0,m-m1))
mu2=c(rep(1*sqrt(log(m)/n2),floor(m1/2)),rep(-0.5*sqrt(log(m)/n2),m1-floor(m1/2)),rep(0,m-m1))
#XY比较同方差情形
X=matrix(mu1,m,n1)+matrix(rnorm(m*n1,0,1),m,n1)#生成X样本
Y=matrix(mu2,m,n2)+matrix(rnorm(m*n2,0,1),m,n2)#生成Y样本
#XZ比较异方差情形
X=matrix(mu1,m,n1)+matrix(rnorm(m*n1,0,1),m,n1)#生成X样本
Z=matrix(mu2,m,n2)+matrix(rnorm(m*n2,0,0.5),m,n2)#生成Z样本


##Model 3 asymptotically sparse
n1=100#样本容量
n2=100#样本容量
m=2000#变量数/检验数
rep=500#重复试验次数
m1=floor(sqrt(m))#nonnull的个数
alpha=seq(0.05,1,length.out = 20)#预设FDR水平

#设置均值向量mu1和mu2
mu1=0#初始化
mu2=0#初始化
for (i in 1:m){
  if (i <= m1){
    mu1[i]=3*sqrt(log(m)/n1)
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


##Model 4 non-sparse
n1=100#样本容量
n2=100#样本容量
m=2000#变量数/检验数
rep=500#重复试验次数
m1=floor(sqrt(m))#nonnull的个数
alpha=seq(0.05,1,length.out = 20)#预设FDR水平

#生成均值向量mu1和mu2
mu1=0#初始化
mu2=0#初始化
for (i in 1:m){
  if (i<=m1){
    mu1[i]=3*sqrt(log(m)/n1)
    mu2[i]=2*sqrt(log(m)/n2)
  }
  else if (i <=m1+sqrt(m)){
    mu1[i]=1
    mu2[i]=1
  }
  else {
    mu1[i]=0.2
    mu2[i]=0.2
  }
}

#XY比较同方差情形
X=matrix(mu1,m,n1)+matrix(rnorm(m*n1,0,1),m,n1)#生成X样本
Y=matrix(mu2,m,n2)+matrix(rnorm(m*n2,0,1),m,n2)#生成Y样本
#XZ比较异方差情形
X=matrix(mu1,m,n1)+matrix(rnorm(m*n1,0,1),m,n1)#生成X样本
Z=matrix(mu2,m,n2)+matrix(rnorm(m*n2,0,0.5),m,n2)#生成Z样本


#Model 5 CARS setting 1
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

#Model 6 CARS setting 2
n1=50
n2=60
varx=1
vary=4
m=5000
alpha=0.05
rep=500

k1=2000#k1表示非零元的下标终点
k2=100#k2=seq(100,1000,100) 表示差别元的下标终点
mu1=0
mu2=0
for (i in 1:m){
  if (i %in% 1:k2){
    mu1[i]=5/sqrt(30)
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

X=matrix(mu1,m,n1)+matrix(rnorm(m*n1,0,sdx),m,n1)#生成X样本
Y=matrix(mu2,m,n2)+matrix(rnorm(m*n2,0,sdy),m,n2)#生成Y样本


#Model 7 CARS setting 3
n1=50
n2=60
varx=1
vary=4
m=5000
alpha=0.05
rep=500
mu0=3.5#mu0=3.5:5
k=750
mu1=0
mu2=0
for (i in 1:m){
  if (i %in% 1:k){
    mu1[i]=mu0/sqrt(30)
    mu2[i]=1/sqrt(30)
  }
  else if (i %in% (k+1):(2*k)){
    mu1[i]=3/sqrt(30)
    mu2[i]=3/sqrt(30)
  }
  else {
    mu1[i]=0
    mu2[i]=0
  }
}

X=matrix(mu1,m,n1)+matrix(rnorm(m*n1,0,sdx),m,n1)#生成X样本
Y=matrix(mu2,m,n2)+matrix(rnorm(m*n2,0,sdy),m,n2)#生成Y样本


#Model 8 GAP setting 1 
#以上模型可能是独立的检验，所以可以任意换下标，以下考虑相关情形
#这里样本量不明，要多试验
n1=80
n2=100
p=2000
rep=200
#生成多元正态分布随机向量，每行一个观测
library(MASS)
Sigma=matrix(NA,p,p)
for (i in 1:p){
  for (j in 1:p){
    Sigma[i,j]=0.8^(abs(i-j))
  }
}

beta=6.5#beta=6.5:9
#猜测后面还有1900个相同的，因为H1要稀疏，文中p=2000,GAP要求s1和s2稀疏
mu1=c(rep(3,50),rep(-3,50),rep(0,1900))
mu2=c(rep(beta,50),rep(-beta,50),rep(0,1900))

X=mvrnorm(n=n1, mu1, Sigma)#生成X样本
Y=mvrnorm(n=n2, mu2, Sigma)#生成Y样本

#Model 9 GAP setting 2
n1=80
n2=100
p=2000
rep=200
k0=floor(p/3)
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




beta=6.5#beta=6.5:9
mu1=c(rep(3,50),rep(-3,50),rep(0,1900))
mu2=c(rep(beta,50),rep(-beta,50),rep(0,1900))

X=mvrnorm(n=n1, mu1, Sigma)#生成X样本
Y=mvrnorm(n=n2, mu2, Sigma)#生成Y样本

#Model 10 GAP setting 3
n1=80
n2=100
p=2000
rep=200
#生成多元正态分布随机向量，每行一个观测
library(MASS)
Sigma=diag(1,p)#省去对 对角元的赋值
for (i in 1:p){
  for (j in 1:p){
    if (i<j) {
      Sigma[i,j]=0.5*rbinom(1,1,0.05)
    }
  }
}
for (i in 1:p){
  for (j in 1:p){
    if (i>j) {
      Sigma[i,j]=Sigma[j,i]
    }
  }
}
Delta=eigen(Sigma)
delta=abs(min(Delta$values))+0.05
Sigma=(Sigma+diag(delta,p))/(1+delta)

beta=6.5#beta=6.5:9
#猜测后面还有1900个相同的，因为H1要稀疏，文中p=2000,GAP要求s1和s2稀疏
mu1=c(rep(3,50),rep(-3,50),rep(0,1900))
mu2=c(rep(beta,50),rep(-beta,50),rep(0,1900))

X=mvrnorm(n=n1, mu1, Sigma)#生成X样本
Y=mvrnorm(n=n2, mu2, Sigma)#生成Y样本
