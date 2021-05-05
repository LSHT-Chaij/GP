adaptZ.func<-function(zv, q)
{
  # the input
  # zv is the z-values transformed from m tests
  # q is the desired FDR level
  # the output is a list with
  # the first element (st.lfdr) the sorted local fdr values
  # the second element (k) the number of hypotheses to be rejected
  # the third element (lfdrk) the threshold for the local fdr values
  # the fourth element (reject) the set of indices of the rejected hypotheses
  # the fifth element (accept) the set of indices of the accepted hypotheses
  ## the estimates for the local fdr statistics
  # density estimates 
  zv.ds<-density(zv, from=min(zv)-10, to=max(zv)+10, n=1000)
  # linear interpolation
  zv.ds<-lin.itp(zv, zv.ds$x, zv.ds$y)
  # estimating the null distribution
  zv.MuSigma<-EstNull.func(zv)
  mu<-zv.MuSigma$mu
  s<-zv.MuSigma$s
  # mu<-0; s<-1
  zv.p0<-1-epsest.func(zv, mu, s)
  zv.lfdr<-zv.p0*dnorm(zv, mu, s)/zv.ds
  y<-adpt.cutz(zv.lfdr, q)
  return (y)
}

adpt.cutz<-function(lfdr, q)
{
  # the input
  # lfdr the vector of local fdr statistics
  # q the desired FDR level
  # the output is a list with
  # the first element (st.lfdr) the sorted local fdr values
  # the second element (k) the number of hypotheses to be rejected
  # the third element (lfdrk) the threshold for the local fdr values
  # the fourth element (reject) the set of indices of the rejected hypotheses
  # the fifth element (accept) the set of indices of the accepted hypotheses
  
  m=length(lfdr)
  st.lfdr<-sort(lfdr)
  k=1
  while(k<m && (1/k)*sum(st.lfdr[1:k])<q){
    k=k+1
  }
  k<-k-1
  lfdrk<-st.lfdr[k]
  reject<-which(lfdr<=lfdrk)
  accept<-which(lfdr>lfdrk)
  y<-list(sf=st.lfdr, nr=k, thr=lfdrk, re=reject, ac=accept)
  return (y)
}

epsest.func <- function(x,u,sigma)
{
  # x is a vector
  # u is the mean
  # sigma is the standard deviation
  
  z  = (x - u)/sigma
  xi = c(0:100)/100
  tmax=sqrt(log(length(x)))
  tt=seq(0,tmax,0.1)
  
  epsest=NULL
  
  for (j in 1:length(tt)) { 
    
    t=tt[j]
    f  = t*xi
    f  = exp(f^2/2)
    w  = (1 - abs(xi))
    co  = 0*xi
    
    for (i in 1:101) {
      co[i] = mean(cos(t*xi[i]*z));
    } 
    epshat = 1 - sum(w*f*co)/sum(w)
    epsest=c(epsest,epshat)
  }
  return(epsest=max(epsest))
}

EstNull.func<-function (x,gamma=0.1)
{
  # x is a vector of z-values
  # gamma is a parameter, default is 0.1
  # output the estimated mean and standard deviation
  
  n = length(x)
  t = c(1:1000)/200
  
  gan    = n^(-gamma)
  that   = 0 
  shat   = 0
  uhat   = 0
  epshat = 0
  
  phiplus   = rep(1,1000)
  phiminus  = rep(1,1000)
  dphiplus  = rep(1,1000)
  dphiminus = rep(1,1000)
  phi       = rep(1,1000)
  dphi      = rep(1,1000)
  
  for (i in 1:1000) {
    s = t[i]
    phiplus[i]   = mean(cos(s*x))
    phiminus[i]  = mean(sin(s*x))
    dphiplus[i]  = -mean(x*sin(s*x))
    dphiminus[i] = mean(x*cos(s*x))
    phi[i]       = sqrt(phiplus[i]^2 + phiminus[i]^2)
  }
  
  ind = min(c(1:1000)[(phi - gan) <= 0])
  tt = t[ind]
  a  = phiplus[ind]
  b  = phiminus[ind]
  da = dphiplus[ind]
  db = dphiminus[ind]
  c  = phi[ind]
  
  that   = tt
  shat   = -(a*da + b*db)/(tt*c*c)
  shat   = sqrt(shat) 
  uhat   = -(da*b - db*a)/(c*c)
  epshat = 1 - c*exp((tt*shat)^2/2)
  
  return(musigma=list(mu=uhat,s=shat))
}

lin.itp<-function(x, X, Y){
  ## x: the coordinates of points where the density needs to be interpolated
  ## X: the coordinates of the estimated densities
  ## Y: the values of the estimated densities
  ## the output is the interpolated densities
  x.N<-length(x)
  X.N<-length(X)
  y<-rep(0, x.N)
  for (k in 1:x.N){
    i<-max(which((x[k]-X)>=0))
    if (i<X.N)
      y[k]<-Y[i]+(Y[i+1]-Y[i])/(X[i+1]-X[i])*(x[k]-X[i])
    else 
      y[k]<-Y[i]
  }
  return(y)
}

n1=100#样本容量
n2=100#样本容量
m=2000#变量数/检验数
rep=500#重复试验次数
m1=floor(sqrt(m))#nonnull的个数
m1=44
alpha=seq(0.05,1,length.out = 20)#预设FDR水平
mu1=c(rep(4*sqrt(log(m)/n1),m1),rep(0,m-m1))
mu2=c(rep(2*sqrt(log(m)/n2),m1),rep(0,m-m1))

#采用zvalues方法
true=c(rep(1,m1),rep(0,m-m1))
fdr=0
alpha=0.1
rep=500
FDR<-0.1
for(i in 1:rep){
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
  
   
  
  flag=rep(0,1300)
  a=c(rnorm(1000),rnorm(300,40,1))
  b=fdrtool(a)
  q2=adazprocedure(p1,0.1,"pvalue")
  
  
  adaptiveZ<-adaptZ.func(a, FDR)
  flag[adaptiveZ$re]=1
  falsediscovery=sum(flag[1:1000])
  fdp=falsediscovery/max(adaptiveZ$nr,1)
  fdr=fdr+fdp
  print(fdp)
}
fdr=fdr/218
fdr
# the threshold
threshold<-adaptiveZ$th
threshold
# number of rejected hypotheses
k<-adaptiveZ$nr
k
# the rejected hypotheses
rh<-adaptiveZ$re
rh



a<-fdrtool(T1,"normal")
b<-a$lfdr
hist(a$lfdr)
summary(a$lfdr)
cdf<-ecdf(a$lfdr)
plot(cdf)
length(a$lfdr[a$lfdr<1])
