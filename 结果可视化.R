#可视化yyds
library(ggplot2)
#install.packages("Cairo")
library(Cairo);


###独立情形差异程度
fdr1=read.csv("D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/fdr1.csv")
etp1=read.csv("D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/etp1.csv")
fnr1=read.csv("D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/fnr1.csv")


rownames(fdr1) <- fdr1[, 1]
fdr1<- fdr1[, -1]
fdr1<- fdr1[-c(2,3,6),]
method=rownames(fdr1)

colnames(fdr1)<- seq(3,6,0.2)
b=NULL
for (i in 1:6){
  a=cbind(rep(rownames(fdr1)[i],length(fdr1[i,])),rownames(t(fdr1[i,])),t(fdr1[i,]))
  b=rbind(b,a)
}
fdr1=b
colnames(fdr1)<-c("group","u0","FDR")
rownames(fdr1)<-NULL
fdr1=as.data.frame(fdr1)
fdr1[,2]=as.numeric(fdr1[,2])
fdr1[,3]=as.numeric(fdr1[,3])

rownames(fnr1) <- fnr1[, 1]
fnr1<- fnr1[, -1]
fnr1<- fnr1[-c(2,3,6),]
colnames(fnr1)<- seq(3,6,0.2)
b=NULL
for (i in 1:6){
  a=cbind(rep(rownames(fnr1)[i],length(fnr1[i,])),rownames(t(fnr1[i,])),t(fnr1[i,]))
  b=rbind(b,a)
}
fnr1=b
colnames(fnr1)<-c("group","u0","FNR")
rownames(fnr1)<-NULL
fnr1=as.data.frame(fnr1)
fnr1[,2]=as.numeric(fnr1[,2])
fnr1[,3]=as.numeric(fnr1[,3])

rownames(etp1) <- etp1[, 1]
etp1<- etp1[, -1]
etp1<- etp1[-c(2,3,6),]
colnames(etp1)<- seq(3,6,0.2)
b=NULL
for (i in 1:6){
  a=cbind(rep(rownames(etp1)[i],length(etp1[i,])),rownames(t(etp1[i,])),t(etp1[i,]))
  b=rbind(b,a)
}
etp1=b
colnames(etp1)<-c("group","u0","ETP")
rownames(etp1)<-NULL
etp1=as.data.frame(etp1)
etp1[,2]=as.numeric(etp1[,2])
etp1[,3]=as.numeric(etp1[,3])

####分组折线图
#####FDR
groupnum <- length(method)
xrange <- range(fdr1$u0)
yrange <- range(fdr1$FDR)
pdf(family="GB1",file='D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/独立_差异水平_FDR.pdf' )
plot(xrange, yrange,
     type="n",
     xlab=expression(u[0]),
     ylab="FDR"
)
colors <- rainbow(groupnum)
linetype <- c(1:groupnum)
plotchar <- seq(1, 1+groupnum, 1)
for (i in 1:groupnum) {
  group <- subset(fdr1, group==method[i])
  lines(group$u0, group$FDR,
        type="b",
        lwd=2,
        lty=linetype[i],
        col=colors[i],
        pch=plotchar[i]
  )
}
title("不同差异水平下各FDR控制方法的FDR","Model 1")
legend=method
legend[2]<-"QV"
legend[3]<-"AZ"
legend("topright",cex=1,legend=legend,col=colors, pch=plotchar,lty=linetype,ncol=2)   
dev.off()

#####FNR
xrange <- range(fnr1$u0)
yrange <- range(fnr1$FNR)
pdf(family="GB1",file='D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/独立_差异水平_FNR.pdf' )
plot(xrange, yrange,
     type="n",
     xlab=expression(u[0]),
     ylab="FNR"
)
colors <- rainbow(groupnum)
linetype <- c(1:groupnum)
plotchar <- seq(1, 1+groupnum, 1)
for (i in 1:groupnum) {
  group <- subset(fnr1, group==method[i])
  lines(group$u0, group$FNR,
        type="b",
        lwd=2,
        lty=linetype[i],
        col=colors[i],
        pch=plotchar[i]
  )
}
title("不同差异水平下各FDR控制方法的FNR","Model 1")
legend=method
legend[2]<-"QV"
legend[3]<-"AZ"
legend("topright",cex=1,legend=legend,col=colors, pch=plotchar,lty=linetype,ncol=2)   
dev.off()

#####ETP
xrange <- range(etp1$u0)
yrange <- range(etp1$ETP)
yrange[2]=yrange[2]+0.2
pdf(family="GB1",file='D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/独立_差异水平_ETP.pdf' )
plot(xrange, yrange,
     type="n",
     xlab=expression(u[0]),
     ylab="ETP"
)
colors <- rainbow(groupnum)
linetype <- c(1:groupnum)
plotchar <- seq(1, 1+groupnum, 1)
for (i in 1:groupnum) {
  group <- subset(etp1, group==method[i])
  lines(group$u0, group$ETP,
        type="b",
        lwd=2,
        lty=linetype[i],
        col=colors[i],
        pch=plotchar[i]
  )
}
title("不同差异水平下各FDR控制方法的ETP","Model 1")
legend=method
legend[2]<-"QV"
legend[3]<-"AZ"
legend("topright",cex=1,legend=legend,col=colors, pch=plotchar,lty=linetype,ncol=2)   
dev.off()


###独立情形稀疏水平k2/k1
fdr2=read.csv("D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/fdr2.csv")
etp2=read.csv("D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/etp2.csv")
fnr2=read.csv("D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/fnr2.csv")


rownames(fdr2) <- fdr2[, 1]
fdr2<- fdr2[, -1]
fdr2<- fdr2[-c(2,3,6),]
method=rownames(fdr2)

colnames(fdr2)<- seq(100,1000,100)
b=NULL
for (i in 1:6){
  a=cbind(rep(rownames(fdr2)[i],length(fdr2[i,])),rownames(t(fdr2[i,])),t(fdr2[i,]))
  b=rbind(b,a)
}
fdr2=b
colnames(fdr2)<-c("group","k2","FDR")
rownames(fdr2)<-NULL
fdr2=as.data.frame(fdr2)
fdr2[,2]=as.numeric(fdr2[,2])
fdr2[,3]=as.numeric(fdr2[,3])

rownames(fnr2) <- fnr2[, 1]
fnr2<- fnr2[, -1]
fnr2<- fnr2[-c(2,3,6),]
colnames(fnr2)<- seq(100,1000,100)
b=NULL
for (i in 1:6){
  a=cbind(rep(rownames(fnr2)[i],length(fnr2[i,])),rownames(t(fnr2[i,])),t(fnr2[i,]))
  b=rbind(b,a)
}
fnr2=b
colnames(fnr2)<-c("group","k2","FNR")
rownames(fnr2)<-NULL
fnr2=as.data.frame(fnr2)
fnr2[,2]=as.numeric(fnr2[,2])
fnr2[,3]=as.numeric(fnr2[,3])

rownames(etp2) <- etp2[, 1]
etp2<- etp2[, -1]
etp2<- etp2[-c(2,3,6),]
colnames(etp2)<- seq(100,1000,100)
b=NULL
for (i in 1:6){
  a=cbind(rep(rownames(etp2)[i],length(etp2[i,])),rownames(t(etp2[i,])),t(etp2[i,]))
  b=rbind(b,a)
}
etp2=b
colnames(etp2)<-c("group","k2","ETP")
rownames(etp2)<-NULL
etp2=as.data.frame(etp2)
etp2[,2]=as.numeric(etp2[,2])
etp2[,3]=as.numeric(etp2[,3])

####分组折线图
#####FDR
groupnum <- length(method)
xrange <- range(fdr2$k2)
yrange <- range(fdr2$FDR)
pdf(family="GB1",file='D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/独立_稀疏水平_FDR.pdf' )
plot(xrange, yrange,
     type="n",
     xlab=expression(k[2]),
     ylab="FDR"
)
colors <- rainbow(groupnum)
linetype <- c(1:groupnum)
plotchar <- seq(1, 1+groupnum, 1)
for (i in 1:groupnum) {
  group <- subset(fdr2, group==method[i])
  lines(group$k2, group$FDR,
        type="b",
        lwd=2,
        lty=linetype[i],
        col=colors[i],
        pch=plotchar[i]
  )
}
title("不同稀疏水平下各FDR控制方法的FDR","Model 2")
legend=method
legend[2]<-"QV"
legend[3]<-"AZ"
legend("topright",cex=1,legend=legend,col=colors, pch=plotchar,lty=linetype,ncol=2)   
dev.off()

#####FNR
xrange <- range(fnr2$k2)
yrange <- range(fnr2$FNR)
pdf(family="GB1",file='D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/独立_稀疏水平_FNR.pdf' )
plot(xrange, yrange,
     type="n",
     xlab=expression(k[2]),
     ylab="FNR"
)
colors <- rainbow(groupnum)
linetype <- c(1:groupnum)
plotchar <- seq(1, 1+groupnum, 1)
for (i in 1:groupnum) {
  group <- subset(fnr2, group==method[i])
  lines(group$k2, group$FNR,
        type="b",
        lwd=2,
        lty=linetype[i],
        col=colors[i],
        pch=plotchar[i]
  )
}
title("不同稀疏水平下各FDR控制方法的FNR","Model 2")
legend=method
legend[2]<-"QV"
legend[3]<-"AZ"
legend("topright",cex=1,legend=legend,col=colors, pch=plotchar,lty=linetype,ncol=2)   
dev.off()

#####ETP
xrange <- range(etp2$k2)
yrange <- range(etp2$ETP)
yrange[2]=yrange[2]+0.1
pdf(family="GB1",file='D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/独立_稀疏水平_ETP.pdf' )
plot(xrange, yrange,
     type="n",
     xlab=expression(k[2]),
     ylab="ETP"
)
colors <- rainbow(groupnum)
linetype <- c(1:groupnum)
plotchar <- seq(1, 1+groupnum, 1)
for (i in 1:groupnum) {
  group <- subset(etp2, group==method[i])
  lines(group$k2, group$ETP,
        type="b",
        lwd=2,
        lty=linetype[i],
        col=colors[i],
        pch=plotchar[i]
  )
}
title("不同稀疏水平下各FDR控制方法的ETP","Model 2")
legend=method
legend[2]<-"QV"
legend[3]<-"AZ"
legend("topright",cex=1,legend=legend,col=colors, pch=plotchar,lty=linetype,ncol=2)   
dev.off()


###相依情形差异程度beta
fdr3=read.csv("D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/fdr3.csv")
etp3=read.csv("D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/etp3.csv")
fnr3=read.csv("D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/fnr3.csv")


rownames(fdr3) <- fdr3[, 1]
fdr3<- fdr3[, -1]
fdr3<- fdr3[-c(3),]
method=rownames(fdr3)
N=length(method)

colnames(fdr3)<- seq(3.5,5.5,0.2)
b=NULL
for (i in 1:N){
  a=cbind(rep(rownames(fdr3)[i],length(fdr3[i,])),rownames(t(fdr3[i,])),t(fdr3[i,]))
  b=rbind(b,a)
}
fdr3=b
colnames(fdr3)<-c("group","beta","FDR")
rownames(fdr3)<-NULL
fdr3=as.data.frame(fdr3)
fdr3[,2]=as.numeric(fdr3[,2])
fdr3[,3]=as.numeric(fdr3[,3])

rownames(fnr3) <- fnr3[, 1]
fnr3<- fnr3[, -1]
fnr3<- fnr3[-c(3),]
colnames(fnr3)<- seq(3.5,5.5,0.2)
b=NULL
for (i in 1:N){
  a=cbind(rep(rownames(fnr3)[i],length(fnr3[i,])),rownames(t(fnr3[i,])),t(fnr3[i,]))
  b=rbind(b,a)
}
fnr3=b
colnames(fnr3)<-c("group","beta","FNR")
rownames(fnr3)<-NULL
fnr3=as.data.frame(fnr3)
fnr3[,2]=as.numeric(fnr3[,2])
fnr3[,3]=as.numeric(fnr3[,3])

rownames(etp3) <- etp3[, 1]
etp3<- etp3[, -1]
etp3<- etp3[-c(3),]
colnames(etp3)<- seq(3.5,5.5,0.2)
b=NULL
for (i in 1:N){
  a=cbind(rep(rownames(etp3)[i],length(etp3[i,])),rownames(t(etp3[i,])),t(etp3[i,]))
  b=rbind(b,a)
}
etp3=b
colnames(etp3)<-c("group","beta","ETP")
rownames(etp3)<-NULL
etp3=as.data.frame(etp3)
etp3[,2]=as.numeric(etp3[,2])
etp3[,3]=as.numeric(etp3[,3])

####分组折线图
fdr3=fdr3[-c(72:77),]
fnr3=fnr3[-c(72:77),]
etp3=etp3[-c(72:77),]
fdr3=subset(fdr3,beta<4.5)
fnr3=subset(fnr3,beta<4.5)
etp3=subset(etp3,beta<4.5)
#####FDR
groupnum <- length(method)
xrange <- range(fdr3$beta)
yrange <- range(fdr3$FDR)
pdf(family="GB1",file='D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/相依_差异程度_FDR.pdf' )
plot(xrange, yrange,
     type="n",
     xlab=expression(beta),
     ylab="FDR"
)
colors <- rainbow(groupnum)
linetype <- c(1:groupnum)
plotchar <- seq(1, 1+groupnum, 1)
for (i in 1:groupnum) {
  group <- subset(fdr3, group==method[i])
  lines(group$beta, group$FDR,
        type="b",
        lwd=2,
        lty=linetype[i],
        col=colors[i],
        pch=plotchar[i]
  )
}
title("不同差异程度下各FDR控制方法的FDR","Model 3")
legend=method
legend[2]<-"ABH"
legend[3]<-"QV"
legend[4]<-"AZ"
legend("topright",cex=1,legend=legend,col=colors, pch=plotchar,lty=linetype,ncol=2)   
dev.off()

#####FNR
xrange <- range(fnr3$beta)
yrange <- range(fnr3$FNR)
pdf(family="GB1",file='D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/相依_差异程度_FNR.pdf' )
plot(xrange, yrange,
     type="n",
     xlab=expression(beta),
     ylab="FNR"
)
colors <- rainbow(groupnum)
linetype <- c(1:groupnum)
plotchar <- seq(1, 1+groupnum, 1)
for (i in 1:groupnum) {
  group <- subset(fnr3, group==method[i])
  lines(group$beta, group$FNR,
        type="b",
        lwd=2,
        lty=linetype[i],
        col=colors[i],
        pch=plotchar[i]
  )
}
title("不同差异程度下各FDR控制方法的FNR","Model 3")
legend=method
legend[2]<-"ABH"
legend[3]<-"QV"
legend[4]<-"AZ"
legend("topright",cex=1,legend=legend,col=colors, pch=plotchar,lty=linetype,ncol=2)   
dev.off()

#####ETP
xrange <- range(etp3$beta)
yrange <- range(etp3$ETP)
yrange[2]=yrange[2]+0.3
pdf(family="GB1",file='D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/相依_差异程度_ETP.pdf' )
plot(xrange, yrange,
     type="n",
     xlab=expression(beta),
     ylab="ETP"
)
colors <- rainbow(groupnum)
linetype <- c(1:groupnum)
plotchar <- seq(1, 1+groupnum, 1)
for (i in 1:groupnum) {
  group <- subset(etp3, group==method[i])
  lines(group$beta, group$ETP,
        type="b",
        lwd=2,
        lty=linetype[i],
        col=colors[i],
        pch=plotchar[i]
  )
}
title("不同差异程度下各FDR控制方法的ETP","Model 3")
legend=method
legend[2]<-"ABH"
legend[3]<-"QV"
legend[4]<-"AZ"
legend("topright",cex=1,legend=legend,col=colors, pch=plotchar,lty=linetype,ncol=2)   
dev.off()


###相依情形稀疏程度k
fdr4=read.csv("D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/fdr4.csv")
etp4=read.csv("D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/etp4.csv")
fnr4=read.csv("D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/fnr4.csv")


rownames(fdr4) <- fdr4[, 1]
fdr4<- fdr4[, -1]
fdr4<- fdr4[-c(3),]
method=rownames(fdr4)
N=length(method)

colnames(fdr4)<- seq(50,150,10)
b=NULL
for (i in 1:N){
  a=cbind(rep(rownames(fdr4)[i],length(fdr4[i,])),rownames(t(fdr4[i,])),t(fdr4[i,]))
  b=rbind(b,a)
}
fdr4=b
colnames(fdr4)<-c("group","k","FDR")
rownames(fdr4)<-NULL
fdr4=as.data.frame(fdr4)
fdr4[,2]=as.numeric(fdr4[,2])
fdr4[,3]=as.numeric(fdr4[,3])

rownames(fnr4) <- fnr4[, 1]
fnr4<- fnr4[, -1]
fnr4<- fnr4[-c(3),]
colnames(fnr4)<- seq(50,150,10)
b=NULL
for (i in 1:N){
  a=cbind(rep(rownames(fnr4)[i],length(fnr4[i,])),rownames(t(fnr4[i,])),t(fnr4[i,]))
  b=rbind(b,a)
}
fnr4=b
colnames(fnr4)<-c("group","k","FNR")
rownames(fnr4)<-NULL
fnr4=as.data.frame(fnr4)
fnr4[,2]=as.numeric(fnr4[,2])
fnr4[,3]=as.numeric(fnr4[,3])

rownames(etp4) <- etp4[, 1]
etp4<- etp4[, -1]
etp4<- etp4[-c(3),]
colnames(etp4)<- seq(50,150,10)
b=NULL
for (i in 1:N){
  a=cbind(rep(rownames(etp4)[i],length(etp4[i,])),rownames(t(etp4[i,])),t(etp4[i,]))
  b=rbind(b,a)
}
etp4=b
colnames(etp4)<-c("group","k","ETP")
rownames(etp4)<-NULL
etp4=as.data.frame(etp4)
etp4[,2]=as.numeric(etp4[,2])
etp4[,3]=as.numeric(etp4[,3])

####分组折线图
#####FDR
groupnum <- length(method)
xrange <- range(fdr4$k)
yrange <- range(fdr4$FDR)
yrange[2]=yrange[2]+0.02
pdf(family="GB1",file='D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/相依_稀疏水平_FDR.pdf' )
plot(xrange, yrange,
     type="n",
     xlab=expression(k[2]),
     ylab="FDR"
)
colors <- rainbow(groupnum)
linetype <- c(1:groupnum)
plotchar <- seq(1, 1+groupnum, 1)
for (i in 1:groupnum) {
  group <- subset(fdr4, group==method[i])
  lines(group$k, group$FDR,
        type="b",
        lwd=2,
        lty=linetype[i],
        col=colors[i],
        pch=plotchar[i]
  )
}
title("不同稀疏水平下各FDR控制方法的FDR","Model 4")
legend=method
legend[2]<-"ABH"
legend[3]<-"QV"
legend[4]<-"AZ"
legend("topright",cex=1,legend=legend,col=colors, pch=plotchar,lty=linetype,ncol=2)   
dev.off()

#####FNR
xrange <- range(fnr4$k)
yrange <- range(fnr4$FNR)
yrange[2]=yrange[2]+0.015
pdf(family="GB1",file='D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/相依_稀疏水平_FNR.pdf' )
plot(xrange, yrange,
     type="n",
     xlab=expression(k[2]),
     ylab="FNR"
)
colors <- rainbow(groupnum)
linetype <- c(1:groupnum)
plotchar <- seq(1, 1+groupnum, 1)
for (i in 1:groupnum) {
  group <- subset(fnr4, group==method[i])
  lines(group$k, group$FNR,
        type="b",
        lwd=2,
        lty=linetype[i],
        col=colors[i],
        pch=plotchar[i]
  )
}
title("不同稀疏水平下各FDR控制方法的FNR","Model 4")
legend=method
legend[2]<-"ABH"
legend[3]<-"QV"
legend[4]<-"AZ"
legend("topright",cex=1,legend=legend,col=colors, pch=plotchar,lty=linetype,ncol=2)   
dev.off()

#####ETP
xrange <- range(etp4$k)
yrange <- range(etp4$ETP)
yrange[2]=yrange[2]+0.1
pdf(family="GB1",file='D:/蔡健/毕设安排/数值模拟代码/数据结果/结果/结果/相依_稀疏水平_ETP.pdf' )
plot(xrange, yrange,
     type="n",
     xlab=expression(k[2]),
     ylab="ETP"
)
colors <- rainbow(groupnum)
linetype <- c(1:groupnum)
plotchar <- seq(1, 1+groupnum, 1)
for (i in 1:groupnum) {
  group <- subset(etp4, group==method[i])
  lines(group$k, group$ETP,
        type="b",
        lwd=2,
        lty=linetype[i],
        col=colors[i],
        pch=plotchar[i]
  )
}
title("不同稀疏水平下各FDR控制方法的ETP","Model 4")
legend=method
legend[2]<-"ABH"
legend[3]<-"QV"
legend[4]<-"AZ"
legend("topright",cex=1,legend=legend,col=colors, pch=plotchar,lty=linetype,ncol=2)   
dev.off()


