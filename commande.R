source("MCMCestimation.R")
test <- read.table("cleandata.dta",header=TRUE)
test <- as.matrix(test)
###we start now 2nd January OO:OO
vari <- apply(test,1,var)
plot(sort(vari),type="h")
ind0 <- which(vari==0) ##no variation
ind <- c(ind0,which(vari>1))
test2 <- test[-ind,]
###count the number of 0
nbzero <- function(x){sum(x==0)}
nb0 <- apply(test2,1,nbzero)
plot(sort(nb0),type="h")
ind0 <- which(nb0==0)
test2 <- test2[ind0,]
###count the number of 0 of diff
nbzero <- function(x){sum(diff(x)==0)}
nb0 <- apply(test2,1,nbzero)
plot(sort(nb0),type="h")
ind0 <- which(nb0>1500)
test2 <- test2[-ind0,]
###count the number of 0 of diff lag 2 
nbzero <- function(x){sum(diff(x,lag=2)==0)}
nb0 <- apply(test2,1,nbzero)
plot(sort(nb0),type="h")
ind0 <- which(nb0>1500)
test2 <- test2[-ind0,]

###########
n.A <- nrow(test2)
set.seed(8835)
# make synthetic data (aggregate)
# sample
smp <- sort(sample( 1:n.A, 400 ))
DR1 <- test2[smp,] ##the test2
DR2 <- test2[smp,] ##the ones with NA
IL <- test2[-smp,] ##the data mart
day <- seq(1,ncol(test2),by=96)
## we could sample idx <- sample(colnames(test2),10)
start <- sort(sample(1:ncol(test2),10))
period <- unlist(lapply(start,FUN=function(x) x:(x+9)))
DR2[,period] <- 0
matplot(t(DR2[1:10,1:96]),type="l",lwd=3,lty=1,main="some consumption Januray 2nd",xlab="times",ylab="consumption")
dd <- 16
pdf("exofconsumptioncurves.pdf")
matplot(t(DR2[1:10,day[dd]:(day[dd]+96)]),type="l",lwd=3,lty=1,main="some consumption Januray 17th",xlab="times",ylab="consumption")
dev.off()
dd <- 2
matplot(t(DR2[1:10,day[dd]:(day[dd]+96)]),type="l",lwd=3,lty=1,main="some consumption Januray 3rd",xlab="times",ylab="consumption")


#  aggregate
aggDRsum <- apply(DR2,2,sum)
aggDRmean <- apply(DR2,2,mean)
pdf("baseline.pdf")
plot(aggDRsum[1:192],type="l",xlab="time",ylab="consumption",main="2nd and 3rd january")
abline(v=97)
dev.off()
pdf("baseline.pdf")
plot(aggDRsum[1:192],type="l",xlab="time",ylab="consumption",main="2nd and 3rd january")
abline(v=97)
dev.off()

pdf("DR_gp.pdf")
matplot(t(DR2[1:8,1:192]),type="l",lwd=1,lty=1,main="some consumptions Januray 2nd and 3rd",xlab="times",ylab="consumption",col=2:9)
lines(aggDRmean[1:192],lwd=4)
abline(v=97)
dev.off()

pdf("DR_gp.pdf")
matplot(t(DR2[1:8,1:192]),type="l",lwd=1,lty=1,main="Comsommation et coupure",xlab="temps",ylab="consommation",col=2:9)
lines(aggDRmean[1:192],lwd=4)
abline(v=97)
dev.off()

###DR period1
dr1 <- start[1]:(start[1]+9)
period1 <- 1:(dr1[1]-1)
estimation1 <- period1[1]:(dr1[10]+10)
#########################################
MCMCperiod1 <- mcmcS(X=IL[,1:156],Y=aggDRsum[1:156],n.iter=2000)
base100 <- baseline(X=IL[,estimation1],A=MCMCperiod1$COEF,nn=100)
leslie <- sequential(X=IL[,period1],Y=aggDRmean[period1])
nn <- nrow(DR1[,period1])
ind <- which(rownames(IL)%in%leslie$sel)
baseleslie <- apply(IL[ind,estimation1],2,mean)*nn
######################################
pdf("estimation.pdf")
#matplot(t(base100),type="l",col=1)
plot(Ytrue[estimation1],col=2,lwd=2,xlab="",ylab="",type="l")
lines(baseleslie,col=3,lwd=2)
#lines(base100[100,],col=4,lwd=2)
lines(apply(base100,2,mean),col=4,lwd=2)
abline(v=c(157,166))
legend("topleft",legend=c("vraie","est1","est2"),lty=1,lwd=2,col=c(2,3,4))
dev.off()

pdf("estimationIC.pdf")
matplot(t(base100),type="l",col=1,xlim=c(150,176),ylim=c(80,165),ylab="")
lines(Ytrue[estimation1],col=2,lwd=2)
lines(baseleslie,col=3,lwd=2)
#lines(base100[100,],col=4,lwd=2)
lines(apply(base100,2,mean),col=4,lwd=2)
abline(v=c(157,166))
legend("topleft",legend=c("vraie","est1","est2"),lty=1,lwd=2,col=c(2,3,4))
dev.off()

pdf("IC.pdf")
par(mfrow=c(1,2))
ICbas <- apply(base100,2,quantile,probs=0.025)
IChaut <- apply(base100,2,quantile,probs=0.975)
plot(Ytrue[estimation1],col=2,lwd=2,type="l",xlim=c(150,176),ylim=c(80,165),xlab="",ylab="")
matlines(157:166,t(base100[,157:166]),type="l",col=1)
abline(v=c(157,166))
lines(157:166,ICbas[157:166],col=4,lwd=3)
lines(157:166,IChaut[157:166],col=4,lwd=3)
lines(Ytrue[estimation1],col=2,lwd=2)
#legend("topleft",legend=c("vraie","IC","IC"),lty=1,lwd=2,col=c(2,4,4))

base500 <- baseline(X=IL[,estimation1],A=MCMCperiod1$COEF,nn=500)
ICbas <- apply(base500,2,quantile,probs=0.025)
IChaut <- apply(base500,2,quantile,probs=0.975)
plot(Ytrue[estimation1],col=2,lwd=2,type="l",xlim=c(150,176),ylim=c(80,165),
     xlab="",ylab="")
matlines(157:166,t(base500[,157:166]),type="l",col=1)
abline(v=c(157,166))
lines(157:166,ICbas[157:166],col=4,lwd=3)
lines(157:166,IChaut[157:166],col=4,lwd=3)
lines(Ytrue[estimation1],col=2,lwd=2)
#legend("topleft",legend=c("vraie","IC","IC"),lty=1,lwd=2,col=c(2,4,4))
dev.off()



matplot(t(base100),type="l",col=1,xlim=c(150,176),ylim=c(80,165))
lines(Ytrue[estimation1],col=2,lwd=2)
lines(baseleslie,col=3,lwd=2)
#lines(base100[100,],col=4,lwd=2)
lines(apply(base100,2,mean),col=4,lwd=2)
abline(v=c(157,166))
legend("topleft",legend=c("vraie","est1","est2"),lty=1,lwd=2,col=c(2,3,4))
dev.off()



plot(aggDRsum[estimation1],col=3,lwd=2)
###################
Ytrue <- apply(DR1,2,sum)
plot(Ytrue[estimation1],col=2,type="l",lwd=2)
lines(aggDRsum[dr1],col=2)
lines(aggDRsum[1:196],col=2)
lines(baseleslie,col=3,lwd=2)
#lines(base100[100,],col=4,lwd=2)
lines(apply(base100,2,mean),col=5,lwd=2)
abline(v=c(157,166))

###graphs
par(mfrow=c(2,2))
plot(aggDRsum[estimation1],type="l",col=1,lwd=2)
lines(apply(DR1,2,sum)[estimation1],col=2)
###################################################
MCMCperiod1 <- mcmcS(X=IL[,1:156],Y=aggDRsum[1:156],n.iter=2000)
base <- baseline(X=IL[,estimation1],A=MCMCperiod1$COEF)
intC <- IC(IL[,estimation1],A=MCMCperiod1$COEF)
intC2 <- IC2(IL[,estimation1],A=MCMCperiod1$COEF)
leslie <- sequential(X=IL[,period1],Y=aggDRmean[period1])
nn <- nrow(DR1[,period1])
ind <- which(rownames(IL)%in%leslie$sel)
baseleslie <- apply(IL[ind,estimation1],2,mean)*nn
####################################################
matplot(t(intC2),type="l",col=1)
lines(apply(intC,2,mean),col=2,lwd=2)
abline(v=c(157,166))
lines(aggDRsum[estimation1],col=3,lwd=2)
####################################################
plot(apply(DR1,2,sum)[estimation1],col=4,type="l")
lines(aggDRsum[estimation1],col=3)
lines(base,lwd=2,col=1)
lines(baseleslie,col=2)
#####################################################
Ytrue <- apply(DR1,2,sum)
plot(Ytrue[dr1],col=4,type="l")
lines(aggDRsum[dr1],col=3)
lines(base[dr1],lwd=2,col=1)
lines(baseleslie[dr1],col=2)

mean((base[dr1]-Ytrue[dr1])^2)
mean((baseleslie[dr1]-Ytrue[dr1])^2)





###DR period2
period2 <- mcmc(X=IL[,167:515],Y=aggDRsum[167:515],n.iter=1000)
base <- baseline(X=IL[,1:580],A=period2$COEF)
leslie <- sequential(X=IL[,167:515],Ym=aggDRmean[167:515])
nn <- nrow(DR1)
ind <- which(rownames(IL)%in%leslie$sel)
baseleslie <- apply(IL[ind,1:580],2,mean)*nn

plot(base,type="l")#,ylim=c(0,max(base)))
lines(baseleslie,col=2)
lines(aggDRsum[167:580],col=3)
lines(apply(DR1,2,sum)[167:580],col=4)

plot(base[500:540],type="l")
lines(baseleslie[500:540],col=2)
Ytrue <- apply(DR1,2,sum)
lines(Ytrue[500:540],col=4,lwd=3)

period <- 516:525
mean((base[period]-Ytrue[period])^2)
mean((baseleslie[period]-Ytrue[period])^2)



##A POSTERIORI
period2 <- mcmc(X=IL[,c(167:515,525:580)],Y=aggDRsum[c(167:515,525:580)],n.iter=500)
base <- baseline(X=IL[,167:580],A=period2$COEF)
plot(base,type="l")#,ylim=c(0,max(base)))
lines(period2$yfit,col=1,lty=2)
lines(aggDRsum[167:580],col=2)
lines(apply(DR1,2,sum)[167:580],col=4)



mm <- apply(CI,2,mean)
vv <- apply(CI,2,var) + mean(1/gg[401:500])
lb <- mm-2*sqrt(vv)
ub <- mm+2*sqrt(vv)

plot(c(1,100), c(min(lb),max(ub)), type="n",
     xlab="time", ylab="predicted consumption")
for ( k in 1:100 ){
  lines(c(k,k),c(lb[k],ub[k]),lwd=2,col="gray")
}

points(1:100,Y2$V1)
abline(v=seq(10.5,99,by=10), lty=3 )
