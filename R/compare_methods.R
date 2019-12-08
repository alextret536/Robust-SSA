library(Rssa)
library(pcaL1)
library(pcaMethods)
library(matrixcalc)

M <- 1

MSE.SSA <- rep(0,M)
MSE.L1svd <- rep(0,M)
MSE <- rep(0,M)
MSE1 <- rep(0,M)
MSE2 <- rep(0,M)

MAD.SSA <- rep(0,M)
MAD.L1svd <- rep(0,M)
MAD <- rep(0,M)
MAD1 <- rep(0,M)
MAD2 <- rep(0,M)


for(k in 1:M)
{
  
  
  sig <- (1:N)*exp(4*(1:N)/N)*sin(2*pi*(1:N)/30)
  
  sig.outl<-sig
  outlier.seq<-sample(1:(N),N*0.05)
  sig.outl[outlier.seq]<-sig.outl[outlier.seq]+1.5*sig.outl[outlier.seq]
  ser<-sig.outl+30*exp(4*(1:N)/N)*rnorm(N)
  
  #plot(ser,type='l')
  X<-hankel(ser, L=120)
  
  ############################### #SSA120.3 ##############################
  
  s <- ssa(ser, L = 120)
  rec <- reconstruct(s, groups = list(c(1:2)))
  
  # plot(s, type = "vectors",idx = 1:10)
  # plot(rec, add.residuals = TRUE, add.original = TRUE,
  #  plot.method = "xyplot", superpose = TRUE, auto.key = list(columns = 2))
  #plot(wcor(s))
  trend.season <- rec$F1
  
  ############################# #pcaL1(120.3)###############################
  X<-hankel(ser,L=120)
  s.L1svd<-l1pca(X,center=FALSE,projections="l1",projDim=2)
  
  Pr<-s.L1svd$projPoints
  Pr.L1svd<-hankL1(Pr)
  
  ###################### WLS ########################
  Pr<-QR_WLS(X,2,'loess')
  Pr0<-hankL1(Pr)
  
  Pr<-QR_WLS(X,2,'median')
  Pr1<-hankL1(Pr)
  
  Pr<-QR_WLS(X,2,'real')
  Pr2<-hankL1(Pr)
  
  
  #MSE
  MSE.SSA[k] <- mean((sig - trend.season)[1:N]^2)
  MSE.L1svd[k] <- mean((sig - Pr.L1svd)[1:N]^2)
  MSE[k] <- mean((sig - Pr0)[1:N]^2)
  MSE1[k] <- mean((sig - Pr1)[1:N]^2)
  MSE2[k] <- mean((sig - Pr2)[1:N]^2)
  
  #MAD
  MAD.SSA[k] <- mean(abs((sig - trend.season)[1:N]))
  MAD.L1svd[k] <- mean(abs((sig - Pr.L1svd)[1:N]))
  MAD[k] <- mean(abs((sig - Pr0)[1:N]))
  MAD1[k] <- mean(abs((sig - Pr1)[1:N]))
  MAD2[k] <- mean(abs((sig - Pr2)[1:N]))
  
}

RMSE.SSA<-sqrt(mean(MSE.SSA))
RMSE.L1svd<-sqrt(mean(MSE.L1svd))
RMSE<-sqrt(mean(MSE))
RMSE1<-sqrt(mean(MSE1))
RMSE2<-sqrt(mean(MSE2))

MAD.SSA<-mean(MAD.SSA)
MAD.L1svd<-mean(MAD.L1svd)
MAD<-mean(MAD)
MAD1<-mean(MAD1)
MAD2<-mean(MAD2)

RMSE.SSA #ошибка SSA
RMSE.L1svd #ошибка последовательного метода L1svd
RMSE #ошибка IRLS с выделением тренда из модуля остатков с помощью loess
RMSE1 #ошибка IRLS с выделением тренда из модуля остатков с помощью скользящей медианы
RMSE2 #ошибка IRLS со знанием реального тренда из модуля остатков



plot(sig,type='l')
lines(ser,col='gray')
lines(trend.season,type='l',col='blue',lw=2)
lines(Pr.L1svd,type='l',col='red',lw=2)
lines(Pr0,type='l',col='violet',lw=2)
lines( Pr1,type='l',col='green',lw=2)
lines( Pr2,type='l',col='orange',lw=2)
legend('topleft', c("SSA","L1svd","IRLS.loess", "IRLS.median", "IRLS.real trend"),
       col=c("blue","red","violet","green","orange"), lty=1, cex=0.8, lw=2)