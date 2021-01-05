# Code to estimate time-varying Ricker log(a), beta, and time-varying SMSY, and associated SD

#---------------------------------------------------------
# Libaries

library(rsample)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(TMB)
library(zoo)
library(viridis)
library(hrbrthemes)
library(gsl)

#---------------------------------------------------------
# Source code
source("R/KFfunctions.R")

# Initial values for all stocks
initial<-list()
initial$mean.a<-1#2
initial$var.a<-1#2
initial$b <- 2
initial$ln.sig.e<-1
initial$ln.sig.w<-1
initial$Ts<-0
initial$EstB<-TRUE


SRDat <- read.csv("DataIn/WCVIrunreconstruction.csv")
# Steps to do:
# 1. plot SR data with Ricker curve
# 2. PLot time-series of R/S
# 3. Plot KF time-series (straight line)

#png(paste("DataOut/KFSMSY.png", sep=""), width=7, height=7, units="in", res=500)


SDKFlSmsy <- NA
SDKFlSrep <- NA
#plot time-varying Ric A



Rec <- SRDat$Rec
Sp <- SRDat$Sp
Yr <- SRDat$Yr
  
KFa <- kf.rw(initial=initial, x=Sp, y=log(Rec/Sp))$smoothe.mean.a
KFb <- kf.rw(initial=initial, x=Sp, y=log(Rec/Sp))$b
KFsmsy <- (1 - gsl::lambert_W0(exp(1 - (KFa)))) / -(rep(KFb,length(KFa)))
KFsrep <- KFa/-rep(KFb,length(KFa))
SDKFlSmsy <- sd(log(KFsmsy))
SDKFlSrep <- sd(log(KFsrep))
  
kfAICc <- kf.rw(initial=initial, x=Sp, y=log(Rec/Sp))$AICc # uses concentrated nLL
  
  
linearnLL <- -logLik( lm ( log (Rec/Sp) ~ Sp ))[1] #full nLL
rss <- sum(lm ( log (Rec/Sp) ~ Sp )$resid^2)
sigEst <- sd(lm ( log (Rec/Sp) ~ Sp )$resid)
lin_param <- 3
N <- length(Sp) #assuming complete data
linearnLLconc <- - ( -(N/2)*log(sigEst^2) - (1/(2*sigEst^2)) * rss ) #Concentrated nLL
linAICc <- 2 * linearnLLconc + 2 * lin_param * (N/(N - lin_param - 1))
linAICc
  
Use <- NA
if( linAICc > (kfAICc-2) ) use <- 1
  
par(cex=1.5)
plot(x=1987:2011, y=rep(0,length(1987:2011)), type="l", col="white", ylim=c(0,3), xlab="Brood year", ylab="Time varying log(SMSY)", xlim=c(1965,2000))
lines(x=Yr, y=KFa, lwd=2)

  


#dev.off()

