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
source("KFfunctions.R")



# Initial values for all stocks
initial<-list()
initial$mean.a<-1
initial$var.a<-1
initial$b <- 2
initial$ln.sig.e<-1
initial$ln.sig.w<-1
initial$Ts<-0
initial$EstB<-TRUE


SRDatwNA <- read.csv("DataIn/SRinputfile.csv")
SRDatwNA <- SRDatwNA %>% filter(Name != "Hoko" & Name != "Hoh") #remove two stocks not used in Parken et al, and not documented in Liermann et al.
# Which stocks have NAs?
stockwNA <- SRDatwNA %>% filter (is.na(Rec) == TRUE) %>% select (Stocknumber) %>% unique() %>% unlist() 
#Do not use Kf model on  stocks with NAs, Humptulips and Queets (20 and 21)
SRDat <- SRDatwNA %>% filter(Name != "Humptulips" & Name != "Queets")

png(paste("DataOut/KFSMSY.png", sep=""), width=7, height=7, units="in", res=500)


stk <- unique(SRDat$Stocknumber)
ind <- 0
nstk <- length(stk)
SDKFlSmsy <- NA
SDKFlSrep <- NA
use <- rep(0,nstk) #need to check that this is working
#plot time-varying Ric A
KFsmsy.ls <- list()

par(cex=1.5)
plot(x=1939:2000, y=rep(0,length(1939:2000)), type="l", col="white", ylim=c(4,12), xlab="Brood year", ylab="Time varying log(SMSY)", xlim=c(1965,2000))
#plot(x=1939:2000, y=rep(0,length(1939:2000)), type="l", col="white", ylim=c(4,12), xlab="Brood year", ylab="Time varying log(SMSY)")
#plot(x=1939:2000, y=rep(0,length(1939:2000)), type="l", col="white", ylim=c(4,12), xlab="Brood year", ylab="Time varying log(SREP)")

for (i in stk){
  ind <- ind + 1
  Rec <- SRDat %>% filter (Stocknumber==i) %>% pull (Rec)
  Sp <- SRDat %>% filter (Stocknumber==i) %>% pull (Sp)
  Yr <- SRDat %>% filter (Stocknumber==i) %>% pull (Yr)
  
  KFa <- kf.rw(initial=initial, x=Sp, y=log(Rec/Sp))$smoothe.mean.a
  KFb <- kf.rw(initial=initial, x=Sp, y=log(Rec/Sp))$b
  KFsmsy <- (1 - gsl::lambert_W0(exp(1 - (KFa)))) / -(rep(KFb,length(KFa)))
  KFsrep <- KFa/-rep(KFb,length(KFa))
  KFsmsy.ls[[ind]] <- KFsmsy
  KFsmsy[which(KFsmsy<0)] <- NA
  KFsrep[which(KFsrep<0)] <- NA
  SDKFlSmsy[ind] <- sd(log(KFsmsy))
  SDKFlSrep[ind] <- sd(log(KFsrep))
  
  kfAICc <- kf.rw(initial=initial, x=Sp, y=log(Rec/Sp))$AICc # uses concentrated nLL
  
  
  linearnLL <- -logLik( lm ( log (Rec/Sp) ~ Sp ))[1] #full nLL
  rss <- sum(lm ( log (Rec/Sp) ~ Sp )$resid^2)
  sigEst <- sd(lm ( log (Rec/Sp) ~ Sp )$resid)
  lin_param <- 3
  N <- length(Sp) #assuming complete data
  linearnLLconc <- - ( -(N/2)*log(sigEst^2) - (1/(2*sigEst^2)) * rss ) #Concentrated nLL
  linAICc <- 2 * linearnLLconc + 2 * lin_param * (N/(N - lin_param - 1))
  linAICc
  
  if( linAICc > (kfAICc-2) ) use[ind] <- 1
  
  if(use[ind]>0){
    lines(x=Yr, y=log(KFsmsy), lwd=2)
  }
  #lines(x=Yr, y=log(KFsrep))
  
}

dev.off()

summary(SDKFlSmsy, na.rm=T) # 4 straight lines
summary(SDKFlSrep, na.rm=T) # 4 straight lines

# Include only those stocks where log(KFSmsy) values vary over time
summary(SDKFlSmsy[which(SDKFlSmsy >0.01)])
hist(SDKFlSmsy[which(SDKFlSmsy >0.01)], breaks=50, xlab="SD of time-varying ln(SMSY)", main ="")
medSDlogSmsy <- median(SDKFlSmsy[which(SDKFlSmsy >0.01)]) # use this a min bound on sigmaDelta

# Include only those stocks where log(KFSrep) values vary over time
summary(SDKFlSrep[which(SDKFlSrep >0.01)])
hist(SDKFlSrep[which(SDKFlSrep >0.01)], breaks=50, xlab="SD of time-varying ln(SMSY)", main ="")
medSDlogSrep <- median(SDKFlSrep[which(SDKFlSrep >0.01)]) # use this a min bound on sigmaDelta

# Include on those stocks where AICc of KF model is < or within 2 units of AICc of linear model
summary(SDKFlSmsy[which(use>0)])
png(paste("DataOut/SDSMSYhist.png", sep=""), width=7, height=7, units="in", res=500)
par(cex.axis=1.5, cex.lab=1.5)
hist(SDKFlSmsy[which(use>0)], breaks=5, main="", xlab="SD of time-varying SMSY", col="light blue")
dev.off()
medSDlogSmsy <- median(SDKFlSmsy[which(use >0)], na.rm=T) # OR better yet use this a min bound on sigmaDelta

# Include on those stocks where AICc of KF model is < or within 2 units of AICc of linear model
summary(SDKFlSrep[which(use>0)])
hist(SDKFlSrep[which(use>0)], breaks=20)
medSDlogSrep <- median(SDKFlSrep[which(use >0)], na.rm=T) # OR better yet use this a min bound on sigmaDelta

# For max bounds, use SD of log(SMSY) among all 25 stocks.
SDlSMSYParken <- sd(log(read.csv("DataIn/ParkenSMSY.csv")$SMSY))
# For max bounds, use SD of log(SREP) among all 25 stocks.
SDlSREPParken <- sd(log(read.csv("DataIn/ParkenSREP.csv")$SREP))

meanSigmaDelta <- (SDlSMSYParken - medSDlogSmsy)/2 + medSDlogSmsy
meanSigmaNu <- (SDlSREPParken - medSDlogSrep)/2 + medSDlogSrep

test<- seq(medSDlogSmsy,SDlSMSYParken,length=50)
plot(x=test, y=dnorm(test, meanSigmaDelta,1), type="l", xlab="sigmaDelta", ylab="Probability Density", ylim=c(0,5))
lines(x=test, y=dnorm(test, meanSigmaDelta,0.5))
lines(x=test, y=dnorm(test, meanSigmaDelta,0.28))# With this sigma, 95% of probablity density is within bounds mean +/- 0.55 (where (SDlSMSYParken - medSDlogSmsy)/2 = 0.55)
#0.28*1.96 = 0.55

test<- seq(medSDlogSrep,SDlSREPParken,length=50)
plot(x=test, y=dnorm(test, meanSigmaNu,1), type="l", xlab="sigmaNu", ylab="Probability Density", ylim=c(0,5))
lines(x=test, y=dnorm(test, meanSigmaNu,0.5))
lines(x=test, y=dnorm(test, meanSigmaNu,0.275))# With this sigma, 95% of probablity density is within bounds mean +/- 0.55 (where (SDlSREPParken - medSDlogSrep)/2 = 0.54)
#0.275*1.96 = 0.54
