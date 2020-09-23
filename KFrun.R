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


stk <- unique(SRDat$Stocknumber)
ind <- 0
nstk <- length(stk)
SDKFSmsy <- NA
use <- 0

for (i in stk){
  ind <- ind + 1
  Rec <- SRDat %>% filter (Stocknumber==i) %>% pull (Rec) #/10000
  Sp <- SRDat %>% filter (Stocknumber==i) %>% pull (Sp) #/10000

  KFa <- kf.rw(initial=initial, x=Sp, y=log(Rec/Sp))$smoothe.mean.a
  KFb <- kf.rw(initial=initial, x=Sp, y=log(Rec/Sp))$b
  KFsmsy <- (1 - gsl::lambert_W0(exp(1 - (KFa)))) / -(rep(KFb,length(KFa)))
  SDKFSmsy[ind] <- sd(KFsmsy)
  
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
  
}



KFsmsy

#data <- read.csv("data/Stellacko.csv") #From KF funcs package
#x <- data$ETS
#y <- log(data$Rec/data$ETS)


KFa <- kf.rw(initial=initial, x=data$ETS, y=log(data$Rec/data$ETS))$smoothe.mean.a
KFb <- kf.rw(initial=initial, x=data$ETS, y=log(data$Rec/data$ETS))$b

kfAICc <- kf.rw(initial=initial, x=data$ETS, y=log(data$Rec/data$ETS))$AICc # uses concentrated nLL
kfAICc


linearnLL <- -logLik( lm ( log (data$Rec/data$ETS) ~ data$ETS ))[1] #full nLL
rss <- sum(lm ( log (data$Rec/data$ETS) ~ data$ETS )$resid^2)
sigEst <- sd(lm ( log (data$Rec/data$ETS) ~ data$ETS )$resid)
linearnLLconc <- - ( -(N/2)*log(sigEst^2) - (1/(2*sigEst^2)) * rss ) #Concentrated nLL
lin_param <- 3
N <- length(data$ETS) #assuming complete data
linAICc <- 2 * linearnLLconc + 2 * lin_param * (N/(N - lin_param - 1))
linAICc

if( linAICc > (kfAIC-2) ) use <- 1

