#---------------------------------------------------------
# Liermann's version of Integrated Watershed Area Model
# Steps
# 1. Read in stock-recruitment data
# 2. Create data and parameter lists for TMB
# 3. Estimate SR parameters and SMSY & SREP for synoptic data sets in TMB
# 4. Caldulate diagnostics for SR models and plot SR curves, etc.
# 5. Read in watershed areas

#---------------------------------------------------------
# Libaries

library(rsample)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(TMB)
library(zoo)

# Functions
count.dig <- function(x) {floor(log10(x)) + 1}
'%not in%' <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

plot <- FALSE
removeSkagit <- FALSE

source ("PlotSR.r")# Plotting functions


#---------------------------------------------------------
# 1. Read in data

SRDatwNA <- read.csv("DataIn/SRinputfile.csv")
SRDatwNA <- SRDatwNA %>% filter(Name != "Hoko" & Name != "Hoh") #remove two stocks not used in Parken et al, and not documented in Liermann et al.
if (removeSkagit==TRUE) {
  SRDatwNA <- SRDatwNA %>% filter(Name != "Skagit")#Stocknumber=22. Need to re-align stock numbers of last two stocks, 23 and 24
  SRDatwNA [which(SRDatwNA$Stocknumber==23),2] = 22
  SRDatwNA [which(SRDatwNA$Stocknumber==24),2] = 23
}

# Which stocks have NAs?
stockwNA <- SRDatwNA %>% filter (is.na(Rec) == TRUE) %>% select (Stocknumber) %>% unique() %>% unlist() 
#Do not use AR(1) model on  stocks with NAs, Humptulips and Queets (20 and 21)

# Remove years with NAs
SRDat <- SRDatwNA %>% filter(Rec != "NA") #%>% filter(Stocknumber <= 24)

# Revise yr_num list where NAs have been removed
test <- SRDat %>% filter(Stocknumber == stockwNA[1]| Stocknumber == stockwNA[2])
if( max(SRDat$Stocknumber) >= stockwNA[1]) {
  for (i in 1:length(stockwNA)) {
    len <- length (SRDat [which (SRDat$Stocknumber == stockwNA[i]), ]$yr_num) - 1
    SRDat [which (SRDat$Stocknumber == stockwNA[i]), ]$yr_num <- c (0:len)
  }
}
test <- SRDat %>% filter(Stocknumber == stockwNA[1]| Stocknumber == stockwNA[2])


# Calculate scale for each stock
digits <- SRDat %>% group_by(Stocknumber) %>% summarize(maxDigits = count.dig(max(Sp)))
SRDat <- left_join(SRDat, digits)
SRDat <- SRDat %>% mutate(Scale = 10^(maxDigits-1))



#remove years 1981-1984, 1986-1987  from Cowichan (Stocknumber 23) as per Tompkins et al. 2005
SRDat_Cow <- SRDat %>% filter(Name == "Cowichan" & Yr >= 1985 & Yr !=1986 & Yr != 1987) 
n_Cow <- length(SRDat_Cow$Yr)
SRDat_Cow$yr_num <- 0:(n_Cow-1)
#if(stksNum_surv == 23) SRDat_surv <- SRDat_surv_Cow
SRDat <- SRDat %>% filter(Name != "Cowichan") %>% bind_rows(SRDat_Cow) %>% arrange(Stocknumber)

# Read in watershed area data and life-history type (stream vs ocean)
WA <- read.csv("DataIn/WatershedArea.csv")
names <- SRDat %>% select (Stocknumber, Name) %>% distinct() 
WA <- WA %>% full_join(names, by="Name") %>% arrange(Stocknumber)
if (removeSkagit==TRUE) {WA <- WA %>% filter(Name !="Skagit")}
Stream <- SRDat %>% select(Stocknumber, Name, Stream)  %>% distinct()

# 2. Create data and parameter lists for TMB

TMB_Inputs <- list(rho_Start = 0.0, logDelta1_start=3.00, logDelta2_start =log(0.72), logDeltaSigma_start = -0.412, 
                   logMuDelta1_mean= 5, logMuDelta1_sig= 10, logMuDelta2_mean=-0.5, logMuDelta2_sig= 10, 
                   Tau_Delta1_dist= 0.1, Tau_Delta2_dist= 0.1) 


# Data 
data <- list()
Scale <- SRDat$Scale 
data$S <- SRDat$Sp/Scale 
data$logRS <- log( (SRDat$Rec/Scale) / (SRDat$Sp/Scale) )
data$stk <- SRDat$Stocknumber
data$yr <- SRDat$yr_num
N_Stocks <- length(unique(SRDat$Name))

# Read in wateshed area data and life-history type....
#data$WA <- WA$WA
data$Scale <- SRDat %>% select(Stocknumber, Scale) %>% distinct() %>% select(Scale) %>% pull(Scale)
#data$Stream <- Stream$Stream
#data$N_stream <-length(which(data$Stream==0))
#data$N_ocean <- length(which(data$Stream==1))

# What does inv gamma prior look like? library(invgamma); plot(x=seq(0,1,0.001), y=dinvgamma(seq(0,1,0.001),0.01,0.01), type="l")

# Data needed for hierarchical Lietrmann WA model
#data$logMuDelta1_mean <- TMB_Inputs$logMuDelta1_mean
#data$logMuDelta1_sig <- TMB_Inputs$logMuDelta1_sig
#data$logMuDelta2_mean <- TMB_Inputs$logMuDelta2_mean
#data$logMuDelta2_sig <- TMB_Inputs$logMuDelta2_sig
#data$Tau_Delta1_dist <- TMB_Inputs$Tau_Delta1_dist
#data$Tau_Delta2_dist <- TMB_Inputs$Tau_Delta2_dist

#Data need to hierarchial Ricker a
data$Tau_dist <- 0.1
data$logMuA_mean <- 1.5
data$logMuA_sig <- 1.5
data$Tau_A_dist <- 0.1


  
# Parameters
param <- list()
param$logA <- SRDat %>% group_by (Stocknumber) %>% summarise(yi = lm(log( Rec / Sp) ~ Sp )$coef[1] ) %>% pull(yi) 
B <- SRDat %>% group_by(Stocknumber) %>% summarise( m = - lm(log( Rec / Sp) ~ Sp )$coef[2] )
param$logB <- log ( 1/ ( (1/B$m)/data$Scale ))
param$logSigma <- rep(-2, N_Stocks)
param$logMuA <- 1.5
param$logSigmaA <- 1

# Parameters for WA model

#param$logDelta1 <- 3.00# with skagit 2.881
##param$logDelta2 <- log(0.72)#log(0.72/(1-0.72)) #logit 0f 0.72 #with skagit logDelta2 = -0.288
#param$Delta2 <- log(0.72/(1-0.72)) #logit 0f 0.72 #with skagit logDelta2 = -0.288
#param$logDeltaSigma <- -0.412 #from Parken et al. 2006 where sig=0.662
# without Skagit lnDelta1_start <- 2.999911
# without Skagit lnDelta2_start <- -0.3238648, or Delta2 = 0.723348

## Separate Stream and Ocean type models
#param$slogDelta1 <- 2.744 #best estimates from run of stream-specific WAregression TMB model run
#param$sDelta2 <- 0.857 
#param$slogDeltaSigma <- -0.709 

#param$ologDelta1 <- 3.00#1.519 #best estimates from run of stream-specific WAregression TMB model run
#param$ologDelta2 <- log(0.94)#0#21.2 
#param$ologDeltaSigma <-  -0.412#-0.94 

## Liermann model
#param$logDelta1 <- 10# with skagit 2.881
#param$logDelta1ocean <- 0# with skagit 2.881
#param$logDelta2 <- log(0.72)#log(0.72/(1-0.72)) #logit 0f 0.72 #with skagit logDelta2 = -0.288
#param$logDelta2ocean <- 0#log(0.72/(1-0.72)) #logit 0f 0.72 #with skagit logDelta2 = -0.288
#param$logDeltaSigma <- -0.412 #from Parken et al. 2006 where sig=0.662

## Hierarchcial WA model

#param$logDelta1 <- TMB_Inputs$logDelta1_start#rep(TMB_Inputs$logDelta1_start, 2)
#param$logDelta2 <- rep(TMB_Inputs$logDelta2_start, 2)
#param$logDeltaSigma <-TMB_Inputs$logDeltaSigma_start 
#param$logMuDelta1 <- TMB_Inputs$logDelta1_start
#param$SigmaDelta1 <- 10
#param$logMuDelta2 <- TMB_Inputs$logDelta2_start
#param$SigmaDelta2 <- 1



# 3. Estimate SR parameters from synoptic data set and SMSY and SREPs

# Compile model if changed:
#dyn.unload(dynlib("TMB_Files/Liermann"))
#compile("TMB_Files/Liermann.cpp")

dyn.load(dynlib("TMB_Files/Liermann"))

# For Phase 1, fix Delta parameters. Do not need to fix Delta's beccause initlal values are lm fits, so very close
#map <- list(logDelta1=factor(NA), Delta2=factor(NA), logDeltaSigma=factor(NA)) 
#obj <- MakeADFun(data, param, DLL="Ricker_AllMod", silent=TRUE, map=map)

obj <- MakeADFun(data, param, DLL="Liermann", random = "logA", silent=TRUE)
#upper <- c(rep(Inf, 80), 5.00, rep(Inf,2))
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))#, upper=upper)
pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1
summary(sdreport(obj), p.value=TRUE)

#library(tmbstan)
#fitmcmc <- tmbstan(obj, chains=3, iter=1000, init=list(opt$par), control = list(adapt_delta = 0.95))

# Create Table of outputs
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)


# Put together readable data frame of values
All_Ests$Param <- sapply(All_Ests$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))

# By first spliting out stocks modelled with standard Ricker, RickerAR(1), and Ricker-survival models
All_Est <- All_Ests %>% filter (Param %in% c("logA", "logB", "logSigma",  "SMSY", "SREP"))
SN <- unique(SRDat[, c("Stocknumber")])
All_Est$Stocknumber <- rep(SN)
All_Est <- left_join(All_Est, unique(SRDat[, c("Stocknumber", "Name")]))


# 4. Calculate diagnostics and plot SR curves, etc.

# Calculate AIC

nLL <- data.frame(nLL=obj$report()$nLL) %>% add_column(Stocknumber=SRDat$Stocknumber) %>% group_by(Stocknumber) %>% summarize(CnLL=sum(nLL))
aic <- nLL %>% mutate(aic = 2 * 3 + 2*CnLL) # 

# Get predicted values and calculate r2
Pred <- data.frame()
Pred <- All_Ests %>% filter (Param %in% c("LogRS_Pred"))
Preds <- SRDat %>% select("Stocknumber","yr_num", "Sp", "Rec", "Scale", "Name") %>% add_column(Pred=Pred$Estimate)
Preds <- Preds %>% mutate(ObsLogRS = log ( (Rec / Scale) / (Sp/Scale) ) )
r2 <- Preds %>% group_by(Stocknumber) %>% summarize(r2=cor(ObsLogRS,Pred)^2)


# Calculate standardized residuals
SRes <- Preds %>% arrange (Stocknumber)
SRes <- SRes %>% mutate ( Res = ObsLogRS- Pred) #%>% mutate (StdRes = Res/??)
sigma <- All_Est %>% filter(Param=="logSigma") %>% select(Stocknumber, Estimate, Name)
SRes <- SRes %>% left_join(sigma) %>% rename(logSig = Estimate)
SRes <- SRes %>% mutate (StdRes = Res/exp(logSig))

SMSY <- All_Est %>% filter(Param=="SMSY") #%>% mutate(ModelOrder=0:(length(unique(All_Est$Stocknumber))-1))
# what is scale of SMSY?
Sc <- SRDat %>% select(Stocknumber, Scale) %>% distinct()
SMSY <- SMSY %>% left_join(Sc) %>% mutate(rawSMSY=Estimate*Scale)

#Plot SR curves. linearized model, standardized residuals, autocorrleation plots for synoptic data set
if (plot==TRUE){ #Need to fix these for straight standard model
  PlotSRCurve(SRDat=SRDat, All_Est=All_Est, SMSY_std=SMSY, stksNum_ar=NA, stksNum_surv=NA, r2=r2, removeSkagit=removeSkagit)
  PlotSRLinear(SRDat=SRDat, All_Est=All_Est, SMSY_std=SMSY, stksNum_ar=NA, stksNum_surv=NA, r2=r2, removeSkagit=removeSkagit) 
  PlotStdResid(SRes)
  Plotacf(SRes)
  
}


# What initial values to use for WA model parameters?

lnSMSY <- log(SMSY$rawSMSY)
lnWA <- log(WA$WA)
ParkenSMSY <- as.data.frame(read.csv("DataIn/ParkenSMSY.csv")) %>% mutate(lnSMSY=log(SMSY))
lnPSMSY <- ParkenSMSY$lnSMSY

# lm(lnPSMSY ~ lnWA) #Get same coefficients as Parken et al. Table 4 for pooled data
lnDelta1_start <- coef(lm(lnSMSY ~ lnWA))[1]
lnDelta2_start <- log(coef(lm(lnSMSY ~ lnWA))[2])
# without Skagit lnDelta1_start <- 2.999911
# without Skagit lnDelta2_start <- -0.3238648, or Delta2 = 0.723348
# With Skagit lnDelta1_start <- 2.881
# with Skagit nDelta2_start <- -0.288

plot(y=lnSMSY, x=lnWA)
plot(y=exp(lnSMSY), x=exp(lnWA))
#pdf("ParkenSMSYWA.pdf", width=4)
#  par(mfcol=c(2,1))
#  plot(y=exp(lnPSMSY), x=exp(lnWA), xlab="Watershed Area, km2", ylab="SMSY, Parken et al. 2006")
#  plot(y=exp(lnPSMSY), x=exp(lnWA), xlim=c(0,2000), ylim=c(0,6000), xlab="Watershed Area, km2", ylab="SMSY, Parken et al. 2006")
#dev.off()

#-----------------------------------------------------------------------------------------
test <- FALSE
if(test==TRUE){

#TEST: Watereshed-Area Regression with data inputs
SMSY <- read.csv("DataIn/SMSY_3mods.csv")
Sca <- SMSY %>% select(Name, Scale)
SMSY <- SMSY %>% left_join(Stream) %>% left_join(WA)
SMSY_stream <- SMSY %>% filter(lh==0)
SMSY_ocean <- SMSY %>% filter(lh==1)

ParkenSMSY_stream <- read.csv("DataIn/ParkenSMSY.csv") %>% right_join(Stream) %>% filter(lh==0) %>% mutate(Estimate=SMSY/SMSY_stream$Scale)

data <- list()
data$SMSY <- SMSY_ocean$Estimate#SMSY$Estimate
data$WA <- SMSY_ocean$WA#SMSY$WA#exp(lnWA)
data$Scale <- SMSY_ocean$Scale#SMSY$Scale
#data$Tau_dist <- 0.1

param <- list()
param$logDelta1 <- 3.00# with skagit 2.881
param$Delta2 <- log(0.72/(1-0.72)) #logit 0f 0.72 #with skagit logDelta2 = -0.288
param$logDeltaSigma <- -0.412 #from Parken et al. 2006 where sig=0.662

#dyn.unload(dynlib("TMB_Files/WAregression"))
#compile("TMB_Files/WAregression.cpp")
plot(x=log(data$WA), y=log(data$SMSY*data$Scale), xlab="ln(WA)", ylab="ln(SMSY)")
abline(lm(log(data$SMSY*data$Scale) ~ log(data$WA)))
plot(x=(data$WA), y=(data$SMSY*data$Scale), xlab="WA, km2", ylab="SMSY")

dyn.load(dynlib("TMB_Files/WAregression"))
obj <- MakeADFun(data, param, DLL="WAregression", silent=TRUE)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))#, upper=upper)
summary(sdreport(obj), p.value=TRUE)
# I get same answers as Parken when I use his SMSY data, for streams (haven't checked ocean, as Skagit is ocean type and numbering is wonky)
# For my SMSY_stream data slogDelta1 <- 2.744
# For my SMSY_stream data sDelta2 <- 0.857
# For my SMSY_stream data slogDeltaSigma <- -0.709

# make sure my output from the TMB model SMSY_stream and WA_stream match the values from:
#All_Est %>% filter(Param=="SMSY") %>% left_join(WA) %>% left_join(Stream) %>% filter(lh==1)
# Yes, they match

# WHen I run ocean and stream-type specific regressions, the SMSY values come out as a line! 
test <- All_Est %>% left_join(Sca) %>% filter(Param=="SMSY") %>% left_join(WA) %>% left_join(Stream) %>% filter(lh==1)#switch to lh==0 for stream
lm(log(test$Estimate*test$Scale) ~ log(test$WA))
plot( y= log(test$Estimate*test$Scale), x=log(test$WA))

}
#-----------------------------------------------------------------------------------------



