#---------------------------------------------------------
# Integrated Watershed Area Model
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
removeSkagit <- TRUE

if( plot== TRUE) {
  source ("PlotSR.r")# Plotting functions
  if(removeSkagit==FALSE) source ("CheckAR1.r")# For SR plotting purposes below, need to estimate std Ricker SMSY for AR1 stocks, "SMSY_std"
}

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


stks_ar <- c("Chikamin", "Keta", "Blossom", "Situk", "Siletz", "Columbia Sp")#Cowichan, stk-23, not included here becuase modelled as per Tompkins with a surival covariate
stksNum_ar <- c(4,5,6,10,11,16)

# Cowichan modeled with Ricker with a surival co-variate. 
# Harrison was modeled with survival co-variate, but gives very poor fit with very high gamma and so excluded 
stksNum_surv <- c(0,23)
stks_surv <- c("Harrison", "Cowichan")
if (removeSkagit==TRUE) {stksNum_surv <- c(0,22)}

len_stk <- length(unique(SRDat$Stocknumber))
stksNum_std <- which(0:(len_stk-1) %not in%c(stksNum_ar, stksNum_surv)==TRUE)-1 # Assuming there are only 25 stocks (0:24 StockNumber)

# When aggregated standard, ar1, surv, this is the order of stocks
stksOrder <- data.frame(Stocknumber =  c(stksNum_std, stksNum_ar, stksNum_surv), ModelOrder = 0:(len_stk-1))

# What is the scale of S,R and SMSY,SREP data, ordered when aggregated by std, AR1, surv
SRDat_Scale <- SRDat %>% select(Stocknumber, Scale) %>% distinct() 
SRDat_Scale <- SRDat_Scale %>% left_join(stksOrder) %>% arrange(ModelOrder)
SRDat_Scale <- SRDat_Scale$Scale

SRDat_std <- SRDat %>% filter(Stocknumber %not in% c(stksNum_ar,stksNum_surv)) 
SRDat_ar <- SRDat %>% filter(Stocknumber %in% stksNum_ar) 
SRDat_surv <- SRDat %>% filter(Stocknumber %in% stksNum_surv) 


# Assign new stock numbers to each stock so that they are sequential. 
ind_std <- tibble(ind_std= 0:(length(unique(SRDat_std$Name))-1))
ind_std <- add_column(ind_std, Stocknumber = (unique(SRDat_std$Stocknumber)))
SRDat_std <- SRDat_std %>% left_join(ind_std)

ind_ar <- tibble(ind_ar= 0:(length(unique(SRDat_ar$Name))-1))
ind_ar <- add_column(ind_ar, Stocknumber = (unique(SRDat_ar$Stocknumber)))
SRDat_ar <- SRDat_ar %>% left_join(ind_ar)

ind_surv <- tibble(ind_surv= 0:(length(unique(SRDat_surv$Name))-1))
#ind_surv <- add_column(ind_surv, Stocknumber = (unique(SRDat_surv$Stocknumber)))
ind_surv <- add_column(ind_surv, Name = (unique(SRDat_surv$Name)))
Surv <- as.data.frame(read.csv("DataIn/Surv.csv")) #%>% filter(Yr<=1999)
SRDat_surv <- SRDat_surv %>% left_join(ind_surv, by="Name") %>% left_join(Surv)


#remove years 1981-1984, 1986-1987  from Cowichan (Stocknumber 23) as per Tompkins et al. 2005
SRDat_surv_Cow <- SRDat_surv %>% filter(Name == "Cowichan" & Yr >= 1985 & Yr !=1986 & Yr != 1987) 
n_surv_Cow <- length(SRDat_surv_Cow$Yr)
SRDat_surv_Cow$yr_num <- 0:(n_surv_Cow-1)
#if(stksNum_surv == 23) SRDat_surv <- SRDat_surv_Cow
if("Cowichan" %in% stks_surv) SRDat_surv <- SRDat_surv %>% filter(Name != "Cowichan") %>% bind_rows(SRDat_surv_Cow)

# Read in watershed area data and life-history type (stream vs ocean)
WA <- read.csv("DataIn/WatershedArea.csv")
names <- SRDat %>% select (Stocknumber, Name) %>% distinct() 
WA <- WA %>% full_join(names, by="Name") %>% full_join (stksOrder, by="Stocknumber") %>% arrange(ModelOrder)
if (removeSkagit==TRUE) {WA <- WA %>% filter(Name !="Skagit")}
Stream <- SRDat %>% select(Stocknumber, Name, Stream) %>% group_by(Stocknumber) %>% summarize(lh=max(Stream))
Stream <- Stream %>% full_join(stksOrder, by="Stocknumber") %>% arrange(ModelOrder)

# 2. Create data and parameter lists for TMB

TMB_Inputs <- list(rho_Start = 0.0, logDelta1_start=3.00, logDelta2_start =log(0.72), logDeltaSigma_start = -0.412, 
                   logMuDelta1_mean= 5, logMuDelta1_sig= 10, logMuDelta2_mean=-0.5, logMuDelta2_sig= 10, 
                   Tau_Delta1_dist= 0.1, Tau_Delta2_dist= 0.1) 


# Data 
data <- list()
Scale_std <- SRDat_std$Scale 
data$S_std <- SRDat_std$Sp/Scale_std 
data$logRS_std <- log( (SRDat_std$Rec/Scale_std) / (SRDat_std$Sp/Scale_std) )
data$stk_std <- as.numeric(SRDat_std$ind_std)
N_Stocks_std <- length(unique(SRDat_std$Name))
data$yr_std <- SRDat_std$yr_num

Scale_ar <- SRDat_ar$Scale 
data$S_ar <- SRDat_ar$Sp/Scale_ar 
data$logRS_ar <- log( (SRDat_ar$Rec/Scale_ar) / (SRDat_ar$Sp/Scale_ar) ) 
data$stk_ar <- as.numeric(SRDat_ar$ind_ar)
N_Stocks_ar <- length(unique(SRDat_ar$Name))
data$yr_ar <- SRDat_ar$yr_num

Scale_surv <- SRDat_surv$Scale 
data$S_surv <- SRDat_surv$Sp/Scale_surv
data$logRS_surv <- log( (SRDat_surv$Rec/Scale_surv) / (SRDat_surv$Sp/Scale_surv) )
data$stk_surv <- as.numeric(SRDat_surv$ind_surv)
N_Stocks_surv <- length(unique(SRDat_surv$Name))
data$yr_surv <- SRDat_surv$yr_num
data$Surv_surv <- log(SRDat_surv$Surv) #Tompkins et al. used Ln(Surv+1)
meanLogSurv <- SRDat_surv %>% group_by(Stocknumber) %>% summarize(meanLogSurv = mean(log(Surv))) %>% 
  select(meanLogSurv) 
data$MeanLogSurv_surv <- meanLogSurv$meanLogSurv
#data$model <- rep(0,N_Stocks)
#data$Sgen_sig <- TMB_Inputs$Sgen_sig


# Read in wateshed area data and life-history type....
data$WA <- WA$WA
data$Scale <- SRDat_Scale #ordered by std, AR1, surv
##data$Tau_dist <- TMB_Inputs$Tau_dist
# What does inv gamma prior look like? library(invgamma); plot(x=seq(0,1,0.001), y=dinvgamma(seq(0,1,0.001),0.01,0.01), type="l")
data$Stream <- Stream$lh
data$N_stream <-length(which(data$Stream==0))
data$N_ocean <- length(which(data$Stream==1))


#data$logMuDelta1_mean <- TMB_Inputs$logMuDelta1_mean
#data$logMuDelta1_sig <- TMB_Inputs$logMuDelta1_sig
data$logMuDelta2_mean <- TMB_Inputs$logMuDelta2_mean
data$logMuDelta2_sig <- TMB_Inputs$logMuDelta2_sig
#data$Tau_Delta1_dist <- TMB_Inputs$Tau_Delta1_dist
data$Tau_Delta2_dist <- TMB_Inputs$Tau_Delta2_dist
  
# Parameters
param <- list()
Scale.stock_std <- (SRDat %>% group_by(Stocknumber) %>% filter(Stocknumber %not in% c(stksNum_ar,stksNum_surv)) %>% 
                      summarize(Scale.stock_std = max(Scale)))$Scale.stock_std
Scale.stock_ar <- (SRDat %>% group_by(Stocknumber) %>% filter(Stocknumber %in% stksNum_ar) %>% 
                     summarize(Scale.stock_ar = max(Scale)))$Scale.stock_ar
Scale.stock_surv <- (SRDat %>% group_by(Stocknumber) %>% filter(Stocknumber %in% stksNum_surv) %>% 
                     summarize(Scale.stock_surv = max(Scale)))$Scale.stock_surv
#Scale.stock <- 10^(digits$maxDigits-1)

# Parameters for stocks without AR1
param$logA_std <- ( SRDat_std %>% group_by (Stocknumber) %>% summarise(yi = lm(log( Rec / Sp) ~ Sp )$coef[1] ) )$yi
B_std <- SRDat_std %>% group_by(Stocknumber) %>% summarise( m = - lm(log( Rec / Sp) ~ Sp )$coef[2] )
param$logB_std <- log ( 1/ ( (1/B_std$m)/Scale.stock_std ))#log(B_std$m/Scale.stock)
param$logSigma_std <- rep(-2, N_Stocks_std)

# Parameters for stocks with AR1
param$logA_ar <- ( SRDat_ar %>% group_by(Stocknumber) %>% summarise(yi = lm(log( Rec / Sp) ~ Sp )$coef[1] ) )$yi
B_ar <- SRDat_ar %>% group_by(Stocknumber) %>% summarise( m = - lm(log( Rec / Sp) ~ Sp )$coef[2] )
param$logB_ar <- log ( 1/ ( (1/B_ar$m)/Scale.stock_ar ))#Take inverse of B (=Smax and apply scale), the take the inverse again and log to get logB of scaled Smax
param$rho <- rep(TMB_Inputs$rho_Start, N_Stocks_ar)
param$logSigma_ar <- rep (-2, N_Stocks_ar)

# Parameters for stock with survival covariate
param$logA_surv <- ( SRDat_surv %>% group_by(Stocknumber) %>% summarise(yi = lm(log( Rec / Sp) ~ Sp )$coef[1] ) )$yi
B_surv <- SRDat_surv %>% group_by(Stocknumber) %>% summarise( m = - lm(log( Rec / Sp) ~ Sp )$coef[2] )
param$logB_surv <- log ( 1/ ( (1/B_surv$m)/Scale.stock_surv ))#Take inverse of B (=Smax and apply scale), the take the inverse again and log to get logB of scaled Smax
param$logSigma_surv <- rep (-2, N_Stocks_surv)
param$gamma <- rep (0, N_Stocks_surv)

#param$logSgen <- log((SRDat %>% group_by(CU_Name) %>%  summarise(x=quantile(Spawners, 0.5)))$x/Scale) 

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

## Lierman model
#param$logDelta1 <- 10# with skagit 2.881
#param$logDelta1ocean <- 0# with skagit 2.881
#param$logDelta2 <- log(0.72)#log(0.72/(1-0.72)) #logit 0f 0.72 #with skagit logDelta2 = -0.288
#param$logDelta2ocean <- 0#log(0.72/(1-0.72)) #logit 0f 0.72 #with skagit logDelta2 = -0.288
#param$logDeltaSigma <- -0.412 #from Parken et al. 2006 where sig=0.662

## Hierarchcial model

param$logDelta1 <- TMB_Inputs$logDelta1_start#rep(TMB_Inputs$logDelta1_start, 2)
param$logDelta2 <- rep(TMB_Inputs$logDelta2_start, 2)
param$logDeltaSigma <-TMB_Inputs$logDeltaSigma_start 
#param$logMuDelta1 <- TMB_Inputs$logDelta1_start
#param$SigmaDelta1 <- 10
param$logMuDelta2 <- TMB_Inputs$logDelta2_start
param$SigmaDelta2 <- 1



# 3. Estimate SR parameters from synoptic data set and SMSY and SREPs

# Compile model if changed:
#dyn.unload(dynlib("TMB_Files/Ricker_AllMod"))
#compile("TMB_Files/Ricker_AllMod.cpp")

dyn.load(dynlib("TMB_Files/Ricker_AllMod"))

# For Phase 1, fix Delta parameters. Do not need to fix Delta's beccause initlal values are lm fits, so very close
#map <- list(logDelta1=factor(NA), Delta2=factor(NA), logDeltaSigma=factor(NA)) 
#obj <- MakeADFun(data, param, DLL="Ricker_AllMod", silent=TRUE, map=map)

obj <- MakeADFun(data, param, DLL="Ricker_AllMod", silent=TRUE)
#upper <- c(rep(Inf, 80), 5.00, rep(Inf,2))
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))#, upper=upper)
pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1
#summary(sdreport(obj), p.value=TRUE)

#library(tmbstan)
#fitmcmc <- tmbstan(obj, chains=3, iter=1000, init=list(opt$par), control = list(adapt_delta = 0.95))

# Create Table of outputs
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)


# Put together readable data frame of values
All_Ests$Param <- sapply(All_Ests$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))

# By first spliting out stocks modelled with standard Ricker, RickerAR(1), and Ricker-survival models
All_Ests_std <- data.frame()
All_Ests_std <- All_Ests %>% filter (Param %in% c("logA_std", "logB_std", "logSigma_std",  "SMSY_std", "SREP_std"))
SN_std <- unique(SRDat_std[, c("Stocknumber")])
All_Ests_std$Stocknumber <- rep(SN_std)
All_Ests_std <- left_join(All_Ests_std, unique(SRDat_std[, c("Stocknumber", "Name")]))

All_Ests_ar <- data.frame()
All_Ests_ar<- All_Ests %>% filter (Param %in% c("logA_ar", "logB_ar", "rho",  "logSigma_ar", "SMSY_ar", "SREP_ar" ))
SN_ar <- unique(SRDat_ar[, c("Stocknumber")])
All_Ests_ar$Stocknumber <- rep(SN_ar)
All_Ests_ar <- left_join(All_Ests_ar, unique(SRDat_ar[, c("Stocknumber", "Name")]))

All_Ests_surv <- data.frame()
All_Ests_surv<- All_Ests %>% filter (Param %in% c("logA_surv", "logB_surv", "gamma", "logSigma_surv", "SMSY_surv", "SREP_surv" ))#"gamma"
SN_surv <- unique(SRDat_surv[, c("Stocknumber")])
All_Ests_surv$Stocknumber <- rep(SN_surv)
All_Ests_surv <- left_join(All_Ests_surv, unique(SRDat_surv[, c("Stocknumber", "Name")]))

# Combine again
All_Est <- bind_rows(All_Ests_std, All_Ests_ar, All_Ests_surv) 
All_Est$ar <- All_Est$Stocknumber %in% stksNum_ar
All_Est$surv <- All_Est$Stocknumber %in% stksNum_surv
All_Est$Param <- sapply(All_Est$Param, function(x) (unlist(strsplit(x, "[_]"))[[1]]))

# 4. Calculate diagnostics and plot SR curves, etc.

# Calculate AIC

nLL_std <- data.frame(nLL_std=obj$report()$nLL_std) %>% add_column(Stocknumber=SRDat_std$Stocknumber) %>% group_by(Stocknumber) %>% summarize(CnLL=sum(nLL_std))
aic_std <- nLL_std %>% mutate(aic = 2 * 3 + 2*CnLL) # 
nLL_ar <- data.frame(nLL_ar=obj$report()$nLL_ar) %>% add_column(Stocknumber=SRDat_ar$Stocknumber) %>% group_by(Stocknumber) %>% summarize(CnLL=sum(nLL_ar))
aic_ar <- nLL_ar %>% mutate(aic = 2 * 4 + 2*CnLL)
nLL_surv <- data.frame(nLL_surv=obj$report()$nLL_surv) %>% add_column(Stocknumber=SRDat_surv$Stocknumber) %>% group_by(Stocknumber) %>% summarize(CnLL=sum(nLL_surv))
aic_surv <- nLL_surv %>% mutate(aic = 2 * 4 + 2*CnLL)

# Get predicted values and calculate r2
Pred_std <- data.frame()
#Pred_std <- All_Ests %>% filter (Param %in% c("LogR_Pred_std"))
Pred_std <- All_Ests %>% filter (Param %in% c("LogRS_Pred_std"))
#Preds_std <- SRDat_std %>% select("Stocknumber","yr_num", "Rec", "Scale") %>% add_column(Pred=Pred_std$Estimate)
Preds_std <- SRDat_std %>% select("Stocknumber","yr_num", "Sp", "Rec", "Scale", "Name") %>% add_column(Pred=Pred_std$Estimate)
#Preds_std <- Preds_std %>% mutate(ObsLogR = log (Rec / Scale)) 
Preds_std <- Preds_std %>% mutate(ObsLogRS = log ( (Rec / Scale) / (Sp/Scale) ) )
#r2_std <- Preds_std %>% group_by(Stocknumber) %>% summarize(r2=cor(ObsLogR,Pred)^2)
r2_std <- Preds_std %>% group_by(Stocknumber) %>% summarize(r2=cor(ObsLogRS,Pred)^2)

Pred_ar <- data.frame()
#Pred_ar <- All_Ests %>% filter (Param %in% c("LogR_Pred_ar"))
Pred_ar <- All_Ests %>% filter (Param %in% c("LogRS_Pred_ar"))
#Preds_ar <- SRDat_ar %>% select("Stocknumber","yr_num", "Rec", "Scale") %>% add_column(Pred=Pred_ar$Estimate)
Preds_ar <- SRDat_ar %>% select("Stocknumber","yr_num", "Sp", "Rec", "Scale", "Name") %>% add_column(Pred=Pred_ar$Estimate)
#Preds_ar <- Preds_ar %>% mutate(ObsLogR = log (Rec / Scale)) 
Preds_ar <- Preds_ar %>% mutate(ObsLogRS = log ( (Rec / Scale) / (Sp / Scale) ) ) 
#r2_ar <- Preds_ar %>% group_by(Stocknumber) %>% summarize(r2=cor(ObsLogR,Pred)^2)
r2_ar <- Preds_ar %>% group_by(Stocknumber) %>% summarize(r2=cor(ObsLogRS,Pred)^2)

Pred_surv <- data.frame()
#Pred_surv <- All_Ests %>% filter (Param %in% c("LogR_Pred_surv"))
Pred_surv <- All_Ests %>% filter (Param %in% c("LogRS_Pred_surv"))
#Preds_surv <- SRDat_surv %>% select("Stocknumber","yr_num", "Rec", "Scale") %>% add_column(Pred=Pred_surv$Estimate)
Preds_surv <- SRDat_surv %>% select("Stocknumber","yr_num", "Sp", "Rec", "Scale", "Name") %>% add_column(Pred=Pred_surv$Estimate)
#Preds_surv <- Preds_surv %>% mutate(ObsLogR = log (Rec / Scale)) 
Preds_surv <- Preds_surv %>% mutate(ObsLogRS = log ( (Rec / Scale) / (Sp / Scale) ) ) 
#r2_surv <- Preds_surv %>% group_by(Stocknumber) %>% summarize(r2=cor(ObsLogR,Pred)^2)
r2_surv <- Preds_surv %>% group_by(Stocknumber) %>% summarize(r2=cor(ObsLogRS,Pred)^2)

r2 <- bind_rows(r2_std, r2_ar, r2_surv) %>% arrange(Stocknumber)

# Calculate standardized residuals
SRes <- bind_rows(Preds_std, Preds_ar, Preds_surv) %>% arrange (Stocknumber)
SRes <- SRes %>% mutate ( Res = ObsLogRS- Pred) #%>% mutate (StdRes = Res/??)
sigma <- All_Est %>% filter(Param=="logSigma") %>% select(Stocknumber, Estimate, Name)
SRes <- SRes %>% left_join(sigma) %>% rename(logSig = Estimate)
SRes <- SRes %>% mutate (StdRes = Res/exp(logSig))


#Plot SR curves. linearized model, standardized residuals, autocorrleation plots for synoptic data set
if (plot==TRUE){
  PlotSRCurve(SRDat=SRDat, All_Est=All_Est, SMSY_std=SMSY_std, stksNum_ar=stksNum_ar, stksNum_surv=stksNum_surv, r2=r2, removeSkagit=removeSkagit)
  PlotSRLinear(SRDat=SRDat, All_Est=All_Est, SMSY_std=SMSY_std, stksNum_ar=stksNum_ar, stksNum_surv=stksNum_surv, r2=r2, removeSkagit=removeSkagit) 
  PlotStdResid(SRes)
  Plotacf(SRes)
  
}


# What initial values to use for WA model parameters?

SMSY <- All_Est %>% filter(Param=="SMSY") %>% mutate(ModelOrder=0:(length(unique(All_Est$Stocknumber))-1))
# what is scale of SMSY?
Sc <- SRDat %>% select(Stocknumber, Scale) %>% distinct()
SMSY <- SMSY %>% left_join(Sc) %>% mutate(rawSMSY=Estimate*Scale)
lnSMSY <- log(SMSY$rawSMSY)
lnWA <- log(WA$WA)
order <- SMSY %>% select(Stocknumber, ModelOrder)
ParkenSMSY <- as.data.frame(read.csv("DataIn/ParkenSMSY.csv"))
ParkenSMSY <- ParkenSMSY %>% left_join(order) %>% arrange(ModelOrder) %>% mutate(lnSMSY=log(SMSY))
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



