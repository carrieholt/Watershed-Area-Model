
#---------------------------------------------------------
# Libaries

library(rsample)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(TMB)
library(zoo)

#Functions
count.dig <- function(x) {floor(log10(x)) + 1}
'%not in%' <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

#------------------------------------------------------------------------------------------------------------------------------------------------
# Are residuals of Ricker model autocorrelated? Run Ricker_CheckAr1.cpp TMB code below to check

SRDatwNA <- read.csv("DataIn/SRinputfile.csv")
SRDatwNA <- SRDatwNA %>% filter(Name != "Hoko" & Name != "Hoh") #remove two stocks not used in Parken et al, and not documented in Liermann et al.
stockwNA <- SRDatwNA %>% filter (is.na(Rec) == TRUE) %>% select (Stocknumber) %>% unique() %>% unlist() #Do not use AR(1) model on 3  stocks with NAs

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

TMB_Inputs <- list(logA_Start = 2, rho_Start = 0.1, Sgen_sig = 1) #Scale = 1000, 

# Data 
data <- list()
Scale <- SRDat$Scale # Scale <- TMB_Inputs$Scale
data$S <- SRDat$Sp/Scale 
#data$logR <- log(SRDat$Rec/Scale)
data$logRS <- log( (SRDat$Rec/Scale) / (SRDat$Sp/Scale) )
data$stk <- as.numeric(SRDat$Stocknumber)#All stocks are included, so no need to index stock numbers
data$yr <- SRDat$yr_num
N_Stks <- length(unique(SRDat$Name))
data$N_Stks <- N_Stks
#data$model <- rep(0,N_Stocks)


# Parameters
param <- list()
Scale.stock <- 10^(digits$maxDigits-1)
param$logA_std <- ( SRDat %>% group_by (Stocknumber) %>% summarise(yi = lm(log( Rec / Sp) ~ Sp )$coef[1] ) )$yi
B_std <- SRDat %>% group_by(Stocknumber) %>% summarise( m = - lm(log( Rec / Sp) ~ Sp )$coef[2] )
param$logB_std <- log ( 1/ ( (1/B_std$m)/Scale.stock ))#log(B_std$m/Scale.stock)
param$logSigma_std <- rep(-2, N_Stks)


# Compile model if changed:
#dyn.unload(dynlib("TMB_Files/Ricker_CheckAr1"))
#compile("TMB_Files/Ricker_CheckAr1.cpp")

dyn.load(dynlib("TMB_Files/Ricker_CheckAr1"))

obj <- MakeADFun(data, param, DLL="Ricker_CheckAr1", silent=TRUE)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) 
#summary(sdreport(obj))

All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
All_Ests$Param <- sapply(All_Ests$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))


#Preds <- All_Ests %>% filter(Param == "LogR_Pred_std")
Preds <- All_Ests %>% filter(Param == "LogRS_Pred_std")
#Preds <- Preds %>% add_column(yr=data$yr, Stocknumber=data$stk, logR=data$logR) %>% mutate(Resid=Estimate-logR) 
Preds <- Preds %>% add_column(yr=data$yr, Stocknumber=data$stk, logRS=data$logRS, Name=SRDat$Name) %>% mutate(Res=Estimate-logRS) 
ac <- Preds %>% group_by(Stocknumber) %>% summarise (autocorr=acf(Res, plot=F)$acf[2])# provides AR(1) autocorrelation
len <- Preds %>% group_by(Stocknumber) %>% summarise (count=length(Res))
ac.CI <- function(n) {qnorm((1 + 0.95)/2)/sqrt(n)} #95% CI for acf assuming white noise, see https://stackoverflow.com/questions/14266333/extract-confidence-interval-values-from-acf-correlogram
len <- len %>% mutate (CI=ac.CI(count))
ac <- ac %>% left_join(len) %>% left_join(unique(SRDat[, c("Stocknumber", "Name")])) %>% filter(abs(autocorr)>CI)
ac # 6 stocks have significant lag-1 autocorrelation: Chikamin, Keta, Blossom, Situk, Siletz, and Columbia Sp

# Plot acfs: 
if (plot==TRUE) {
  png("DataOut/ACFstd.png", width=7, height=7, units="in", res=500)
  Plotacf(Preds)
  dev.off()
  
}




A_std <- All_Ests %>% filter(Param=="logA_std") %>% add_column(Stocknumber=unique(data$stk)) %>% mutate(A=exp(Estimate))
B_std <- All_Ests %>% filter(Param=="logB_std") %>% add_column(Stocknumber=unique(data$stk)) %>% mutate(B=exp(Estimate)/Scale.stock) 
SMSY_std <- All_Ests %>% filter(Param=="SMSY") %>% add_column(Stocknumber=unique(data$stk)) 
nLL_All_std <- data.frame(nLL=obj$report()$nLL) %>% add_column(Stocknumber=data$stk) %>% group_by(Stocknumber) %>% summarize(CnLL_std=sum(nLL))

aic_All_std <- nLL_All_std %>% mutate(aic_std = 2 * 3 + 2 *CnLL_std) 

#Create data frame for plotting SR curves
All_Est <- All_Ests %>% filter (Param %in% c("logA_std", "logB_std", "logSigma_std",  "SMSY"))
All_Est$ParamShort <- sapply(All_Est$Param, function(x) (unlist(strsplit(x, "_"))[1]))
All_Est <- All_Est %>% select(Estimate, Std..Error, ParamShort) %>% rename(Param=ParamShort)
SN_std <- unique(SRDat[, c("Stocknumber")])
All_Est$Stocknumber <- rep(SN_std)
All_Est <- left_join(All_Est, unique(SRDat[, c("Stocknumber", "Name")]))

# Plot SR curves:

if(plot ==TRUE) {
  png("DataOut/SRstd.png", width=7, height=7, units="in", res=500)
  PlotSRCurve(SRDat, All_Est, SMSY_std, stksNum_ar=NA, stksNum_surv=NA, r2=NA, removeSkagit=FALSE, mod="IWAM_FixedSep") #Mod needed just to get Ricker model types
  dev.off()
  
}


#------------------------------------------------------------------------------------
# WHat is AIC of AR1 model (for those stocks that converge)

SRDatwNA <- read.csv("DataIn/SRinputfile.csv")
SRDatwNA <- SRDatwNA %>% filter(Name != "Hoko" & Name != "Hoh") #remove two stocks not used in Parken et al, and not documented in Liermann et al.
stockwNA <- SRDatwNA %>% filter (is.na(Rec) == TRUE) %>% select (Stocknumber) %>% unique() %>% unlist() #Do not use AR(1) model on 3  stocks with NAs

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


#Choose only those stocks for which AR1 model converges (Not Salcha or Chena)
stks_ar <- c( "Harrison", "Stikine", "Taku", "Unuk", "Chikamin", "Keta", "Blossom", 
              "KSR", "Andrew Cr", "Siuslaw", "Situk", "Siletz", "Niehalem", "Lewis", 
              "Klukshu", "Kitsumkalem", "Columbia Sp", "Salcha","Chena", "Chehalis", 
              "Humptulips", "Queets", "Skagit", "Cowichan", "Quillayute")
stksNum_ar <- c(0:16,19:24)
SRDat_ar <- SRDat %>% filter(Stocknumber %in% stksNum_ar) 

ind_ar <- tibble(ind_ar= 0:(length(unique(SRDat_ar$Name))-1))
ind_ar <- add_column(ind_ar, Stocknumber = (unique(SRDat_ar$Stocknumber)))
SRDat_ar <- SRDat_ar %>% left_join(ind_ar)


# Calculate scale for each stock
digits <- SRDat_ar %>% group_by(Stocknumber) %>% summarize(maxDigits = count.dig(max(Sp)))
SRDat_ar <- left_join(SRDat_ar, digits)
SRDat_ar <- SRDat_ar %>% mutate(Scale = 10^(maxDigits-1))

TMB_Inputs <- list(logA_Start = 2, rho_Start = 0.1, Sgen_sig = 1) #Scale = 1000, 

# Data 
data <- list()
Scale <- SRDat_ar$Scale # Scale <- TMB_Inputs$Scale
data$S_ar <- SRDat_ar$Sp/Scale 
#data$logR <- log(SRDat$Rec/Scale)
data$logRS_ar <- log( (SRDat_ar$Rec/Scale) / (SRDat_ar$Sp/Scale) )
data$stkInd_ar <- as.numeric(SRDat_ar$ind_ar)#as.numeric(SRDat_ar$Stocknumber)
data$yr_ar <- SRDat_ar$yr_num
N_Stks <- length(unique(SRDat_ar$Name))
#data$N_Stks <- N_Stks
#data$model <- rep(0,N_Stocks)


# Parameters
param <- list()
Scale.stock <- 10^(digits$maxDigits-1)
param$logA_ar <- ( SRDat_ar %>% group_by (Stocknumber) %>% summarise(yi = lm(log( Rec / Sp) ~ Sp )$coef[1] ) )$yi
B_ar <- SRDat_ar %>% group_by(Stocknumber) %>% summarise( m = - lm(log( Rec / Sp) ~ Sp )$coef[2] )
param$logB_ar <- log ( 1/ ( (1/B_ar$m)/Scale.stock ))#log(B_std$m/Scale.stock)
param$rho <- rep(TMB_Inputs$rho_Start, N_Stks)
param$logSigma_ar <- rep(-2, N_Stks)


# Compile model if changed:
#dyn.unload(dynlib("TMB_Files/Ricker_ar1"))
#compile("TMB_Files/Ricker_ar1.cpp")

dyn.load(dynlib("TMB_Files/Ricker_ar1"))

obj <- MakeADFun(data, param, DLL="Ricker_ar1", silent=TRUE)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) 
#summary(sdreport(obj))

All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
All_Ests$Param <- sapply(All_Ests$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))
#Preds <- All_Ests %>% filter(Param == "LogR_Pred_ar")
Preds <- All_Ests %>% filter(Param == "LogRS_Pred_ar")
#Preds <- Preds %>% add_column(yr=data$yr, Stocknumber=data$stk, logR=data$logR) %>% mutate(Resid=Estimate-logR) 


#my indexing is wrong here. I need to math up indices 0:22 with stocks 0:24 (missing Chena and Salhca)
#Preds <- Preds %>% add_column(yr=data$yr, Stocknumber=data$stk, logRS=data$logRS) %>% mutate(Resid=Estimate-logRS) 
Preds <- Preds %>% add_column(yr=data$yr, Stocknumber=SRDat_ar$Stocknumber, logRS=data$logRS) %>% mutate(Resid=Estimate-logRS) 
ac <- Preds %>% group_by(Stocknumber) %>% summarise (autocorr=acf(Resid, plot=F)$acf[2])# provides AR(1) autocorrelation
len <- Preds %>% group_by(Stocknumber) %>% summarise (count=length(Resid))
ac.CI <- function(n) {qnorm((1 + 0.95)/2)/sqrt(n)} #95% CI for acf assuming white noise, see https://stackoverflow.com/questions/14266333/extract-confidence-interval-values-from-acf-correlogram
len <- len %>% mutate (CI=ac.CI(count))
ac <- ac %>% left_join(len) %>% left_join(unique(SRDat[, c("Stocknumber", "Name")])) %>% filter(abs(autocorr)>CI)
ac # 0 stocks have significant lag-1 autocorrelation from AR1 model. #6 stocks have significant lag-1 autocorrelation in Ricker resids: Chikamin, Keta, Blossom, Situk, Siletz, and Columbia Sp


#A_ar <- All_Ests %>% filter(Param=="logA_ar") %>% add_column(Stocknumber=unique(data$stk)) %>% mutate(A=exp(Estimate))
#B_ar <- All_Ests %>% filter(Param=="logB_ar") %>% add_column(Stocknumber=unique(data$stk)) %>% mutate(B=exp(Estimate)/Scale.stock) 
#SMSY_ar <- All_Ests %>% filter(Param=="SMSY_ar") %>% add_column(Stocknumber=unique(data$stk)) 
A_ar <- All_Ests %>% filter(Param=="logA_ar") %>% add_column(Stocknumber=unique(SRDat_ar$Stocknumber)) %>% mutate(A=exp(Estimate))
B_ar <- All_Ests %>% filter(Param=="logB_ar") %>% add_column(Stocknumber=unique(SRDat_ar$Stocknumber)) %>% mutate(B=exp(Estimate)/Scale.stock) 
SMSY_ar <- All_Ests %>% filter(Param=="SMSY_ar") %>% add_column(Stocknumber=unique(SRDat_ar$Stocknumber)) 
nLL_All_ar <- data.frame(nLL=obj$report()$nLL) %>% add_column(Stocknumber=SRDat_ar$Stocknumber) %>% group_by(Stocknumber) %>% summarize(CnLL_ar=sum(nLL))

aic_All_ar <- nLL_All_ar %>% mutate(aic_ar = 2 * 4 + 2 *CnLL_ar) 
aic_All_ar %>% full_join(aic_All_std) %>% filter ( aic_ar <= aic_std)
#AIC is lower AR1 model for stocks 3,4,5,6,10,11,16,19,20,23. stksNum_ar <- c(4,5,6,10,11,16)

