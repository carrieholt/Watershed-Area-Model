#---------------------------------------------------------
# Integrated Watershed Area Model

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


#---------------------------------------------------------
# Data

SRDatwNA <- read.csv("DataIn/SRinputfile.csv")
SRDatwNA <- SRDatwNA %>% filter(Name != "Hoko" & Name != "Hoh") #remove two stocks not used in Parken et al, and not documented in Liermann et al.

# Which stocks have NAs?
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


stks_ar <- c("Chikamin", "Keta", "Blossom", "Situk", "Siletz", "Columbia Sp")
stksNum_ar <- c(4,5,6,10,11,16)

SRDat_std <- SRDat %>% filter(Stocknumber %not in% stksNum_ar) 
SRDat_ar <- SRDat %>% filter(Stocknumber %in% stksNum_ar) 

# Assign new stock numbers to each stock so that they are sequential. 
ind_std <- tibble(ind_std= 0:(length(unique(SRDat_std$Name))-1))
ind_std <- add_column(ind_std, Stocknumber = (unique(SRDat_std$Stocknumber)))
SRDat_std <- SRDat_std %>% left_join(ind_std)

ind_ar <- tibble(ind_ar= 0:(length(unique(SRDat_ar$Name))-1))
ind_ar <- add_column(ind_ar, Stocknumber = (unique(SRDat_ar$Stocknumber)))
SRDat_ar <- SRDat_ar %>% left_join(ind_ar)

TMB_Inputs <- list(logA_Start = 2, rho_Start = 0.1, Sgen_sig = 1) #Scale = 1000, 

# Data 
data <- list()
Scale_std <- SRDat_std$Scale # Scale <- TMB_Inputs$Scale
data$S_std <- SRDat_std$Sp/Scale_std 
data$logR_std <- log(SRDat_std$Rec/Scale_std)
data$stk_std <- as.numeric(SRDat_std$ind_std)#as.numeric(SRDat_std$Stocknumber)
N_Stocks_std <- length(unique(SRDat_std$Name))
data$yr_std <- SRDat_std$yr_num

Scale_ar <- SRDat_ar$Scale # Scale <- TMB_Inputs$Scale
data$S_ar <- SRDat_ar$Sp/Scale_ar 
data$logR_ar <- log(SRDat_ar$Rec/Scale_ar)
data$stk_ar <- as.numeric(SRDat_ar$ind_ar)#as.numeric(SRDat_ar$Stocknumber)
N_Stocks_ar <- length(unique(SRDat_ar$Name))
data$yr_ar <- SRDat_ar$yr_num
#data$model <- rep(0,N_Stocks)


#data$model[1] <- 1
#data$Sgen_sig <- TMB_Inputs$Sgen_sig

# Parameters
param <- list()
Scale.stock_std <- (SRDat %>% group_by(Stocknumber) %>% filter(Stocknumber %not in% stksNum_ar) %>% 
                      summarize(Scale.stock_std = max(Scale)))$Scale.stock_std
Scale.stock_ar <- (SRDat %>% group_by(Stocknumber) %>% filter(Stocknumber %in% stksNum_ar) %>% 
                     summarize(Scale.stock_ar = max(Scale)))$Scale.stock_ar
#Scale.stock <- 10^(digits$maxDigits-1)

# Parameters for stocks without AR1
#param$logA_std <- rep(TMB_Inputs$logA_Start, N_Stocks)
#param$logB_std <- log(1/( (SRDat %>% group_by(Stocknumber) %>% summarise(x=quantile(Sp, 0.8)))$x/Scale.stock) )
param$logA_std <- ( SRDat_std %>% group_by (Stocknumber) %>% summarise(yi = lm(log( Rec / Sp) ~ Sp )$coef[1] ) )$yi
B_std <- SRDat_std %>% group_by(Stocknumber) %>% summarise( m = - lm(log( Rec / Sp) ~ Sp )$coef[2] )
param$logB_std <- log ( 1/ ( (1/B_std$m)/Scale.stock_std ))#log(B_std$m/Scale.stock)
param$logSigma_std <- rep(-2, N_Stocks_std)

# Parameters for stocks with AR1
#param$logA_ar <- rep(TMB_Inputs$logA_Start, N_Stocks)
#param$logB_ar <- log(1/( (SRDat %>% group_by(Stocknumber) %>% summarise(x=quantile(Sp, 0.8)))$x/Scale.stock) )
param$logA_ar <- ( SRDat_ar %>% group_by(Stocknumber) %>% summarise(yi = lm(log( Rec / Sp) ~ Sp )$coef[1] ) )$yi
B_ar <- SRDat_ar %>% group_by(Stocknumber) %>% summarise( m = - lm(log( Rec / Sp) ~ Sp )$coef[2] )
param$logB_ar <- log ( 1/ ( (1/B_ar$m)/Scale.stock_ar ))#Take inverse of B (=Smax and apply scale), the take the inverse again and log to get logB of scaled Smax
param$rho <- rep(TMB_Inputs$rho_Start, N_Stocks_ar)
param$logSigma_ar <- rep (-2, N_Stocks_ar)


#param$logSgen <- log((SRDat %>% group_by(CU_Name) %>%  summarise(x=quantile(Spawners, 0.5)))$x/Scale) 

# Compile model if changed:
#dyn.unload(dynlib("TMB_Files/Ricker_ar1"))
#compile("TMB_Files/Ricker_ar1.cpp")

#dyn.load(dynlib("TMB_Files/Ricker_ar1"))

# For Phase 1, fix Sgen
#map <- list(logSgen=factor(rep(NA, N_Stocks))) # Determine which parameters to fix
#obj <- MakeADFun(data, param, DLL="Ricker_ar1", silent=TRUE, map=map)
obj <- MakeADFun(data, param, DLL="Ricker_ar1", silent=TRUE)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) 
summary(sdreport(obj))

# Create Table of outputs
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)


# Put together readable data frame of values
All_Ests$Param <- sapply(All_Ests$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))

# By first spliting out stocks modelled with standard Ricker and those with aR(1)
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

# Combine again
All_Est <- bind_rows(All_Ests_std, All_Ests_ar) 
All_Est$ar <- All_Est$Stocknumber %in% stksNum_ar
All_Est$Param <- sapply(All_Est$Param, function(x) (unlist(strsplit(x, "[_]"))[[1]]))

Stks <- unique(SRDat$Stocknumber)
NStks <- length(Stks)
par(mfrow=c(5,5), mar=c(3, 2, 2, 1) + 0.1)

# Plot SR curves
for (i in Stks){
  R <- SRDat %>% filter (Stocknumber==i) %>% select(Rec) 
  S <- SRDat %>% filter (Stocknumber==i) %>% select(Sp) 
  # what is the scale of Ricker b estimate?
  Sc <- SRDat %>% filter (Stocknumber==i) %>% select(Scale) %>% distinct() %>% as.numeric()
  if(i !=22) plot(x=S$Sp, y=R$Rec, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)), ylim=c(0,max(R$Rec) ) )
  if(i ==22) plot(x=S$Sp, y=R$Rec, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)*3), ylim=c(0,max(R$Rec) ) )
  
  a <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="logA") %>% 
    summarise(A=exp(Estimate)) %>% as.numeric()
  # Divide b by scale
  b <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="logB") %>% 
    summarise(B=exp(Estimate)/Sc) %>% as.numeric()
  
  #Parken values for skagit
  skagit_alpha <- 7.74
  skagit_beta <- 0.0000657
  RR_skagit <- NA
  #RR_std <- NA

  for (j in 1:100){
    if (i!=22) SS[j] <- j*(max(S$Sp)/100)
    if (i==22) SS[j] <- j*(max(S$Sp*3)/100)
    RR[j] <- a * SS[j] * exp(-b * SS[j])
    if(i==22) {RR_skagit[j] <- skagit_alpha * SS[j] * exp(-skagit_beta * SS[j])}
    #if (i %in% stks_ar) {RR_std[j] <- A_std$A[which(A_std$Stocknumber==i)] * SS[j] *  exp(-B_std$B[which(B_std$Stocknumber==i)] * SS[j])}
  }
  lines(x=SS, y=RR, col="black")
  if(i==22) lines(x=SS, y=RR_skagit, lty="dashed")
  name <- All_Est %>% filter (Stocknumber==i) %>% select ("Name") %>% distinct()
  mtext(name$Name, side=3)

  # Plot SMSYs (black for std, red for AR(1), and dashed for Parken et al. 2006)
  smsy <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="SMSY") %>% 
    summarise(SMSY=Estimate*Sc) %>% as.numeric()
  
  if(i %in% stksNum_ar) abline(v=smsy, col="red")
  else abline(v=smsy, col="black")
  if(i %in% stksNum_ar) abline(v=SMSY_std$Estimate[which(SMSY_std$Stocknumber==i)]*Scale.stock[i+1] , col="black")
  
  ParkenSMSY <- as.tibble(read.csv("DataIn/ParkenSMSY.csv"))
  ParkenSMSY <- ParkenSMSY %>% filter(Stocknumber==i) %>% select (SMSY) %>% as.numeric()
  abline(v=ParkenSMSY, lty="dashed")
  

}


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
data$logR <- log(SRDat$Rec/Scale)
data$stk <- as.numeric(SRDat$Stocknumber)
N_Stocks <- length(unique(SRDat$Name))
data$yr <- SRDat$yr_num
data$model <- rep(0,N_Stocks)


# Parameters
param <- list()
Scale.stock <- 10^(digits$maxDigits-1)
param$logA_std <- ( SRDat %>% group_by (Stocknumber) %>% summarise(yi = lm(log( Rec / Sp) ~ Sp )$coef[1] ) )$yi
B_std <- SRDat %>% group_by(Stocknumber) %>% summarise( m = - lm(log( Rec / Sp) ~ Sp )$coef[2] )
param$logB_std <- log ( 1/ ( (1/B_std$m)/Scale.stock ))#log(B_std$m/Scale.stock)
param$logSigma_std <- rep(-2, N_Stocks)


# Compile model if changed:
#dyn.unload(dynlib("TMB_Files/Ricker_CheckAr1"))
#compile("TMB_Files/Ricker_CheckAr1.cpp")

#dyn.load(dynlib("TMB_Files/Ricker_CheckAr1"))

obj <- MakeADFun(data, param, DLL="Ricker_CheckAr1", silent=TRUE)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) 
summary(sdreport(obj))

All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
All_Ests$Param <- sapply(All_Ests$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))
Preds <- All_Ests %>% filter(Param == "LogR_Pred_std")
Preds <- Preds %>% add_column(yr=data$yr, Stocknumber=data$stk, logR=data$logR) %>% mutate(Resid=Estimate-logR) 
ac <- Preds %>% group_by(Stocknumber) %>% summarise (autocorr=acf(Resid, plot=F)$acf[2])# provides AR(1) autocorrelation
len <- Preds %>% group_by(Stocknumber) %>% summarise (count=length(Resid))
ac.CI <- function(n) {qnorm((1 + 0.95)/2)/sqrt(n)} #95% CI for acf assuming white noise, see https://stackoverflow.com/questions/14266333/extract-confidence-interval-values-from-acf-correlogram
len <- len %>% mutate (CI=ac.CI(count))
ac <- ac %>% left_join(len) %>% left_join(unique(SRDat[, c("Stocknumber", "Name")])) %>% filter(abs(autocorr)>CI)
ac # 6 stocks have significant lag-1 autocorrelation: Chikamin, Keta, Blossom, Situk, Siletz, and Columbia Sp


A_std <- All_Ests %>% filter(Param=="logA_std") %>% add_column(Stocknumber=unique(data$stk)) %>% mutate(A=exp(Estimate))
B_std <- All_Ests %>% filter(Param=="logB_std") %>% add_column(Stocknumber=unique(data$stk)) %>% mutate(B=exp(Estimate)/Scale.stock) 
SMSY_std <- All_Ests %>% filter(Param=="SMSY") %>% add_column(Stocknumber=unique(data$stk)) 


