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

# Cowichan modeled with Ricker with a surival co-variate. 
# Harrison was modeled with survival co-variate, but gives very poor fit with very high gamma and 
stksNum_surv <- 23#c(0,23)

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
ind_surv <- add_column(ind_surv, Stocknumber = (unique(SRDat_surv$Stocknumber)))
Surv <- as.data.frame(read.csv("DataIn/Surv.csv")) #%>% filter(Yr<=1999)
SRDat_surv <- SRDat_surv %>% left_join(ind_surv) %>% left_join(Surv)


#remove years 1981-1984, 1986-1987  from Cowichan as per Tompkins et al. 2005
SRDat_surv_Cow <- SRDat_surv %>% filter(Name == "Cowichan" & Yr >= 1985 & Yr !=1986 & Yr != 1987) 
n_surv_Cow <- length(SRDat_surv_Cow$Yr)
SRDat_surv_Cow$yr_num <- 0:(n_surv_Cow-1)
#SRDat_surv_Har <- SRDat_surv %>% filter(Name == "Harrison") 
SRDat_surv <- SRDat_surv_Cow#bind_rows(SRDat_surv_Har, SRDat_surv_Cow)

TMB_Inputs <- list(logA_Start = 2, rho_Start = 0.1, Sgen_sig = 1) #Scale = 1000, 

# Data 
data <- list()
Scale_std <- SRDat_std$Scale 
data$S_std <- SRDat_std$Sp/Scale_std 
data$logR_std <- log(SRDat_std$Rec/Scale_std)
data$stk_std <- as.numeric(SRDat_std$ind_std)
N_Stocks_std <- length(unique(SRDat_std$Name))
data$yr_std <- SRDat_std$yr_num

Scale_ar <- SRDat_ar$Scale 
data$S_ar <- SRDat_ar$Sp/Scale_ar 
data$logR_ar <- log(SRDat_ar$Rec/Scale_ar)
data$stk_ar <- as.numeric(SRDat_ar$ind_ar)
N_Stocks_ar <- length(unique(SRDat_ar$Name))
data$yr_ar <- SRDat_ar$yr_num

Scale_surv <- SRDat_surv$Scale 
data$S_surv <- SRDat_surv$Sp/Scale_surv
data$logR_surv <- log(SRDat_surv$Rec/Scale_surv)
data$stk_surv <- as.numeric(SRDat_surv$ind_surv)
N_Stocks_surv <- length(unique(SRDat_surv$Name))
data$yr_surv <- SRDat_surv$yr_num
data$Surv_surv <- log(SRDat_surv$Surv) #Tompkins et al. used Ln(Surv+1)

#data$model <- rep(0,N_Stocks)
#data$Sgen_sig <- TMB_Inputs$Sgen_sig

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
summary(sdreport(obj), p.value=TRUE)

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

All_Ests_surv <- data.frame()
All_Ests_surv<- All_Ests %>% filter (Param %in% c("logA_surv", "logB_surv", "gamma",  "logSigma_surv", "SMSY_surv", "SREP_surv" ))
SN_surv <- unique(SRDat_surv[, c("Stocknumber")])
All_Ests_surv$Stocknumber <- rep(SN_surv)
All_Ests_surv <- left_join(All_Ests_surv, unique(SRDat_surv[, c("Stocknumber", "Name")]))

# Combine again
All_Est <- bind_rows(All_Ests_std, All_Ests_ar, All_Ests_surv) 
All_Est$ar <- All_Est$Stocknumber %in% stksNum_ar
All_Est$Param <- sapply(All_Est$Param, function(x) (unlist(strsplit(x, "[_]"))[[1]]))


# Get predicted values
Pred_std <- data.frame()
Pred_std <- All_Ests %>% filter (Param %in% c("LogR_Pred_std"))
Preds_std <- SRDat_std %>% select("Stocknumber","yr_num", "Rec", "Scale") %>% add_column(Pred=Pred_std$Estimate)
Preds_std <- Preds_std %>% mutate(ObsLogR = log (Rec / Scale)) 
r2_std <- Preds_std %>% group_by(Stocknumber) %>% summarize(r2=cor(ObsLogR,Pred)^2)

Pred_ar <- data.frame()
Pred_ar <- All_Ests %>% filter (Param %in% c("LogR_Pred_ar"))
Preds_ar <- SRDat_ar %>% select("Stocknumber","yr_num", "Rec", "Scale") %>% add_column(Pred=Pred_ar$Estimate)
Preds_ar <- Preds_ar %>% mutate(ObsLogR = log (Rec / Scale)) 
r2_ar <- Preds_ar %>% group_by(Stocknumber) %>% summarize(r2=cor(ObsLogR,Pred)^2)

Pred_surv <- data.frame()
Pred_surv <- All_Ests %>% filter (Param %in% c("LogR_Pred_surv"))
Preds_surv <- SRDat_surv %>% select("Stocknumber","yr_num", "Rec", "Scale") %>% add_column(Pred=Pred_surv$Estimate)
Preds_surv <- Preds_surv %>% mutate(ObsLogR = log (Rec / Scale)) 
r2_surv <- Preds_surv %>% group_by(Stocknumber) %>% summarize(r2=cor(ObsLogR,Pred)^2)

r2 <- bind_rows(r2_std, r2_ar, r2_surv) %>% arrange(Stocknumber)

# Plot SR curves

# For plotting purposes, need to estimate std Ricker SMSY for AR1 stocks
source ("CheckAR1.r")

Stks <- unique(SRDat$Stocknumber)
NStks <- length(Stks)
par(mfrow=c(5,5), mar=c(3, 2, 2, 1) + 0.1)

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
  SS <- RR<- NA
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
    summarise(SMSY = Estimate * Sc) %>% as.numeric()
  smsy_ul <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="SMSY") %>% 
    summarise(SMSY_ul = Estimate * Sc + 1.96 * Std..Error * Sc ) %>% as.numeric()
  smsy_ll <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="SMSY") %>% 
    summarise(SMSY_ul = Estimate * Sc - 1.96 * Std..Error * Sc ) %>% as.numeric()
  
  
  
  if (i %not in% c(stksNum_ar, stksNum_surv)) abline(v=smsy, col="black") 
  if (i %in% stksNum_ar) abline(v=smsy, col="red") 
  if (i %in% stksNum_surv) abline(v=smsy, col="dark blue") 
  #else abline(v=smsy, col="black")

  if (i %in% stksNum_ar) polygon(x=c(smsy_ul, smsy_ll, smsy_ll, smsy_ul), y=c(-10000,-10000,max(R$Rec),max(R$Rec)), col=rgb(1,0,0, alpha=0.1), border=NA ) 
  if (i %in% stksNum_surv) polygon(x=c(smsy_ul, smsy_ll, smsy_ll, smsy_ul), y=c(-10000,-10000,max(R$Rec),max(R$Rec)), col=rgb(0,0,1, alpha=0.1), border=NA ) 
  if (i %not in% c(stksNum_ar, stksNum_surv))  polygon(x=c(smsy_ul, smsy_ll, smsy_ll, smsy_ul), y=c(-10000,-10000,max(R$Rec),max(R$Rec)), col=grey(0.8, alpha=0.4), border=NA )
  #else polygon(x=c(smsy_ul, smsy_ll, smsy_ll, smsy_ul), y=c(0,0,max(R$Rec),max(R$Rec)), col=grey(0.8, alpha=0.4), border=NA )
  
  if(i %in% stksNum_ar) abline(v=SMSY_std$Estimate[which(SMSY_std$Stocknumber==i)]*Scale.stock[i+1] , col="black")
  if(i %in% stksNum_surv) abline(v=SMSY_std$Estimate[which(SMSY_std$Stocknumber==i)]*Scale.stock[i+1] , col="black")
  
  ParkenSMSY <- as.tibble(read.csv("DataIn/ParkenSMSY.csv"))
  ParkenSMSY <- ParkenSMSY %>% filter(Stocknumber==i) %>% select (SMSY) %>% as.numeric()
  abline(v=ParkenSMSY, lty="dashed")
  lab <-  r2 %>% filter(Stocknumber==i) %>% select(r2) %>% as.numeric() %>% round(3)
  #text(x=max(S$Sp), y= max(R$Rec), labels=paste0("r2=",lab))
  legend("topright", legend = "", title= paste0("r2=",lab), bty="n")
  #legend("topright", legend = "", title= expression(paste(r^2,"=",lab)), bty="n")
  

}

KSR.SMSY <- All_Est %>% filter(Name=="KSR") %>% filter(Param=="SMSY") %>% select(Estimate) %>% as.numeric()
KSR.SMSY.ul <- KSR.SMSY + 1.96 * (All_Est %>% filter(Name=="KSR") %>% filter(Param=="SMSY") %>% select(Std..Error) %>% as.numeric())
KSR.SMSY.ul <- KSR.SMSY.ul * SRDat %>% filter (Name=="KSR") %>% select(Scale) %>% distinct() %>% as.numeric()

Stikine.SMSY <- All_Est %>% filter(Name=="Stikine") %>% filter(Param=="SMSY") %>% select(Estimate) %>% as.numeric()
Stikine.SMSY.ul <- Stikine.SMSY + 1.96 * (All_Est %>% filter(Name=="Stikine") %>% filter(Param=="SMSY") %>% select(Std..Error) %>% as.numeric())
Stikine.SMSY.ul <- Stikine.SMSY.ul * SRDat %>% filter (Name=="Stikine") %>% select(Scale) %>% distinct() %>% as.numeric()

Cow.SMSY <- All_Est %>% filter(Name=="Cowichan") %>% filter(Param=="SMSY") %>% select(Estimate) %>% as.numeric()
Cow.SMSY.ul <- Cow.SMSY + 1.96 * (All_Est %>% filter(Name=="Cowichan") %>% filter(Param=="SMSY") %>% select(Std..Error) %>% as.numeric())
Cow.SMSY.ul <- Cow.SMSY.ul * SRDat %>% filter (Name=="Cowichan") %>% select(Scale) %>% distinct() %>% as.numeric()
Cow.SMSY.ll <- Cow.SMSY - 1.96 * (All_Est %>% filter(Name=="Cowichan") %>% filter(Param=="SMSY") %>% select(Std..Error) %>% as.numeric())
Cow.SMSY.ll <- Cow.SMSY.ll * SRDat %>% filter (Name=="Cowichan") %>% select(Scale) %>% distinct() %>% as.numeric()
