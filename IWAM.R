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

#---------------------------------------------------------
# Data

SRDatwNA <- read.csv("DataIn/SRinputfile.csv")

# Which stocks have NAs?
stockwNA <- SRDatwNA %>% filter (is.na(Rec) == TRUE) %>% select (Stocknumber) %>% unique() %>% unlist() #Do not use AR(1) model on 3  stocks with NAs

# Remove years with NAs, and stock "Hoko" since n=6 for this stock
SRDat <- SRDatwNA %>% filter(Rec != "NA") %>% filter (Name != "Hoko") %>% filter(Stocknumber <= 16)

# Revise yr_num list where NAs have been removed
test <- SRDat %>% filter(Stocknumber == stockwNA[1]| Stocknumber == stockwNA[2] | Stocknumber == stockwNA[3])
if( max(SRDat$Stocknumber) >= stockwNA[1]) {
  for (i in 1:length(stockwNA)) {
    len <- length (SRDat [which (SRDat$Stocknumber == stockwNA[i]), ]$yr_num) - 1
    SRDat [which (SRDat$Stocknumber == stockwNA[i]), ]$yr_num <- c (0:len)
  }
}
test <- SRDat %>% filter(Stocknumber == stockwNA[1]| Stocknumber == stockwNA[2] | Stocknumber == stockwNA[3])

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
#data$model[1] <- 1
#data$Sgen_sig <- TMB_Inputs$Sgen_sig

# Parameters
param <- list()

# Parameters for stocks with AR1
#param$logA_ar <- rep(TMB_Inputs$logA_Start, N_Stocks)
param$logA_ar <- ( SRDat %>% group_by(Stocknumber) %>% summarise(yi = lm(log( Rec / Sp) ~ Sp )$coef[1] ) )$yi
Scale.stock <- 10^(digits$maxDigits-1)
#param$logB_ar <- log(1/( (SRDat %>% group_by(Stocknumber) %>% summarise(x=quantile(Sp, 0.8)))$x/Scale.stock) )
B_ar <- SRDat %>% group_by(Stocknumber) %>% summarise( m = - lm(log( Rec / Sp) ~ Sp )$coef[2] )
param$logB_ar <- log ( 1/ ( (1/B_ar$m)/Scale.stock ))#Take inverse of B (=Smax and apply scale), the take the inverse again and log to get logB of scaled Smax
param$rho <- rep(TMB_Inputs$rho_Start, N_Stocks)
param$logSigma_ar <- rep (-2, N_Stocks)

# Parameters for stocks without AR1
#param$logA_std <- rep(TMB_Inputs$logA_Start, N_Stocks)
param$logA_std <- ( SRDat %>% group_by (Stocknumber) %>% summarise(yi = lm(log( Rec / Sp) ~ Sp )$coef[1] ) )$yi
#param$logB_std <- log(1/( (SRDat %>% group_by(Stocknumber) %>% summarise(x=quantile(Sp, 0.8)))$x/Scale.stock) )
B_std <- SRDat %>% group_by(Stocknumber) %>% summarise( m = - lm(log( Rec / Sp) ~ Sp )$coef[2] )
param$logB_std <- log ( 1/ ( (1/B_std$m)/Scale.stock ))#log(B_std$m/Scale.stock)
param$logSigma_std <- rep(-2, N_Stocks)
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



