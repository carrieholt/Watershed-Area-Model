#--------------------------------------------------------------------------------------------
# Code to estimate Ricker parameters and Sgen in an integrated TMB model

#--------------------------------------------------------------------------------------------
# Libraries
library(rsample)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(TMB)
library(zoo)

#--------------------------------------------------------------------------------------------
# Get data
SRDat <- read.csv("DataIn/SRDat.csv")

TMB_Inputs <- list(Scale = 1000, logA_Start = 1, rho_Start = 0.1, Sgen_sig = 1)


#Set up data and parameter lists for input into TMB model
Scale <- TMB_Inputs$Scale
data <- list()
data$S <- SRDat$Spawners/Scale 
data$logR <- log(SRDat$Recruits/Scale)
data$stk <- as.numeric(SRDat$CU_ID)
data$yr <- SRDat$yr_num
N_Stocks <- length(unique(SRDat$CU_Name))
data$N_Stks <- N_Stocks
#data$Sgen_sig <- TMB_Inputs$Sgen_sig

param <- list()
param$logA <- rep(TMB_Inputs$logA_Start, N_Stocks)
param$logB <- log(1/( (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/Scale) )
param$logSigma_ar <- rep(-2, N_Stocks)

#param$logSgen <- log((SRDat %>% group_by(CU_Name) %>%  summarise(x=quantile(Spawners, 0.5)))$x/Scale) 
#param$logSgen[3] <- 1.63#

# Compile model if changed:
dyn.unload(dynlib("TMB_Files/Ricker_Sgen"))
compile("TMB_Files/Ricker_Sgen.cpp")

dyn.load(dynlib("TMB_Files/Ricker_Sgen"))

#map <- list(logSgen=factor(rep(NA, N_Stocks))) # Determine which parameters to fix
#obj <- MakeADFun(data, param, DLL="Ricker_ar1", silent=TRUE, map=map)
obj <- MakeADFun(data, param, DLL="Ricker_Sgen", silent=TRUE)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) 
summary(sdreport(obj))

# For Phase 2, pull out SMSY values
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]
upper = c(rep(Inf, 3*N_Stocks), log(SMSYs), Inf, Inf )

pl$logSgen <- log(0.3*SMSYs)
#pl$logSgen[3] <- 1.61

#Phase 2 get Sgen, SMSY etc.

obj <- MakeADFun(data, pl, DLL="Ricker_Sgen", silent=TRUE)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
              upper = upper )


