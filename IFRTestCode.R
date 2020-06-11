# Code by C.Holt, Adapted from LRP code developed by K. Holt and B. Davis

# Purposes: 
# (1) estiamte parameters watershed-area model for Chinook salmon, 
# (2) estimate CU-level benchmarks based on watershed-area model
# (3) estimate SMU-level LRPs based on logistic regression

library(MASS) # dose.p function to get SE around P95
library(rsample)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(TMB)
library(zoo)


sourceAll <- function(){
  source("benchmarkFunctions.r")
  source("helperFunctions.r")
}
sourceAll()

# Load TMB  model (code by B. Davis, adapted by C. Holt)

# Compile model if changed:
dyn.unload(dynlib("TMB_Files/Aggregate_LRPs"))
compile("TMB_Files/Aggregate_LRPs.cpp")

dyn.load(dynlib("TMB_Files/Aggregate_LRPs"))

# ======================================================================
# Read-in example data:  
# =====================================================================
EscDat <- read.csv("DataIn/EscDat.csv")
SRDat <- read.csv("DataIn/SRDat.csv")



# ==================================================================================================================
# Set up inputs
# =====================================================================================================================

TMB_Inputs <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 2.5, 
                   logMuA_sig = 2, Tau_dist = 0.1, Tau_A_dist = 0.1, 
                   gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1)


# Choose p values for logistic model
ps <- c(seq(0.5, 0.95,.05), 0.99)
pp <- 1
p <- ps[pp]

# Other inputs 
genYrs <- 3
Mod <- "Aggregate_LRPs"
LRPmodel <- "BinLogistic"
if (LRPmodel == "BernLogistic") Bern_Logistic <- TRUE
if (LRPmodel == "BinLogistic") Bern_Logistic <- FALSE



#Set up data and parameter lists for input into TMB model
Scale <- TMB_Inputs$Scale
data <- list()
data$S <- SRDat$Spawners/Scale 
data$logR <- log(SRDat$Recruits/Scale)
data$stk <- as.numeric(SRDat$CU_ID)
N_Stocks <- length(unique(SRDat$CU_Name))
data$N_Stks <- N_Stocks
data$yr <- SRDat$yr_num
data$Bern_Logistic <- as.numeric(Bern_Logistic)

# also give model year for which will fit logistic model
# Only fit the logistic model to years when all CUs have escapement data

Num_CUs_Over_Time <- EscDat %>%  filter(is.na(Escp) == F) %>% group_by(yr) %>% summarise(n=length(unique((CU_Name))))
Mod_Yrs <- Num_CUs_Over_Time$yr[Num_CUs_Over_Time$n == max(Num_CUs_Over_Time$n)]
Logistic_Dat <- EscDat %>% filter(yr %in% Mod_Yrs)
# need year as index
Logistic_Dat$yr_num <- group_indices(Logistic_Dat, yr) - 1
# use trigger for whether or not use generational means outside
# of TMB to keep TMB model as simple as possible
Agg_Abund <- Logistic_Dat %>% group_by(yr) %>% summarise(Agg_Esc = sum(Escp)) %>%
  mutate(Gen_Mean = rollapply(Agg_Esc, genYrs, gm_mean, fill = NA, align="right"))

data$LM_S <- Logistic_Dat$Escp / Scale
data$LM_Agg_Abund <- Agg_Abund$Agg_Esc / Scale
data$LM_yr <- Logistic_Dat$yr_num
data$LM_stk <- Logistic_Dat$CU_ID

# range of agg abund to predict from
data$Pred_Abund <- seq(0, max(data$LM_Agg_Abund), length.out = 100)
data$p <- p

# set variance to be used for likelihood for estimating Sgen
data$Sgen_sig <- TMB_Inputs$Sgen_sig

param <- list()
param$logA <- rep(TMB_Inputs$logA_Start, N_Stocks)
param$logB <- log(1/( (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/Scale) )
param$logSigma <- rep(-2, N_Stocks)
param$logSgen <-  log((SRDat %>% group_by(CU_Name) %>%  summarise(x=quantile(Spawners, 0.5)))$x/Scale) 
param$B_0 <- 2
param$B_1 <- 0.1



# ==================================================================================================================
# Run TMB
# =====================================================================================================================

# Phase 1 estimate SR params
map <- list(logSgen=factor(rep(NA, data$N_Stks)), B_0=factor(NA), B_1=factor(NA)) # Determine which parameters to fix

obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, map=map)


opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1
summary(sdreport(obj))


# Are residuals of Ricker model autocorrelated?


# pull out SMSY values
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]
upper = c(rep(Inf, 3*N_Stocks), log(SMSYs), Inf, Inf )

pl$logSgen <- log(0.3*SMSYs)

#Phase 2 get Sgen, SMSY etc.
map2 <- list(B_0=factor(NA), B_1=factor(NA))

obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE, map=map2)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
              upper = upper )
pl2 <- obj$env$parList(opt$par) # Parameter estimate after phase 2
summary(sdreport(obj))

# Phase 3 fit logistic model, holding  other estimates constant
obj <- MakeADFun(data, pl2, DLL=Mod, silent=TRUE)

opt <- tryCatch(
  {nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
          upper = upper)},
  error=function(cond) {
    # Choose a return value in case of error
    return(NA)
  },
  warning=function(cond) {
    return(NA)
  }
)   


# ==================================================================================================================
# Compile outputs
# =====================================================================================================================

# Create Table of outputs
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
# put together readable data frame of values
All_Ests$Param <- sapply(All_Ests$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))
All_Ests$Mod <- Mod
All_Ests$CU_ID[!(All_Ests$Param %in% c("logMuA", "logSigmaA", "B_0", "B_1", "Agg_LRP", "gamma", "logsigma", "prod"))] <- rep(0:(N_Stocks-1)) 
All_Ests <- left_join(All_Ests, unique(SRDat[, c("CU_ID", "CU_Name")]))
# don't want logged param values
# need to unlog A,B,Sigma -- took A off this list for now
All_Ests$Estimate[All_Ests$Param %in% c( "logB", "logSigma")] <- exp(All_Ests$Estimate[All_Ests$Param %in% c( "logB", "logSigma")] )
#All_Ests$Param[All_Ests$Param == "logA"] <- "A"
All_Ests$Param[All_Ests$Param == "logB"] <- "B"
All_Ests$Param[All_Ests$Param == "logSigma"] <- "sigma"
All_Ests$Param[All_Ests$Param == "SMSY"] <- "Smsy"
All_Ests[All_Ests$Param == "B",] <- All_Ests %>% filter(Param == "B") %>% mutate(Estimate = Estimate/Scale) %>% mutate(Std..Error = Std..Error/Scale)
All_Ests[All_Ests$Param %in% c("Sgen", "SMSY", "Agg_LRP"), ] <-  All_Ests %>% filter(Param %in% c("Sgen", "SMSY", "Agg_LRP")) %>% 
  mutate(Estimate = Estimate*Scale) %>% mutate(Std..Error = Std..Error*Scale)
Preds <- All_Ests %>% filter(Param == "Logit_Preds")
All_Ests <- All_Ests %>% filter(!(Param %in% c( "logSgen", "Logit_Preds"))) 

out <- list()
out$All_Ests <- All_Ests

if(Bern_Logistic == T){
  N_CUs <- obj$report()$All_Above_BM
} else {
  N_CUs <- obj$report()$N_Above_BM/N_Stocks
}

Logistic_Data <- data.frame(Mod = Mod, yr = Agg_Abund$yr, 
                            yy = N_CUs, xx = Agg_Abund$Agg_Esc)

out$Logistic_Data <- Logistic_Data

Logistic_Fits <- data.frame(xx = data$Pred_Abund*Scale, fit = inv_logit(Preds$Estimate),
                            lwr = inv_logit(Preds$Estimate - 1.96*Preds$Std..Error),
                            upr = inv_logit(Preds$Estimate + 1.96*Preds$Std..Error))

out$Preds <- Logistic_Fits

out$LRP <- data.frame(fit = All_Ests %>% filter(Param == "Agg_LRP") %>% pull(Estimate), 
                      lwr = All_Ests %>% filter(Param == "Agg_LRP") %>% mutate(xx =Estimate - 1.96*Std..Error) %>% pull(xx),
                      upr = All_Ests %>% filter(Param == "Agg_LRP") %>% mutate(xx =Estimate + 1.96*Std..Error) %>% pull(xx))

out



