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

#---------------------------------------------------------
# Data

SRDatwNA <- read.csv("DataIn/SRinputfile.csv")
stockwNA <- SRDatwNA %>% filter (is.na(Rec) == TRUE) %>% select (Stocknumber) %>% unique() %>% unlist() #Do not use AR(1) model on 3  stocks with NAs
SRDat <- SRDatwNA %>% filter(Rec != "NA") %>% filter (Name != "Hoko") 

# Remove NAs for stock-recruitment modelling
test <- SRDat %>% filter(Stocknumber == stockwNA[1]| Stocknumber == stockwNA[2] | Stocknumber == stockwNA[3])
for (i in 1:length(stockwNA)) {
  len <- length (SRDat [which (SRDat$Stocknumber == stockwNA[i]), ]$yr_num) - 1
  SRDat [which (SRDat$Stocknumber == stockwNA[i]), ]$yr_num <- c (0:len)
}
test <- SRDat %>% filter(Stocknumber == stockwNA[1]| Stocknumber == stockwNA[2] | Stocknumber == stockwNA[3])


TMB_Inputs <- list(Scale = 1000, logA_Start = 1, rho_Start = 0.1, Sgen_sig = 1)

Scale <- TMB_Inputs$Scale
data <- list()
# Data 
data$S <- SRDat$Sp/Scale 
data$logR <- log(SRDat$Rec/Scale)
data$stk <- as.numeric(SRDat$Stocknumber)
N_Stocks <- length(unique(SRDat$Name))
data$yr <- SRDat$yr_num

data$model <- rep(0,N_Stocks)
#data$Sgen_sig <- TMB_Inputs$Sgen_sig

param <- list()
# Parameters for stocks with AR1
param$logA_ar <- rep(TMB_Inputs$logA_Start, N_Stocks)
param$logB_ar <- log(1/( (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/Scale) )
param$rho <- rep(TMB_Inputs$rho_Start, N_Stocks)
param$logSigma_ar <- rep(-2, N_Stocks)

# Parameters for stocks without AR1
param$logA_std <- rep(TMB_Inputs$logA_Start, N_Stocks)
param$logB_std <- log(1/( (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/Scale) )
param$logSigma_std <- rep(-2, N_Stocks)
#param$logSgen <- log((SRDat %>% group_by(CU_Name) %>%  summarise(x=quantile(Spawners, 0.5)))$x/Scale) 


