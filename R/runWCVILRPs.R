
#-------------------------------------------------------------------------------
# Code to run estimation of LRPs for WCVI
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Libraries and Functions
#-------------------------------------------------------------------------------

list.of.packages <- c("tidyverse", "ggplot2", "gsl", "TMB", "viridis", 
                      "rsample", "gridExtra", "reshape2", "zoo", "hrbrthemes")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(tidyverse)
library(ggplot2)
library(gsl)
library(TMB)
library(viridis)
library(rsample)
library(gridExtra)
library(reshape2)
library(zoo)
library(hrbrthemes)


source("R/helperFunctions.r")
source ("R/PlotSR.r")
source ("R/IWAM.R")
source("R/WCVILRPs.R")
source("R/WCVILRPs_bootstrap.R")

#-------------------------------------------------------------------------------
# Step 1
# Run integrated watershed-area model to get SMSY and SREP estiamtes for WCVI 
# Chinook indicator stocks and inlets, with and without enhancement
runIWAM(remove.EnhStocks = TRUE, removeSkagit = FALSE, 
                    mod = "Liermann_PriorRicSig_PriorDeltaSig", plot = FALSE)
  
runIWAM(remove.EnhStocks = FALSE, removeSkagit = FALSE, 
        mod = "Liermann_PriorRicSig_PriorDeltaSig", plot = FALSE)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Step 2
# Run LRP code to derive inlet-level Sgen and SMU level LRPs for WCVI CK, 
# accounting for uncertainty in LRP from the logistic regression
# Keep Bern_logistic = FALSE for now, as for WCVI CK, there are no years ppn=1
# so Bernoulli logistic regression cannot be esimated

# There are two assumptions about productivity, (1) from life-stage model with 
# expert opinion (W. Luedke pers. comm.), and (2) from "RunReconstruction
# which assumes same harvest rates across wCVI Chinook stocks (D. Dobson 
# pers. comm.)  The default is the life-stage model with expert opinion (1)

Get.LRP(remove.EnhStocks = TRUE, Bern_logistic=FALSE, 
        prod="LifeStageModel")$out$LRP

Get.LRP(remove.EnhStocks = FALSE, Bern_logistic=FALSE, 
        prod="LifeStageModel")$out$LRP
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Step 3 
# Run LRP code to derive inlet-level Sgen and SMU level LRPs for WCVI CK,
# further accounting for uncertainty in the underlying benchmarks by 
# bootstrapping the benchmark estimates from the watershed-area model, and 
# productivity estimates from a plausible range derived from expert opinion
# Function to run bootstrapped uncertainty is here: R/WCVILRPs_boostrap.R

# This code estimates LRPs based on Binomial regression. The Bernoulli 
# regression is not possible for WCVI CK because there are no years of data 
# where all inlets are > lower benchmarks
# See RmdReports/WCVI_LRPs_Binomial.Rmd and WCVI_LRPs_Bernoulli.Rmd

# INSTEAD LOOK AT WCVI_bootstrap.R code, bottom
nBS <- 200 # number trials for bootstrapping 
outBench <- list() 

for (k in 1:nBS) {
  out <- Get.LRP.bs(remove.EnhStocks=TRUE, Bern_logistic=FALSE)
  
  # Save LRPs for each bootstrap
  outLRP <- as.data.frame(out$out$LRP) 
  if(k==1) LRP.bs <- data.frame(fit=outLRP$fit, upr=outLRP$upr, lwr=outLRP$lwr)
  if(k>1) LRP.bs <- add_row(LRP.bs, outLRP)
  
  # Save benchmarks for each boostraap  
  outBench[[k]] <- out$bench
}

# # Is 200 enough trials? Yes
# running.mean <- cumsum(LRP.bs$fit) / seq_along(LRP.bs$fit) 
# plot(running.mean)

# Calculate distribution of overall LRPs by integrating bootstrapped LRP 
# values with uncertainty of each LRP value from TMB
LRP.samples <- rnorm(nBS*10, LRP.bs$fit, (LRP.bs$fit - LRP.bs$lwr) / 1.96)
hist(LRP.samples)
LRP.boot <- quantile(LRP.samples, probs=c(0.05, 0.5, 0.95))
names(LRP.boot) <- c("lwr", "LRP", "upr")


nameStocks <- read.csv("DataOut/WCVI_SMSY_noEnh_wBC.csv") %>% 
  filter(Stock != "Cypre") %>% select(Stock)
nameStocks <- unique(nameStocks)



# Compile bootstrapped estimates of Sgen, SMSY, and SREP, and identify 5th and 
# 95th percentiles
SGEN.bs <- select(as.data.frame(outBench), starts_with("SGEN"))
rownames(SGEN.bs) <- nameStocks$Stock
SGEN.boot <- data.frame(SGEN= apply(SGEN.bs, 1, quantile, 0.5), 
                        lwr=apply(SGEN.bs, 1, quantile, 0.05),
                        upr=apply(SGEN.bs, 1, quantile, 0.95) )

SMSY.bs <- select(as.data.frame(outBench), starts_with("SMSY"))
rownames(SMSY.bs) <- nameStocks$Stock
SMSY.boot <- data.frame(SMSY= apply(SMSY.bs, 1, quantile, 0.05), 
                        lwr=apply(SMSY.bs, 1, quantile, 0.05),
                        upr=apply(SMSY.bs, 1, quantile, 0.95) )

SREP.bs <- select(as.data.frame(outBench), starts_with("SREP"))
rownames(SREP.bs) <- nameStocks$Stock
SREP.boot <- data.frame(SREP= apply(SREP.bs, 1, quantile, 0.5), 
                        lwr=apply(SREP.bs, 1, quantile, 0.05),
                        upr=apply(SREP.bs, 1, quantile, 0.95) )

boot <- list(LRP.boot=LRP.boot, SGEN.boot=SGEN.boot, SMSY.boot=SMSY.boot, SREP.boot=SREP.boot)


#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# Step 4
# Plot outputs
# See RmdReports/WCVI_LRPs_Binomial.Rmd and WCVI_LRPs_Bernoulli.Rmd
#-------------------------------------------------------------------------------

