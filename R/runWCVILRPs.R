
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
Get.LRP(remove.EnhStocks = TRUE, Bern_logistic=FALSE)$out$LRP

Get.LRP(remove.EnhStocks = FALSE, Bern_logistic=FALSE)$out$LRP
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Step 3 (see WCVILRPs_boostrap.R)
# Run LRP code to derive inlet-level Sgen and SMU level LRPs for WCVI CK,
# further accounting for uncertainty in the underlying benchmarks by 
# bootstrapping the benchmark estimates from the watershed-area model, and 
# productivity estimates from a plausible range derived from expert opinion

nBS <- 200 # number trials for bootstrapping
outBench <- list() 

for (k in 1:nBS) {
  out <- Get.LRP.bs(remove.EnhStocks = TRUE, Bern_logistic=FALSE)
  outLRP <- as.data.frame(out$out$LRP) 
  if(k==1) LRP.bs <- data.frame(fit=outLRP$fit, upr=outLRP$upr, lwr=outLRP$lwr)
  if(k>1) LRP.bs <- add_row(LRP.bs, outLRP)
  
  outBench[[k]] <- out$bench
}
# hist(LRP.bs$fit)

# # Is 200 enough trials? Yes
# running.mean <- cumsum(LRP.bs$fit) / seq_along(LRP.bs$fit) 
# plot(running.mean)

# Calculate distribution of overall LRPs by integrating bootstrapped LRP 
# values with uncertainty of each LRP value from TMB
LRP.samples <- rnorm(nBS*10, LRP.bs$fit, (LRP.bs$fit - LRP.bs$lwr) / 1.96)
hist(LRP.samples)
LRP.boot <- quantile(LRP.samples, probs=c(0.05, 0.5, 0.95))
names(LRP.boot) <- c("lwr", "LRP", "upr")

SGEN.bs <- select(as.data.frame(outBench), starts_with("SGEN"))
rownames(SGEN.bs) <- stock_SMSY$Stock
SGEN.boot <- data.frame(SGEN= apply(SGEN.bs, 1, quantile, 0.05), 
                        lwr=apply(SGEN.bs, 1, quantile, 0.5),
                        upr=apply(SGEN.bs, 1, quantile, 0.95) )

SMSY.bs <- select(as.data.frame(outBench), starts_with("SMSY"))
rownames(SMSY.bs) <- stock_SMSY$Stock
SMSY.boot <- data.frame(SMSY= apply(SMSY.bs, 1, quantile, 0.05), 
                        lwr=apply(SMSY.bs, 1, quantile, 0.5),
                        upr=apply(SMSY.bs, 1, quantile, 0.95) )

#Print median and upper and lower 95% intervals for LRP, SGEN & SMSY
LRP.boot
SGEN.boot
SMSY.boot


#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# Step 4
# Plot outputs
# See RmdReports/WCVI_LRPs.Rmd
#-------------------------------------------------------------------------------

