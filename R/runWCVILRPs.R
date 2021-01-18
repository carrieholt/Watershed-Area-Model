
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
# Run LRP code to derive inlet-level Sgen and SMU level LRPs for WCVI CK
Get.LRP(remove.EnhStocks = TRUE)$out$LRP

Get.LRP(remove.EnhStocks = FALSE)$out$LRP
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Step 3
# Plot outputs
# See RmdReports/WCVI_LRPs.Rmd
#-------------------------------------------------------------------------------

