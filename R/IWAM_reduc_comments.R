
#### Introduction -------------------------------------------------------------

# This model estimates biological benchmarks for Chinook populations based on 
# watershed area of the spawning habitat
# The underlying model uses the relationship between watershed area and 
# stock-recruitment parameters for a set of Chinook Populations across the NE 
# (synoptic data set) to derive stock-recruitment parameters (and associated 
# benchmarks) from novel populations from their watershed areas. 
# The model is adapted from Parken et al. (2007) and Liermann et al. (2012)
# This version is taken from IWAM.R (developed by C. Holt, adapted by T. 
# Kitching)

# Integrated Watershed Area Model
# Section Descriptions:
# 1. Read in stock-recruitment data, life-history type, and watershed areas
#   for synoptic data set, and watershed area and life-history type for 
#   additional stocks
#   1.a. Data cleaning
#   1.b. Scale calculation
#   1.c. Setup Watershed area sets
# 2. Create data and parameter lists for TMB
# 3. Run TMB model to estimate Ricker parameters and SMSY & SREP for synoptic 
#   data sets, estimate paraemeters of watershed-area regression, and 
#   estimate SMSY and SREP for additional stocks 
# 4. Compile model outputs
# 5. Calculate diagnostics for SR models in synoptic data set and plot SR 
#   curves, etc.
# 6. Calculate prediction intervals for SMSY and SREP estimates for additional 
#   "test" stocks. These are written to a *.csv file

# This reduced code snippet has removed the main wrapper function- To be added

# Future as a function:
#   - Section 1 e.g. data cleaning and setup not included
#   - Section 1 scale calculation should be embedded as either a separate 
#     function or as an input
#   - Section 2 requires a list of inputs both from the "core" data file
#     and a params list from TMB
#   - Section 3 and 4 are basic
#   - Section 5 should be re-worked so that the ouputs are easily then used in 
#     the already created PlotFunctions.R

#### Libraries -----------------------------------------------------------------

library(rsample)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(reshape2)
library(TMB)
library(zoo)
library(viridis)
library(hrbrthemes)

# Both helperFunctions and PlotFunctions are required.

# Tor- based on your experience with COSEWIC R package, can we remove here::here 
# to avoid problems when using in pkg? And in the meantime, if we're running 
# from *.Proj files, a simple source("R/helperFunctions.R") should work.Correct?
source (here::here("R/helperFunctions.R"))
source(here::here("R/PlotFunctions.R"))

# Consider renaming all model.R scripts to include "mod" or some
  # other suffix/prefix
# To do- rename this to IWAM_model.R and rename IWAM.R to IWAMarchived_model.R


#### Remaining wrapper function objects ----------------------------------------

# Originally part of the main wrapper function stated outright
remove.EnhStocks <- TRUE

#### 1. Read in data -------------------------------------------------

# Our data includes: Stock name, stock number, year, spawners, recruits, 
# stream num., and year num.
# NA's are present in this sample data set and will be removed in the 
# following sections.
SRDatwNA <- read.csv(here::here("DataIn/SRinputfile.csv"))


# * Data Removals and Cleaning ----
# First, remove any unused stocks using filter()
# For e.g., two stocks not used in Parken et al, and not documented in Liermann
SRDatwNA <- SRDatwNA %>% filter(Name != "Hoko" & Name != "Hoh") 


# Determine Which stocks have NAs? Below filter returns only stock numbers.
stockwNA <- SRDatwNA %>% filter (is.na(Rec) == TRUE) %>% 
  dplyr::select (Stocknumber) %>%  unique() %>% unlist() 

# Remove years with NAs
SRDat <- SRDatwNA %>% filter(Rec != "NA") 

# Revise yr_num list where NAs have been removed to be continuous
# Create a test df to check the order of stock numbers by yr_num
  # Within the subset of stocks with NA's identified earlier as stocks 20 and 21
  # test_1 is not a required object for the model. It is only for checking
test_1 <- SRDat %>% filter(Stocknumber == stockwNA[1] | 
                             Stocknumber == stockwNA[2])

# if/for loop to adjust main df (SRDat) to have a continuous year list
if( max(SRDat$Stocknumber) >= stockwNA[1]) { # if the max stock number (24)
    # is greater or equal then the stock's identifed (20), then
  for (i in 1:length(stockwNA)) { # for  stocks identified with NAs (2)
    len <- length (SRDat[which (SRDat$Stocknumber == stockwNA[i]), ]$yr_num) - 1
      # Create a single value object based on the length of:
      # the number of year's - 1
    SRDat [which (SRDat$Stocknumber == stockwNA[i]), ]$yr_num <- c (0:len)
      # re-write the year numbering for  SRDat for the selected stock for a new
      # total length calculated in the line before
  }
}

# Check for the correction to yr_num - wanted to have the time series - 
# consistent = re-indexing the stocks - so that there are 
# no gaps in the. E.g. 0, 1, 2, 3, remove 2 - 0, 1, 3 (now has a break-point)
# test_2 is not a required object for the model
test_2 <- SRDat %>% filter(Stocknumber == stockwNA[1] | 
                             Stocknumber == stockwNA[2])

# **Future update: mask or simulate NAN's in future - COSEWIC example


# * Scale Calculation ----------------------------------------------------------
# Desired scale: 1000 - 0.1 to 100 - responsible for scaling the spawners
# Points of scaling application and removal:
  # - This scaling is APPLIED in section 2. Create data and parameter lists
  # - This scaling is REMOVED when plotting within the plot functions defined
  # in the file PlotSR.R
  # - This scaling is REMOVED for the calculation of predicted values, R2,
  # and Standard Residuals

# Calculate scale for each stock as a tibble (tidyverse df)
digits <- SRDat %>% group_by(Stocknumber) %>% 
  summarize(maxDigits = count.dig(max(Sp)))
  # count.dig() Creates a count [numeric] of the max number of digits 
  # of spawners as digits per stock
  # the function count.dig() can be found in the script: helperFunctions.R

# Join main df with digits by Stocknumber and re-write over SRDat
SRDat <- left_join(SRDat, digits)
# Mutate main df to create a new column: Scale
SRDat <- SRDat %>% mutate(Scale = 10^(maxDigits-1)) # Original Scale
# SRDat <- SRDat %>% mutate(Scale = 10^4) # Alternate Scale
  # using mutate; creates a per stock scale by taking the number of digits - 1,
  # as the exponent on a base 10 log scale


# What is the scale of S, R, SMSY, and SREP data,
# Produces df with two columns: stock number, and scale
SRDat_Scale <- SRDat %>% dplyr::select(Stocknumber, Scale) %>% distinct()
# Creates the obj. SRDat_Scale into a vector of only the scales 
SRDat_Scale <- SRDat_Scale$Scale 

# Remove years 1981-1984, 1986-1987  from Cowichan (Stocknumber 23) as per 
  # Tompkins et al. 2005
SRDat_Cow <- SRDat %>% filter(Name == "Cowichan" & 
                                Yr >= 1985 & 
                                Yr !=1986 & 
                                Yr != 1987) # length = 10
n_Cow <- length(SRDat_Cow$Yr)
SRDat_Cow$yr_num <- 0:(n_Cow-1)
SRDat <- SRDat %>%  filter(Name != "Cowichan") %>% bind_rows(SRDat_Cow) %>%
  arrange(Stocknumber)


# * Read in watershed area data and life-history type --------------------------
  # (stream vs ocean)

WA <- read.csv("DataIn/WatershedArea.csv")
# Create a df of names and corresponding stock numbers to use in joining
names <- SRDat %>% dplyr::select (Stocknumber, Name) %>% distinct()

WA <- WA %>% full_join(names, by="Name") %>% arrange(Stocknumber)

Stream <- SRDat %>% dplyr::select(Stocknumber, Name, Stream) %>% 
  group_by(Stocknumber) %>% 
  summarize(lh=max(Stream)) %>% 
  arrange (Stocknumber)


#### 2. Create data and parameter lists for TMB --------------------------------

# Data list
data <- list()
Scale_std <- SRDat$Scale # Scale enters the TMB data
data$S_std <- SRDat$Sp/Scale_std # Spawners / scale 
data$logRS_std <- log( (SRDat$Rec/SRDat$Scale) / (SRDat$Sp/SRDat$Scale) )
  # logged: scaled recruits / scaled spawners
data$stk_std <- as.numeric(SRDat$Stocknumber) 
N_Stocks_std <- length(unique(SRDat$Name))
data$yr_std <- SRDat$yr_num

# Final remaining if statement for mods
data$logMuAs_mean <- 1.5
data$logMuAs_sig <- 2
data$logMuAo_mean <- 0 #1.5
data$logMuAo_sig <- 2
data$HalfNormMean <- 0 #TMB_Inputs$Tau_sigma
data$HalfNormSig <- 1 #TMB_Inputs$Tau_sigma
data$HalfNormMeanA <- 0 #0.44 #TMB_Inputs$Tau_sigma
data$HalfNormSigA <- 1 #0.5 #TMB_Inputs$Tau_sigma
data$SigRicPriorNorm <- as.numeric(F)
data$SigRicPriorGamma <- as.numeric(T)
data$SigRicPriorCauchy <- as.numeric(F)
data$biasCor <- as.numeric(TRUE)
data$Tau_dist <- 0.1
  
data$sigDelta_mean <- 0.80 # See KFrun.R, #For half-normal use N(0,1)
data$sigDelta_sig <- 0.28 # See KFrun.R,
data$sigNu_mean <- 0.84 # See KFrun.R,
data$sigNu_sig <- 0.275 # See KFrun.R,
data$SigDeltaPriorNorm <- as.numeric(F)
data$SigDeltaPriorGamma <- as.numeric(T)
data$SigDeltaPriorCauchy <- as.numeric(F)
data$Tau_D_dist <- 1

data$TestlnWAo <- read.csv("DataIn/WCVIStocks.csv") %>% 
  mutate (lnWA=log(WA)) %>%
  filter(lh==1) %>% 
  pull(lnWA)
# Add aggregated WAs at inlet level
InletlnWA <- data.frame(read.csv("DataIn/WCVIStocks.csv")) %>% 
  filter(Stock != "Cypre") %>% group_by(Inlet) %>%
  summarize(InletlnWA = log(sum(WA))) %>% filter(Inlet != "San Juan") %>%
  filter(Inlet !="Nitinat")
InletlnWAnoEnh <- data.frame(read.csv("DataIn/WCVIStocks.csv")) %>% 
  filter(Stock != "Cypre") %>% filter(Enh==0) %>%
  group_by(Inlet) %>% summarize(InletlnWA = log(sum(WA))) %>% 
  filter(Inlet != "San Juan") %>%
  filter(Inlet !="Nitinat")
CUlnWA <- data.frame(read.csv("DataIn/WCVIStocks.csv")) %>% 
  filter(Stock != "Cypre") %>% group_by(CU) %>%
  summarize(CUlnWA = log(sum(WA)))
CUlnWAnoEnh <- data.frame(read.csv("DataIn/WCVIStocks.csv")) %>% 
  filter(Stock != "Cypre") %>% filter(Enh==0) %>%
  group_by(CU) %>% summarize(CUlnWA = log(sum(WA)))

if(remove.EnhStocks) data$TestlnWAo <- c(data$TestlnWAo, 
                                         InletlnWAnoEnh$InletlnWA,
                                           CUlnWAnoEnh$CUlnWA)
if(!remove.EnhStocks) data$TestlnWAo <- c(data$TestlnWAo, 
                                          InletlnWA$InletlnWA,
                                            CUlnWA$CUlnWA )

# Read in watershed area data and life-history type and scale
data$WA <- WA$WA
data$Stream <- Stream$lh
data$Scale <- SRDat_Scale # Ordered by Stocknumber

# Read in log(watershed area) for additional stocks
# Predicted lnWA for plottig CIs:
data$PredlnWA <- seq(min(log(WA$WA)), max(log(WA$WA)), 0.1)

# Parameters
param <- list()

# Parameters for stocks without AR1
param$logA_std <- ( SRDat %>% group_by (Stocknumber) %>% 
                      summarise(yi = lm(log( Rec / Sp) ~ Sp )$coef[1] ) )$yi
  # SRDat_std: Rec and Sp are not scaled

B_std <- SRDat %>% group_by(Stocknumber) %>% 
  summarise( m = - lm(log( Rec / Sp) ~ Sp )$coef[2] )
  # *Tor* why the negative here?

param$logB_std <- log ( 1/ ( (1/B_std$m)/data$Scale ))
  # *Carrie* Need to apply the scale to the inverse of Beta, and then re-invert 
  # and log it. This way the initial parameter is log(Scaled beta)
  # logB is scaled
  # in the TMB Ricker model - S is scaled
param$logSigma_std <- rep(-2, N_Stocks_std)

param$logMuAs <- 1.5
param$logMuAo <- 0
param$logSigmaA <- -2

## Liermann model
param$logDelta1 <- 3 
param$logDelta1ocean <- 0 #
param$logDelta2 <- log(0.72) 
param$Delta2ocean <- 0 
param$logDeltaSigma <- -0.412 # from Parken et al. 2006 where sig=0.662

param$logNu1 <- 3
param$logNu1ocean <- 0
param$logNu2 <- log(0.72)
param$Nu2ocean <- 0
param$logNuSigma <- -0.412 #from Parken et al. 2006 where sig=0.66


# 3. Estimate SR parameters from synoptic data set and SMSY and SREPs ----------

# mod remains required for the model call to identify the correct file
# It is easier to reduce the name anyways for ease of calling
mod <- "Liermann_PriorRicSig_PriorDeltaSig" 

# Compile model if changed:
#dyn.unload(dynlib(paste("TMB_Files/", mod, sep="")))
#compile(paste("TMB_Files/", mod, ".cpp", sep=""))
dyn.load(dynlib(paste("TMB_Files/", mod, sep="")))

obj <- MakeADFun(data, param, DLL=mod, silent=TRUE, random = c("logA_std"))

upper<-unlist(obj$par)
upper[1:length(upper)]<- Inf

lower<-unlist(obj$par)
lower[1:length(lower)]<- -Inf


#### RUN THE MODEL ---------------------------------------------------------
# Required objects/inputs
  # obj created from MakeADFun function (TMB) that requires:
    # data,
    # parameters,
    # associated TMB file name (dll)
    # misc. information e.g. tracing, random effects parameters 

opt <- nlminb(obj$par, 
              obj$fn, 
              obj$gr, 
              control = list(eval.max = 1e5, iter.max = 1e5), 
              lower=lower, 
              upper=upper)

pl <- obj$env$parList(opt$par) # Gives the parameter estimates from the model
#summary(sdreport(obj), p.value=TRUE)


#### 4. Compile model outputs --------------------------------------------------
  # *Tor*: Rename estimate objects for clarity
  # Rename All_Ests to  Summary_Ests
  # Rename All_Est to ____

# Create Table of outputs
Summary_Ests <- data.frame(summary(sdreport(obj)))
Summary_Ests$Param <- row.names(Summary_Ests)
# Rename parameter names
Summary_Ests$Param <- sapply(Summary_Ests$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))


All_Est <- data.frame()
# Remove all of the _std objects from parameters - will require removal from
  # TMB param's list in advance
All_Est <- Summary_Ests %>% filter (Param %in% c("logA_std", "logB_std", "logSigma_std",  
                                                  "SMSY_std", "SREP_std"))
SN_std <- unique(SRDat[, c("Stocknumber")])
All_Est$Stocknumber <- rep(SN_std)
All_Est <- left_join(All_Est, unique(SRDat[, c("Stocknumber", "Name")]))

logDeltaSigma <- Summary_Ests %>% filter (Param %in% c("logDeltaSigma")) 
DeltaSigmaUCL <- exp(logDeltaSigma$Estimate + logDeltaSigma$Std..Error*1.96)
DeltaSigmaLCL <- exp(logDeltaSigma$Estimate - logDeltaSigma$Std..Error*1.96) 
DeltaSigma <- exp(logDeltaSigma$Estimate)

# Combine again and rename
# All_Est <- All_Ests_std
All_Est$Param <- sapply(All_Est$Param, function(x) (unlist(strsplit(x, "[_]"))[[1]]))
All_Est <- All_Est %>%left_join(Stream, by="Stocknumber")


All_Deltas <- data.frame()
All_Deltas <- Summary_Ests %>% filter (Param %in% c("logDelta1", "logDelta2","sigma_delta", 
                                                "Delta2_bounded", "logDelta1ocean", 
                                                "logDelta2ocean", "Delta2ocean", "logNu1", 
                                                "logNu2", "sigma_nu", "logNu1ocean", 
                                                "Nu2ocean"))


#### 5. Calculate diagnostics and plot SR curves, etc. -------------------------

# Calculate AIC
  # No RE-SCALING
nLL_std <- data.frame(nLL_std=obj$report()$nLL_std) %>% 
  add_column(Stocknumber=SRDat$Stocknumber) %>% group_by(Stocknumber) %>% 
  summarize(CnLL=sum(nLL_std))
aic_std <- nLL_std %>% mutate(aic = 2 * 3 + 2*CnLL) 

# Get predicted values and calculate r2
Pred_std <- data.frame()
Pred_std <- Summary_Ests %>% filter (Param %in% c("LogRS_Pred_std"))
Preds_std <- SRDat %>% dplyr::select("Stocknumber","yr_num", "Sp", "Rec", "Scale", "Name") %>% 
  add_column(Pred=Pred_std$Estimate)
# mutate the predicted values with Scale 
  # RE-SCALED VALUES
  # These Preds_stds are not used for plotting
Preds_std <- Preds_std %>% mutate(ObsLogRS = log ( (Rec / Scale) / (Sp/Scale) ) )
r2 <- Preds_std %>% group_by(Stocknumber) %>% summarize(r2=cor(ObsLogRS,Pred)^2)


# Get predicted values and their SEs to plot CIs
  # *These are not re-scaled*
  # They are used in the plotting functions and scaled within
PredlnSMSY <- data.frame() 
PredlnSMSY <- Summary_Ests %>% filter (Param %in% c("PredlnSMSY_S", "PredlnSMSY_O", "PredlnSMSY_CI", 
                                                "PredlnSMSYs_CI", "PredlnSMSYo_CI"))
PredlnSREP <- data.frame() 
PredlnSREP <- Summary_Ests %>% filter (Param %in% c("PredlnSREP_S", "PredlnSREP_O", "PredlnSREP_CI", 
                                                "PredlnSREPs_CI", "PredlnSREPo_CI"))

# Calculate standardized residuals
  # These are RE-SCALED values
SRes <- Preds_std %>% arrange (Stocknumber)

SRes <- SRes %>% mutate ( Res = ObsLogRS- Pred) #%>% mutate (StdRes = Res/??)
sigma <- All_Est %>% filter(Param=="logSigma") %>% dplyr::select(Stocknumber, Estimate, Name)
SRes <- SRes %>% left_join(sigma) %>% rename(logSig = Estimate)
SRes <- SRes %>% mutate (StdRes = Res/exp(logSig))


#### * Plot SR Curves ----------------------------------------------------------
# Plot SR curves. linearized model, standardized residuals, autocorrleation plots for synoptic data set
# if using a Liermann model, use SRDat=SRDat_std; otherwise SRDat=SRDat
plot <- TRUE
# Plotted values are RE-SCALED either by plotting function or are already
  # scaled e.g., "SRes"

if (plot==TRUE){
  png(paste("DataOut/SR_", mod, ".png", sep=""), width=7, height=7, units="in", res=500)
  PlotSRCurve(SRDat=SRDat, All_Est=All_Est, r2=r2, removeSkagit = FALSE, mod=mod)
  dev.off()
  png(paste("DataOut/SRLin_", mod, ".png", sep=""), width=7, height=7, units="in", res=1000)
  PlotSRLinear(SRDat=SRDat, All_Est=All_Est, r2=r2, removeSkagit = FALSE)
  dev.off()
  png(paste("DataOut/StdResid_", mod, ".png", sep=""), width=7, height=7, units="in", res=1000)
  PlotStdResid(SRes)
  dev.off()
  png(paste("DataOut/ACF_", mod, ".png", sep=""), width=7, height=7, units="in", res=1000)
  Plotacf(SRes)
  dev.off()
}


#### * Plot WA Regression ------------------------------------------------------
# Plotted values are RE-SCALED within plot func()

if(plot==TRUE){
  png(paste("DataOut/WAregSMSY_", mod, "_wBC.png", sep=""), width=7, height=7, units="in", res=500)
  par(mfrow=c(1,1), mar=c(4, 4, 4, 2) + 0.1)
  title_plot <- "Prior Ricker sigma and prior WA regression sigma"
  #title_plot <- "Separate life-histories: n=17\nFixed-effect yi (logDelta1), \nFixed-effect slope (Delta2)"
  plotWAregressionSMSY (All_Est, All_Deltas, SRDat, Stream, WA, PredlnSMSY, 
                        PredlnWA = data$PredlnWA, title1=title_plot, mod)
  dev.off()
  
  png(paste("DataOut/WAregSREP_", mod, "_wBC.png", sep=""), width=7, height=7, units="in", res=500)
  #png(paste("DataOut/WAreg_Liermann_SepRicA_UniformSigmaAPrior.png", sep=""), width=7, height=7, units="in", res=500)
  par(mfrow=c(1,1), mar=c(4, 4, 4, 2) + 0.1)
  title_plot <- "Prior Ricker sigmas and prior on WA regression sigma"
  #title_plot <- "Separate life-histories: n=17\nFixed-effect yi (logDelta1), \nFixed-effect slope (Delta2)"
  plotWAregressionSREP (All_Est, All_Deltas, SRDat, Stream, WA, PredlnSREP, 
                        PredlnWA = data$PredlnWA, title1=title_plot, mod)
  dev.off()
  #plotWAregression (All_Est, All_Deltas, SRDat, Stream, WA, PredlnSMSY, PredlnWA = data$PredlnWA, 
  # title1="Common, fixed yi (logDelta1), \nRandom slope (Delta2)")
}


#### 6. Calculate prediction intervals for SMSY and SREP for additional stocks ----

# Get predicted values to estimate prediction intervals
  # These values are RE-SCALED to raw estimates during outputting in the TMB code
PredlnSMSY_PI <- data.frame()
PredlnSMSY_PI <- Summary_Ests %>% filter (Param %in% c("PredlnSMSY", "lnSMSY"))
PredlnSREP_PI <- data.frame()
PredlnSREP_PI <- Summary_Ests %>% filter (Param %in% c("PredlnSREP", "lnSREP"))

PredlnSMSY_PI$Stocknumber <- rep(SN_std)
PredlnSREP_PI$Stocknumber <- rep(SN_std)


# To calculate prediction intervals, first get predicted and observed logSMSY 
# and logSREP values for synoptic data set
#   (actually only need observed logSMSY and logSREP values)

#  First need to get the scale for each stock
Scale_PI <- SRDat %>% dplyr::select(Stocknumber, Scale) %>% distinct()
PredlnSMSY_PI <- PredlnSMSY_PI %>% left_join(unique(SRDat[, c("Stocknumber", "Name")])) %>% 
  left_join(Scale_PI)
PredlnSREP_PI <- PredlnSREP_PI %>% left_join(unique(SRDat[, c("Stocknumber", "Name")])) %>% 
  left_join(Scale_PI)

# Then need to separate observed stream vs ocean type

# Plsmsys = predicted log SMSY for stream type
# Plsmsyo = predicted log SMSY for ocean type
# Olsmsys = observed log SMSY for stream type
# Olsmsyo = observed log SMSY for ocean type

Plsmsys <- PredlnSMSY_PI %>% filter(Param=="PredlnSMSY") %>% left_join(Stream) %>% 
  filter(lh == 0) %>% pull(Estimate) #Predicted lnSMSY from WA regression- Stream
Plsmsyo <- PredlnSMSY_PI %>% filter(Param=="PredlnSMSY") %>% left_join(Stream) %>% 
  filter(lh == 1) %>% pull(Estimate) #Predicted lnSMSY from WA regression- ocean
Olsmsys <- PredlnSMSY_PI %>% filter(Param=="lnSMSY") %>% left_join(Stream) %>% 
  filter( lh== 0) %>% pull(Estimate) # "observed" lnSMSY data output from SR models- stream
Olsmsyo <- PredlnSMSY_PI %>% filter(Param=="lnSMSY") %>% left_join(Stream) %>% 
  filter(lh == 1) %>% pull(Estimate) # "observed" lnSMSY data output from SR models- ocean
Plsreps <- PredlnSREP_PI %>% filter(Param=="PredlnSREP") %>% left_join(Stream) %>% 
  filter(lh == 0) %>% pull(Estimate) #Predicted lnSMSY from WA regression- Stream
Plsrepo <- PredlnSREP_PI %>% filter(Param=="PredlnSREP") %>% left_join(Stream) %>% 
  filter(lh == 1) %>% pull(Estimate) #Predicted lnSMSY from WA regression- ocean
Olsreps <- PredlnSREP_PI %>% filter(Param=="lnSREP") %>% left_join(Stream) %>% 
  filter(lh == 0) %>% pull(Estimate) # "observed" lnSMSY data output from SR models- stream
Olsrepo <- PredlnSREP_PI %>% filter(Param=="lnSREP") %>% left_join(Stream) %>% 
  filter(lh == 1) %>% pull(Estimate) # "observed" lnSMSY data output from SR models- ocean


# Get watershed areas for synoptic data set to calculate PIs (stream =s and ocean =o)
WAs <- WA %>% left_join(Stream) %>% filter(lh == 0) %>% pull(WA)
WAo <- WA %>% left_join(Stream) %>% filter(lh == 1) %>% pull(WA)

# Get names of WCVI stocks
sn <- read.csv(here::here("DataIn/WCVIStocks.csv"))
StockNames <- c(as.vector(sn$Stock), as.vector(InletlnWA$Inlet), as.vector(CUlnWA$CU))

# Get Predicted SMSY and SREP values for new "test" WCVI stocks and their Prediction Intervals
TestSMSY <- data.frame() 
TestSMSYs <- data.frame() 
TestSMSYo <- data.frame() 
TestSMSY_SREP <- data.frame() 

TestSMSY <- Summary_Ests %>% filter (Param %in% c("TestlnSMSYo")) %>% add_column(Stock = StockNames)
TestSREP <- Summary_Ests %>% filter (Param %in% c("TestlnSREPo")) %>% add_column(Stock = StockNames)
TestSMSYpull <- TestSMSY %>% pull(Estimate)
TestSREPpull <- TestSREP %>% pull(Estimate)

# Use custom function: PredInt() to estimate prediction intervals 
  # PredInt() defined in helperFunctions.R
TestSMSY_PI <- PredInt(x=log(WAo), y=Olsmsyo, Predy=TestSMSYpull, Newx= data$TestlnWAo)
TestSREP_PI <- PredInt(x=log(WAo), y=Olsrepo, Predy=TestSREPpull, Newx= data$TestlnWAo)

# exp() bounds
TestSMSY <- TestSMSY %>% add_column(LL=exp(TestSMSY_PI$lwr), UL=exp(TestSMSY_PI$upr))
TestSREP <- TestSREP %>% add_column(LL=exp(TestSREP_PI$lwr), UL=exp(TestSREP_PI$upr))

TestSMSY <- TestSMSY %>% mutate (Estimate = exp(Estimate)) %>% dplyr::select(-Std..Error, - Param) %>% 
  add_column(Param = "SMSY")
TestSREP <- TestSREP %>% mutate (Estimate = exp(Estimate)) %>% dplyr::select(-Std..Error, - Param) %>% 
  add_column(Param = "SREP")

WCVISMSY <- TestSMSY %>% mutate(Estimate=round(Estimate, 0), LL=round(LL,0), UL=round(UL,0))
WCVISREP <- TestSREP %>% mutate(Estimate=round(Estimate, 0), LL=round(LL,0), UL=round(UL,0))

WCVISMSY <- WCVISMSY %>% bind_rows(WCVISREP)

# Write SMSY and SREP with PIs to file
if(remove.EnhStocks) write.csv(WCVISMSY, "DataOut/WCVI_SMSY_noEnh_wBC.csv")
if(!remove.EnhStocks) write.csv(WCVISMSY, "DataOut/WCVI_SMSY_wEnh_wBC.csv")

#### End -------------------------------------------------
# This is the end of the original runIWAM() function

