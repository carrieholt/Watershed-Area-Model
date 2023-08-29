
#### Introduction ------------------------------------------------------

# This version is a selective copy of IWAM_REDUC.R under work by Tor Kitching

# Integrated Watershed Area Model
# Steps:
# 1. Read in stock-recruitment data, life-history type, and watershed areas
#   for synoptic survey data, and watershed area and life-history type for 
#   additional stocks
# 2. Create data and parameter lists for TMB
# 3. Run TMB model to estimate Ricker parameters and SMSY & SREP for synoptic 
#   data sets, estimate paraemeters of watershed-area regression, and 
#   estimate SMSY and SREP for additional stocks 
# 4. Compile model outputs
# 5. Calculate diagnostics for SR models in synoptic data set and plot SR 
#   curves, etc.
# 6. Calculate prediction intervals for SMSY and SREP estimates for additional 
#   "test" stocks. These are written to a *.csv file

# This reduced code snippet has removed the main wrapper function, please
# refer to IWAM.R for the original complete code repository

#### Set Working Directory -----------------------------------------------------

# Tor's Home Computer
# setwd("~/GitHub/Watershed-Area-Model")
# No longer used under here::here conventions

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

# Both helperFunctions and PlotSR are required to be run upon init.
source (here::here("R/helperFunctions.R"))
  # rename to AccessoryFunctions.R

source(here::here("R/PlotFunction.R"))

# Consider renaming all model.R scripts to include "mod" or some
# other suffix/prefix


#### Opening calls -------------------------------------------------------------

# Otherwise part of the main wrapper function, best to outright state these early
removeSkagit <- FALSE 
remove.EnhStocks <- TRUE

# This is where the original runIWAM function begins

#### 1. Read in data -------------------------------------------------

# Our test data includes: Stock name, stock number, year, spawners, recruits, 
  # stream num., and year num.
# NA's are present in this sample data set and will be removed in the 
  # following sections.
SRDatwNA <- read.csv(here::here("DataIn/SRinputfile.csv"))


# * Data Removals and Cleaning ----
# First, remove any unused stocks using filter()

# For e.g., two stocks not used in Parken et al, and not documented in Liermann
SRDatwNA <- SRDatwNA %>% filter(Name != "Hoko" & Name != "Hoh") 

# Remove Skagit and re-align the stock numbers to be continuous post removal
if (removeSkagit==TRUE) {
  SRDatwNA <- SRDatwNA %>% filter(Name != "Skagit")
  # Skagit is Stocknumber = 22. Need to re-align stock numbers of last 2 
    # stocks: 23 & 24
  SRDatwNA [which(SRDatwNA$Stocknumber==23),2] = 22
  SRDatwNA [which(SRDatwNA$Stocknumber==24),2] = 23
}

# Determine Which stocks have NAs? Below filter returns only stock numbers.
stockwNA <- SRDatwNA %>% filter (is.na(Rec) == TRUE) %>% 
  dplyr::select (Stocknumber) %>%  unique() %>% unlist() 
# Do not use AR(1) model on  stocks with NAs, Humptulips and Queets (20 & 21)

# Remove years with NAs
SRDat <- SRDatwNA %>% filter(Rec != "NA") 

# Revise yr_num list where NAs have been removed to be continuous
  # yr_num is used as a method for indexing within [**INSERT**]
# Create a test df to check the order of stock numbers by yr_num
  # Within the subset of stocks with NA's identified earlier as stocks 20 and 21
  # test_1 is not a required object for the model. It is only for checking
test_1 <- SRDat %>% filter(Stocknumber == stockwNA[1]|Stocknumber == stockwNA[2])

# if/for loop to adjust main df (SRDat) to have a continuous year list
if( max(SRDat$Stocknumber) >= stockwNA[1]) { # if the max stock number (24)
    # is greater or equal then the stock's identifed (20), then
  for (i in 1:length(stockwNA)) { # for the number of stocks identifed with NAs (2)
    len <- length (SRDat [which (SRDat$Stocknumber == stockwNA[i]), ]$yr_num) - 1
      # create a single value object based on the length of:
      # the number of year's - 1
    SRDat [which (SRDat$Stocknumber == stockwNA[i]), ]$yr_num <- c (0:len)
      # re-write the year numbering for  SRDat for the selected stock for a new
      # total length calculated in the line before
  }
}

# Check for the correction to yr_num - wanted to have the time series - 
# consistent = re-indexing the stocks - so that there are 
# no gaps in the. E.g. 0, 1, 2, 3, remove 2 - 0, 1, 3 (now has a break-point)
# test_2 is not a required object for the model. It is only for checking
test_2 <- SRDat %>% filter(Stocknumber == stockwNA[1]| Stocknumber == stockwNA[2])

# Future update: mask or simulate NAN's in future - COSEWIC example


# * Scale Calculation ----------------------------------------------------------

# *Tor*:
  # Scaling in this manner has the following problems:
    # 1. It is skewed in the presence of outliers.
    # 2. It is not a form of standarization. Different stocks with the same -
      # number of digits may vary in distribution and lead to distortion
      # in the scaled values.
    # 3. Dividing by base 10 may lead to a loss in information.
  # Other options to scale - with a targeted range in mind:
    # 1. Min-max scaling (Issues with knowing true min and max ie. clipping)
    # 2. z-score normalization (Can affect interpretation if you want it -
      # in the original units.)
  # Original notes: hoping for around 1000 - 0.1 to 100 - responsible 
    # for scaling the spawners

  # New option: scale with a simple scalar across all populations
    # Next steps:
      # - If we go by digits, the unique options are: 6, 5, 4, 3
      # This will mean 10^ 5, 4, 3, 2
      # and a max range of 10^3?
# Points of scaling application and removal:
  # - This scaling is APPLIED in section 2. Create data and parameter lists
  # - This scaling is REMOVED when plotting within the plot functions defined
  # in the file PlotSR.R
  # *Tor* I recommend we rename the PlotSR file to make sure its importance
  # is recognized in the repository.
  # - This scaling is REMOVED for the calculation of predicted values, R2,
  # and Standard Residuals

# Calculate scale for each stock as a tibble (tidyverse df)
digits <- SRDat %>% group_by(Stocknumber) %>% 
  summarize(maxDigits = count.dig(max(Sp)))
  # count.dig() Creates a count [numeric] of the max number of digits of spawners as 
    # digits per stock
  # the function count.dig() can be found in the script: helperFunctions.R

# join main df with digits by Stocknumber and re-write over SRDat
SRDat <- left_join(SRDat, digits)
# mutate main df to create a new column: Scale
SRDat <- SRDat %>% mutate(Scale = 10^(maxDigits-1)) # Original Scale
# SRDat <- SRDat %>% mutate(Scale = 10^4) # Alternate Scale
  # using mutate; creates a per stock scale by taking the number of digits - 1, as
    # the exponent on a base 10 log scale

# *Tor*: Working on scale
  # First: are we really looking to "Scale" or to "Normalize"
  # I think our goal is to get Spawners to be gaussian, in which cause a log-
  # transform works great
  # Check out the three plots below:

# p1 <- ggplot(SRDat, aes(Sp)) +
#   geom_histogram() +
#   labs(title='Raw Counts of Spawners')
# 
# p2 <- ggplot(SRDat, aes(Sp/Scale))+
#   geom_histogram()+
#   labs(title='Scaled Spawners: Sp/(Scale per stock)')
# 
# p3 <- ggplot(SRDat, aes(log(Sp)))+
#   geom_histogram()+
#   labs(title='Log(Spawners)')
# 
# p4 <- ggplot(SRDat, aes(Sp/10^4))+
#   geom_histogram()+
#   labs(title='Spawners per 10^5')
# 
# plot_grid(p2, p4) # grid of above plots
# plot(p4)
# 
# testing <- SRDat$Sp/10^4
# original <- SRDat$Sp/SRDat$Scale


# * Calculation of Survival Covariates -----------------------------------------
# Pre-calculation of survival covariates 
  # previously differentiated for standard, AR(1), survival
# If using mod=="Liermann_PriorRicSig_PriorDeltaSig" SRDat will overwrite
  # everything done below

# Cowichan, stk-23, not included here because they are modeled as per Tompkins
  # with a survival covariate
# Create a vector that will be used to identify stocks per model type going forward
  # Required for distinguishing stocks when plotting SR curves
# stksNum_ar <- c(4,5,6,10,11,16) # Ricker with a surival co-variate. 
# stksNum_surv <- c(0,23)
# stks_surv <- c("Harrison", "Cowichan") # Stock name

# Remove Skagit
# if (removeSkagit==TRUE) {stksNum_surv <- c(0,22)}
  # only works if removeSkagit is found
  # Changes stksNum_surv - which is no longer used anyways
  # Once _surv is fully removed this line can be removed

# Creates an object with the length of the number of stocks
# len_stk <- length(unique(SRDat$Stocknumber))
# Creates a vector of the stock ID's of those not included in the AR or Survival
  # model (_ar, and _surv)
# stksNum_std <- which(0:(len_stk-1) %not in%c(stksNum_ar, stksNum_surv)==TRUE)-1 
  # Assuming there are only 25 stocks (0:24 StockNumber)
  # Required for creating the stock order object used to arrange models 
    # when the stream and ocean life-history data is added prior to parameter
    # creation (Section 2)

# When aggregated standard, ar1, surv, "ModelOrder" is the order of stocks for 
  # aligning with WA and life-history data
# stksOrder <- data.frame(Stocknumber =  c(stksNum_std, stksNum_ar, stksNum_surv),
#                         ModelOrder = 0:(len_stk-1))
# stksOrder: is used to arrange the stocks when the stream and ocean
  # life-history data is added into the complete df.

# What is the scale of S, R, SMSY, and SREP data,
  # ordered when aggregated by std, AR1, surv?
# Produces df with two columns: stock number, and scale
SRDat_Scale <- SRDat %>% dplyr::select(Stocknumber, Scale) %>% distinct()

# Missing a loop file - check if ModelOrder is required - if so -this is where
  # it would be missing from.

SRDat_Scale <- SRDat_Scale$Scale 
# Creates the obj. SRDat_Scale into a vector of only the scales - now
  # ordered by stock number
  # Not sure if ModelOrder matters in this case

# If using any of the below mods:
  # overwrite SRDat_std with the main df SRDat
  # Othewise, use standard Ricker model
# SRDat_std <- SRDat
  # Now can remove all _std from SRDat?

# Assign new stock numbers to each stock so that they are sequential within each 
  # model form. These are used in TMB When only one model form is use (std Ricker)
  # then ind_std = Stocknumber
# ind_std <- tibble(ind_std= 0:(length(unique(SRDat_std$Name))-1))
# ind_std <- add_column(ind_std, Stocknumber = (unique(SRDat_std$Stocknumber)))
  # this is just two columns with the same numbers?
# SRDat_std <- SRDat_std %>% left_join(ind_std)
  # ind_std is now a duplicate column of stock numbers

# Remove years 1981-1984, 1986-1987  from Cowichan (Stocknumber 23) as per 
  # Tompkins et al. 2005
SRDat_Cow <- SRDat %>% filter(Name == "Cowichan" & Yr >= 1985 & Yr !=1986 & Yr != 1987) # length = 10
n_Cow <- length(SRDat_Cow$Yr)
  #SRDat_Cow$yr_num <- 0:(n_surv_Cow-1) #n_surv_Cow can most likely be exchanged for n_Cow
SRDat_Cow$yr_num <- 0:(n_Cow-1)
SRDat <- SRDat %>%  filter(Name != "Cowichan") %>% bind_rows(SRDat_Cow) %>%
    # arrange(ind_std) 
    arrange(Stocknumber)
  # changed from arrange(ind_std) to arrange(Stocknumber)


# * Read in watershed area data and life-history type --------------------------
  # (stream vs ocean)

WA <- read.csv("DataIn/WatershedArea.csv")
# Create a df of names and corresponding stock numbers to use in joining
names <- SRDat %>% dplyr::select (Stocknumber, Name) %>% distinct()

WA <- WA %>% full_join(names, by="Name") %>% arrange(Stocknumber)
  # Only difference between WA and test_WA is the ModelOrder
  # Checking now for ModelOrder future usage
    # Previously used under the section: Additional code; not currently needed
  # Final decision: Tor: Removal

if (removeSkagit==TRUE) {WA <- WA %>% filter(Name !="Skagit")}

Stream <- SRDat %>% dplyr::select(Stocknumber, Name, Stream) %>% group_by(Stocknumber) %>% 
  summarize(lh=max(Stream)) %>% arrange (Stocknumber)
# Stream <- Stream  %>% arrange(Stocknumber) # Added to above code in single pipe for cleanliness


#### 2. Create data and parameter lists for TMB --------------------------------

# Only rho_Start used. Disregard other TMB_Input
  # TMB_Inputs is not used anywhere in the code regarding the model
  # Tor: Can this be deleted?
# TMB_Inputs <- list(rho_Start = 0.0, logDelta1_start=3.00, logDelta2_start =log(0.72), 
#                    logDeltaSigma_start = -0.412, logMuDelta1_mean= 5, logMuDelta1_sig= 10, 
#                    logMuDelta2_mean=-0.5, logMuDelta2_sig= 10, Tau_Delta1_dist= 0.1, 
#                    Tau_Delta2_dist= 0.1, Tau_sigma = 0.01)

# Data list
data <- list()
Scale_std <- SRDat$Scale # Scale enters the TMB data
data$S_std <- SRDat$Sp/Scale_std # Spawners / scale - check Scale calculation
  # This S_std is not used again in this code
data$logRS_std <- log( (SRDat$Rec/Scale_std) / (SRDat$Sp/Scale_std) )
  # logged: scaled recruits / scaled spawners
# data$stk_std <- as.numeric(SRDat$ind_std)
data$stk_std <- as.numeric(SRDat$Stocknumber) # ind_std and Stocknumber are the same
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

data$TestlnWAo <- read.csv("DataIn/WCVIStocks.csv") %>% mutate (lnWA=log(WA)) %>%
  filter(lh==1) %>% pull(lnWA)
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

if(remove.EnhStocks) data$TestlnWAo <- c(data$TestlnWAo, InletlnWAnoEnh$InletlnWA,
                                           CUlnWAnoEnh$CUlnWA)
if(!remove.EnhStocks) data$TestlnWAo <- c(data$TestlnWAo, InletlnWA$InletlnWA,
                                            CUlnWA$CUlnWA )

# Read in wateshed area data and life-history type and scale
data$WA <- WA$WA
data$Stream <- Stream$lh
data$Scale <- SRDat_Scale #ordered by std, AR1, surv, if all 3 Ricker models uses. 
# Otherwise ordered by Stocknumber

# Read in log(watershed area) for additional stocks
# Predicated lnWA for plottig CIs:
data$PredlnWA <- seq(min(log(WA$WA)), max(log(WA$WA)), 0.1)

# Parameters
param <- list()

Scale.stock_std <- (SRDat %>% group_by(Stocknumber) %>% 
  summarize(Scale.stock_std = max(Scale)))$Scale.stock_std
  # SRDat - main data
  # Groupby Stock Number
  # Filter for specific Stocks
  # summarize into new df the max Scale
    # Scale has always been max per Stock

# Parameters for stocks without AR1
param$logA_std <- ( SRDat %>% group_by (Stocknumber) %>% 
                      summarise(yi = lm(log( Rec / Sp) ~ Sp )$coef[1] ) )$yi
  # SRDat_std: Rec and Sp are not scaled

B_std <- SRDat %>% group_by(Stocknumber) %>% 
  summarise( m = - lm(log( Rec / Sp) ~ Sp )$coef[2] )
  # *Tor* why the negative here?

param$logB_std <- log ( 1/ ( (1/B_std$m)/Scale.stock_std )) # log(B_std$m/Scale.stock)
  # *Carrie* Need to apply the scale to the inverse of Beta, and then re-invert and log it.
    # This way the initial parameter is log(Scaled beta)
  # logB is scaled
  # in the TMB Ricker model - S is scaled
param$logSigma_std <- rep(-2, N_Stocks_std)

param$logMuAs <- 1.5
  #param$logSigmaAs <- 1
param$logMuAo <- 0#1.5
param$logSigmaA <- -2#5
  #param$logSigmaAo <- 1

## Separate Stream and Ocean type models
#param$slogDelta1 <- 2.744 # best estimates from run of stream-specific WAreg TMB
#param$sDelta2 <- 0.857 
#param$slogDeltaSigma <- -0.709 

#param$ologDelta1 <- 3.00 # 1.519 # best estimates from run of ocean-specific WAreg TMB 
#param$ologDelta2 <- log(0.94) # 0 # 21.2 
#param$ologDeltaSigma <-  -0.412 #-0.94 

## Lierman model
param$logDelta1 <- 3 #10 # with skagit 2.881
param$logDelta1ocean <- 0 # with skagit 2.881
param$logDelta2 <- log(0.72) #log(0.72/(1-0.72)) 
param$Delta2ocean <- 0 #log(0.72/(1-0.72)) 
param$logDeltaSigma <- -0.412 # from Parken et al. 2006 where sig=0.662

param$logNu1 <- 3#10# with skagit 2.881
param$logNu1ocean <- 0# with skagit 2.881
param$logNu2 <- log(0.72)#log(0.72/(1-0.72)) 
param$Nu2ocean <- 0#log(0.72/(1-0.72)) 
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
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5), 
              lower=lower, upper=upper)
pl <- obj$env$parList(opt$par) # Gives the parameter estimates from the model
#summary(sdreport(obj), p.value=TRUE)


#library(tmbstan)
# fitmcmc <- tmbstan(obj, chains=3, iter=1000, init=list(opt$par), 
# control = list(adapt_delta = 0.95))
# fitmcmc <- tmbstan(obj, chains=3, iter=1000, init=list("last.par.best"), 
# control = list(adapt_delta = 0.95))


#### 4. Compile model outputs --------------------------------------------------

# Create Table of outputs
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)


# Put together readable data frame of values
All_Ests$Param <- sapply(All_Ests$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))

# By first spliting out stocks modelled with standard Ricker, RickerAR(1), 
# and Ricker-survival models
All_Ests_std <- data.frame()
All_Ests_std <- All_Ests %>% filter (Param %in% c("logA_std", "logB_std", "logSigma_std",  
                                                  "SMSY_std", "SREP_std"))
SN_std <- unique(SRDat[, c("Stocknumber")])
All_Ests_std$Stocknumber <- rep(SN_std)
All_Ests_std <- left_join(All_Ests_std, unique(SRDat[, c("Stocknumber", "Name")]))

logDeltaSigma <- All_Ests %>% filter (Param %in% c("logDeltaSigma")) 
DeltaSigmaUCL <- exp(logDeltaSigma$Estimate + logDeltaSigma$Std..Error*1.96)
DeltaSigmaLCL <- exp(logDeltaSigma$Estimate - logDeltaSigma$Std..Error*1.96) 
DeltaSigma <- exp(logDeltaSigma$Estimate)

# Combine again
# All_Est <- bind_rows(All_Ests_std, All_Ests_ar, All_Ests_surv) 
All_Est <- All_Ests_std
# All_Est$ar <- All_Est$Stocknumber %in% stksNum_ar
# All_Est$surv <- All_Est$Stocknumber %in% stksNum_surv
All_Est$Param <- sapply(All_Est$Param, function(x) (unlist(strsplit(x, "[_]"))[[1]]))
All_Est <- All_Est %>%left_join(Stream, by="Stocknumber")

All_Deltas <- data.frame()
All_Deltas <- All_Ests %>% filter (Param %in% c("logDelta1", "logDelta2","sigma_delta", 
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
Pred_std <- All_Ests %>% filter (Param %in% c("LogRS_Pred_std"))
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
PredlnSMSY <- All_Ests %>% filter (Param %in% c("PredlnSMSY_S", "PredlnSMSY_O", "PredlnSMSY_CI", 
                                                "PredlnSMSYs_CI", "PredlnSMSYo_CI"))
PredlnSREP <- data.frame() 
PredlnSREP <- All_Ests %>% filter (Param %in% c("PredlnSREP_S", "PredlnSREP_O", "PredlnSREP_CI", 
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

# Tor: testing removal of stksNum_ar and stksNum_surv from original call
  # DONE
if (plot==TRUE){
  png(paste("DataOut/SR_", mod, ".png", sep=""), width=7, height=7, units="in", res=500)
  PlotSRCurve(SRDat=SRDat, All_Est=All_Est, r2=r2, removeSkagit=removeSkagit, mod=mod)
  dev.off()
  png(paste("DataOut/SRLin_", mod, ".png", sep=""), width=7, height=7, units="in", res=1000)
  PlotSRLinear(SRDat=SRDat, All_Est=All_Est, r2=r2, removeSkagit=removeSkagit)
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
  #if (mod=="Liermann_HalfNormRicVar_FixedDelta") title_plot <- "Half-normal Ricker sigma and fixed WA regression sigma"
  #title_plot <- "Separate life-histories: n=17\nFixed-effect yi (logDelta1), \nFixed-effect slope (Delta2)"
  plotWAregressionSMSY (All_Est, All_Deltas, SRDat, Stream, WA, PredlnSMSY, 
                        PredlnWA = data$PredlnWA, title1=title_plot, mod)
  dev.off()
  
  png(paste("DataOut/WAregSREP_", mod, "_wBC.png", sep=""), width=7, height=7, units="in", res=500)
  #png(paste("DataOut/WAreg_Liermann_SepRicA_UniformSigmaAPrior.png", sep=""), width=7, height=7, units="in", res=500)
  par(mfrow=c(1,1), mar=c(4, 4, 4, 2) + 0.1)
  title_plot <- "Prior Ricker sigmas and prior on WA regression sigma"
  #if (mod=="Liermann_HalfNormRicVar_FixedDelta") title_plot <- "Half-normal Ricker sigma and fixed WA regression sigma"
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
PredlnSMSY_PI <- All_Ests %>% filter (Param %in% c("PredlnSMSY", "lnSMSY"))
PredlnSREP_PI <- data.frame()
PredlnSREP_PI <- All_Ests %>% filter (Param %in% c("PredlnSREP", "lnSREP"))

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

TestSMSY <- All_Ests %>% filter (Param %in% c("TestlnSMSYo")) %>% add_column(Stock = StockNames)
TestSREP <- All_Ests %>% filter (Param %in% c("TestlnSREPo")) %>% add_column(Stock = StockNames)
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

