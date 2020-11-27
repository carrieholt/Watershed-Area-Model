# Code to estimate WCVI LRPs
# Libraries
library(tidyverse)
library(ggplot2)
library(gsl)
library(TMB)

# Functions
source("helperFunctions.r")

# # Example: Artlish
# SMSY <- 345 
# SREP <- 971
# 
# png(paste("DataOut/Artlish_WCVI_SRcurve.png", sep=""), width=4, height=7, units="in", res=500)
#   par(mfcol=c(3,1),  mar = c(4, 4, 2.5, 2) + 0.1)
#   Sgen.fn(SMSY, SREP, half.a = FALSE, const.SMAX = FALSE, plot=TRUE)
#   Sgen.fn(SMSY, SREP, half.a = TRUE, const.SMAX = FALSE, plot=TRUE)
#   Sgen.fn(SMSY, SREP, half.a = TRUE, const.SMAX = TRUE, plot=TRUE)
# dev.off()

wcviRPs_long <- read.csv("DataOut/WCVI_SMSY.csv")

# Remove Cypre as it's not a core indicator (Diana McHugh, 22 Oct 2020)
stock_SMSY <- wcviRPs_long %>% filter(Stock != "Cypre") %>% filter (Param == "SMSY") %>% rename(SMSY=Estimate, SMSYLL=LL, SMSYUL=UL) %>% select (-Param, -X)#, -CU)
stock_SREP <- wcviRPs_long %>% filter(Stock != "Cypre") %>% filter (Param == "SREP") %>% rename(SREP=Estimate, SREPLL=LL, SREPUL=UL) %>% select (-Param, -X)
wcviRPs <- stock_SMSY %>% left_join(stock_SREP, by="Stock")

# Calculate scale for each stock
digits <- count.dig(stock_SMSY$SMSY)
Scale <- 10^(digits)



# Not sure how to add a function to mutate. I tried vecorizting Sgen.fn, but it doesn't work
# Sgen.fn_v <- Vectorize(Sgen.fn)
# wcviRPs %>% mutate( SGEN = Sgen.fn_v, SMSY, SREP )
# apply doesn work:
# SgenList <- apply (X = wcviRPs[,c('SMSY','SREP')], MARGIN = 1, FUN = Sgen.fn, SMSY=SMSY, SREP=SREP)
# This works, but is clunky:
# Sgen <- unlist(t(mapply(FUN=Sgen.fn, wcviRPs$SMSY, wcviRPs$SREP))[,1])
# SMSY <- unlist(t(mapply(FUN=Sgen.fn, wcviRPs$SMSY, wcviRPs$SREP))[,2])
# SREP <- unlist(t(mapply(FUN=Sgen.fn, wcviRPs$SMSY, wcviRPs$SREP))[,3])
# wcviRPs <- wcviRPs %>% mutate (SGEN=Sgen) %>% mutate(SGEN=round(SGEN,0))

# This is better (using PURRR)
SGENcalcs <- map2_dfr (wcviRPs$SMSY/Scale,wcviRPs$SREP/Scale, Sgen.fn)
wcviRPs <- wcviRPs %>% mutate (SGEN = SGENcalcs$SGEN) %>% mutate(SGEN=round(SGEN*Scale,0))
wcviRPs <- wcviRPs %>% mutate (a.par = SGENcalcs$apar) %>% mutate(a.par=round(a.par,2))

wcviRPs <- wcviRPs[c("Stock", "SGEN", "SMSY", "SMSYLL", "SMSYUL", "SREP", "SREPLL", "SREPUL", "a.par")]#"CU"

# SGENcalcsv2 <- map2_dfr (wcviRPs$SMSY/Scale,wcviRPs$SREP/Scale, Sgen.fn, half.a = TRUE, const.SMAX = FALSE)
# wcviRPs <- wcviRPs %>% mutate (SGENha.cSREP = SGENcalcsv2$SGEN) %>% mutate( SGENha.cSREP = round( SGENha.cSREP*Scale, 0 ) )
# wcviRPs <- wcviRPs %>% mutate (SMSYha.cSREP = SGENcalcsv2$SMSY) %>% mutate( SMSYha.cSREP = round( SMSYha.cSREP*Scale, 0 ) )
# wcviRPs <- wcviRPs %>% mutate (SREPha.cSREP = SGENcalcsv2$SREP) %>% mutate( SREPha.cSREP = round( SREPha.cSREP*Scale, 0 ) )
# ###wcviRPs <- wcviRPs %>% mutate (SMAXrev = 1/SGENcalcs$bpar) %>% mutate(SMAXrev=round(SMAXrev,0))
# 
# SGENcalcsv3 <- map2_dfr (wcviRPs$SMSY/Scale, wcviRPs$SREP/Scale, Sgen.fn, half.a = TRUE, const.SMAX = TRUE)
# wcviRPs <- wcviRPs %>% mutate (SGENha.cSMAX = SGENcalcsv3$SGEN) %>% mutate( SGENha.cSMAX = round( SGENha.cSMAX*Scale, 0 ) )
# wcviRPs <- wcviRPs %>% mutate (SMSYha.cSMAX = SGENcalcsv3$SMSY) %>% mutate( SMSYha.cSMAX = round( SMSYha.cSMAX*Scale, 0 ) )
# wcviRPs <- wcviRPs %>% mutate (SREPha.cSMAX = SGENcalcsv3$SREP) %>% mutate( SREPha.cSMAX = round( SREPha.cSMAX*Scale, 0 ) )

wcviRPs 
#write.csv(wcviRPs, "DataOut/wcviRPs.csv")
#read.csv("DataOut/wcviRPs.csv")


#----------------------------------------------------------------------------------
# Sum escapements across indicators within inlets
# ---------------------------------------------------------------------------------

WCVIEsc <- data.frame(read.csv("DataIn/WCVIEsc.csv", row.names="Yr")) %>% select (-"Little.Zeballos")

# Take "." out of nameas in escapement data
WCVIEsc_names <- sapply(colnames(WCVIEsc), function(x) (gsub(".", " ", x, fixed=TRUE) ) )
WCVIEsc_names <- sapply(WCVIEsc_names, function(x) (gsub("Bedwell Ursus", "Bedwell/Ursus", x, fixed=TRUE) ) )
WCVIEsc_names <- sapply(WCVIEsc_names, function(x) (gsub("Nootka Esperanza", "Nootka/Esperanza", x, fixed=TRUE) ) )
colnames(WCVIEsc) <- WCVIEsc_names 

EnhStocks <- c("Artlish", "Burman",  "Conuma", "Leiner", "Nitinat", "Sarita",  "Somass",  "Zeballos", "San Juan")

remove.EnhStocks <- FALSE#TRUE

if (remove.EnhStocks) {WCVIEsc <- WCVIEsc %>% select(-EnhStocks) }

Years <- rownames(WCVIEsc)

# Remove Cypre as it's not an indicator stocks
WCVIStocks <- read.csv("DataIn/WCVIStocks.csv") %>% filter (Stock != "Cypre")

if (remove.EnhStocks) WCVIStocks <- WCVIStocks %>% filter (Stock %not in% EnhStocks)
Inlet_Names <- unique(WCVIStocks$Inlet)
Inlet_Sum <- matrix(NA, nrow=length(Years), ncol=length(Inlet_Names))
colnames(Inlet_Sum) <- Inlet_Names


# Sum escapements across stocks within inlets
for (i in 1:length(Inlet_Names)) {
  # For each inlet, which are the component indicators
  Ins <- WCVIStocks %>% filter(Inlet==Inlet_Names[i]) %>% pull(Stock)
  WCVIEsc_Inlets <- matrix(NA, nrow= length(Years), ncol= length(Ins))
  
  #  Make a matrix of escapements of component indicators
  for (j in 1:length(Ins)){
    WCVIEsc_Inlets[,j] <- WCVIEsc %>% select(as.character(Ins[j])) %>% pull()
    #xx<- wcviRPs %>% filter(StockNames == as.character(Ins[j])) %>% pull(SGEN)
    
  }
  
  # Sum the escapement of component indicators, removing years where there are any NAs
  Inlet_Sum[,i] <- apply(WCVIEsc_Inlets, 1, sum, na.rm=F)
}
#WCVIEsc <- cbind(WCVIEsc, Inlet_Sum) 


# Sum escapements across stocks within CUs
CU_Sum <- matrix(NA, nrow=length(Years), ncol=length(CU_Names))
colnames(CU_Sum) <- CU_Names
for (k in 1:length(CU_Names)) {
  # For each CU, which are the component indicators
  CUss <- WCVIStocks %>% filter(CU==CU_Names[k]) %>% pull(Stock)
  WCVIEsc_CUs <- matrix(NA, nrow= length(Years), ncol= length(CUss))
  
  #  Make a matrix of escapements of component indicators
  for (j in 1:length(CUss)){
    WCVIEsc_CUs[,j] <- WCVIEsc %>% select(as.character(CUss[j])) %>% pull()
    #xx<- wcviRPs %>% filter(StockNames == as.character(Ins[j])) %>% pull(SGEN)
    
  }
  
  # Sum the escapement of component indicators, removing years where there are any NAs
  CU_Sum[,k] <- apply(WCVIEsc_CUs, 1, sum, na.rm=F)
}
WCVIEsc <- cbind(WCVIEsc, Inlet_Sum, CU_Sum) 

# Assess status for each inlet relative to inlet-level SGEN for each year

Inlet_Status <- matrix(NA, nrow=length(Years), ncol=length(Inlet_Names) )
colnames(Inlet_Status) <- Inlet_Names

for (i in 1:length(Inlet_Names)) {
  Inlet_Status[,i] <- ( Inlet_Sum[,i] > (wcviRPs %>% filter(Stock == as.character(Inlet_Names[i])) %>% pull(SGEN)) )
}

Inlet_Status <- as.data.frame(Inlet_Status)

# Assess status for each CU for each year of the time-series
CU_Names <- unique(WCVIStocks$CU)
CU_Status <- matrix(NA, nrow=length(Years), ncol=length(CU_Names))
colnames(CU_Status) <- CU_Names

stock.LRP <- TRUE
if(stock.LRP){
  for (k in 1:length(CU_Names)) {
    # For each CU, which are the component indicators
    CU_ins <- unique( WCVIStocks %>% filter(CU==CU_Names[k]) %>% pull(Inlet) )
    
    isAbove <- matrix(NA, nrow= length(Years), ncol= length(CU_ins))
    
    #  Make a matrix of status of component inlets. Is each inlet > Sgen values?
    for (i in 1:length(CU_ins)){
      isAbove[,i] <- Inlet_Status %>% select(as.character(CU_ins[i])) %>% pull()
    }
    
    # CU-level status: are ALL inlets above their Sgen values
    isAboveFun <- function(x){ floor(sum( as.numeric(x), na.rm=F) / length(x) ) }
    CU_Status[,k] <- apply(X= isAbove, MARGIN = 1, FUN=isAboveFun)
  }
  CU_Status <- as.data.frame(CU_Status)
  
}

CU.LRP <- FALSE
if(CU.LRP){
  for (k in 1:length(CU_Names)){
    CU_Status[,k] <- ( CU_Sum[,k] > (wcviRPs %>% filter(Stock == as.character(CU_Names[k])) %>% pull(SGEN)) )
  }
  CU_Status <- as.data.frame(CU_Status)
  
}

# Proportion of CUs that are not in the red zone

ppnAboveFun <- function(x) {sum( as.numeric(x), na.rm=F) / length(x) }
SMU_ppn <- apply(X=CU_Status, MARGIN=1, FUN=ppnAboveFun)

#----------------------------------------------------------------------------------
# Logistic regression
#----------------------------------------------------------------------------------
SMU_Esc <- apply(Inlet_Sum, 1, sum, na.rm=F)

SMUlogisticData <- data.frame(SMU_Esc) %>% add_column(ppn=SMU_ppn, Years=as.numeric(Years)) %>% filter(SMU_Esc != "NA")

data <- list()
data$N_Stks <- length(CU_Names)
Scale <- # Calculate scale for each stock
digits <- count.dig(SMU_Esc)
ScaleSMU <- min(10^(digits -1 ), na.rm=T)

data$LM_Agg_Abund <- SMUlogisticData$SMU_Esc/ScaleSMU
data$N_Above_BM <- SMUlogisticData$ppn * data$N_Stks
data$Pred_Abund <- seq(0, max(data$LM_Agg_Abund)*1.1, 0.1)
data$p <- 0.95#0.67
data$Penalty <- as.numeric(TRUE)
# consider adding prior stuch that when p is near zero, i.e., 0.01, the aggegate abundance is lognorm with CIs at size of smallest CU - size of all CUs combined.
# or set B0 such  that when p=0.01, then aggregate abundance = smallest CU, as model wants to go lower than this.

# To find pars for log-normal penalty. First look at distribution in log-space:
# Mean in log space = ave of min and max values. Sd = that which allows the density = 0.05 at min and max values
min <- min(apply(CU_Sum, 2, mean, na.rm=T), na.rm=T)#min(CU_Sum, na.rm=T)
max <- mean(apply(CU_Sum, 1, sum, na.rm=F), na.rm=T)
# Plot in raw space
plot(x=seq(min,max,100), y=dnorm(seq(min,max,100), mean=mean(c(min,max)), sd=22000),type="l")
#sum(dnorm(seq(min,max,1), mean=mean(c(min,max)), sd=22000))# Should give 95% density

# # Alternative priors
# #Plot in log-space. When I do this, I find I need sd=1.2 to get the density at 0.05 at low end, 7.04
# plot(x=seq(log(min),log(max),0.01), y=dnorm(seq(log(min),log(max),0.01), mean=(log(min) + log(max))/2, sd=1.18),type="l")
# 
# # Plot in raw space
# plot(x=seq(min,max,100), y=dlnorm(seq(min,max,100), meanlog=(log(min)+log(max))/2, sdlog=1.18),type="l")


#mean of distribution in raw space: exp((log(min)+log(max))/2+1.18^2/2)
#SD in raw space:  sqrt( (exp(1.18^2)-1) * exp( 2* ( (log(min)+log(max))/2 ) +1.18^2) )
data$B_penalty_mu <- mean(c(min,max))/ScaleSMU#min/ScaleSMU#1#mean(c(log(min),log(max)))
data$B_penalty_sig <- 22000/ScaleSMU#(min/ScaleSMU)*10#1.18


param <- list()
param$B_0 <- -2
param$B_1 <- 0.1
  

#dyn.unload(dynlib(paste("TMB_Files/Logistic_LRPs", sep="")))
#compile(paste("TMB_Files/Logistic_LRPs.cpp", sep=""))
dyn.load(dynlib(paste("TMB_Files/Logistic_LRPs", sep="")))
obj <- MakeADFun(data, param, DLL="Logistic_LRPs", silent=TRUE)#random = c( "logDelta2"), 

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) 
summary(sdreport(obj), p.value=TRUE)

All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
All_Ests$Param <- sapply(All_Ests$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))
Preds <- All_Ests %>% filter(Param == "Logit_Preds")
All_Ests <- All_Ests %>% filter(!(Param %in% c( "Logit_Preds"))) 

out <- list()
out$All_Ests <- All_Ests


Logistic_Data <- data.frame(yr = SMUlogisticData$Years, 
                            yy = SMUlogisticData$ppn, xx = SMUlogisticData$SMU_Esc)

out$Logistic_Data <- Logistic_Data

Logistic_Fits <- data.frame(xx = data$Pred_Abund*ScaleSMU, fit = inv_logit(Preds$Estimate),
                            lwr = inv_logit(Preds$Estimate - 1.96*Preds$Std..Error),
                            upr = inv_logit(Preds$Estimate + 1.96*Preds$Std..Error))

out$Preds <- Logistic_Fits

out$LRP <- data.frame(fit = (All_Ests %>% filter(Param == "Agg_LRP") %>% pull(Estimate))*ScaleSMU, 
                      lwr = (All_Ests %>% filter(Param == "Agg_LRP") %>% mutate(xx =Estimate - 1.96*Std..Error) %>% 
                               pull(xx) ) * ScaleSMU,
                      upr = (All_Ests %>% filter(Param == "Agg_LRP") %>% mutate(xx =Estimate + 1.96*Std..Error) %>% 
                               pull(xx) ) * ScaleSMU)

out

  
#----------------------------------------------------------------------------------------------
# R version of logistic regression
#----------------------------------------------------------------------------------------------

ModDat <- data.frame(xx=data$LM_Agg_Abund, yy=SMUlogisticData$ppn)
Fit_Mod <- glm( yy ~ xx , family = quasibinomial, data=ModDat)#or family=binomial, which gives much larger SEs, and assumes var=1.
summary(Fit_Mod)$coefficients
LRP <- (log(data$p/(1-data$p)) - Fit_Mod$coefficients[[1]])/ Fit_Mod$coefficients[[2]]
# use MASS function to get "dose" 
library(MASS)
Dose <- dose.p(Fit_Mod, p=data$p)
Dose

#  - Make x vector to predict with this model, for plotting
xx <- data.frame(xx = seq(0, max(data$LM_Agg_Abund*1.25), by=(max(data$LM_Agg_Abund*1.25)/1000)))

# - Create model predictions that include standard error fit
preds <- predict.glm(Fit_Mod, newdata = xx, type = "link", se.fit = TRUE)#, dispersion = )

# Create a confidence interval (lwr, upr) on the link scale as the fitted value plus or minus 1.96 times the standard error
critval <- 1.96 ## approx 95% CI
upr <- preds$fit + (critval * preds$se.fit)
lwr <- preds$fit - (critval * preds$se.fit)
fit <- preds$fit

# Transform confidene interval using the inverse of the link function to map the fitted values and the upper and lower limits of the interval back on to the response scale
fit2 <- Fit_Mod$family$linkinv(fit)
upr2 <- Fit_Mod$family$linkinv(upr)
lwr2 <- Fit_Mod$family$linkinv(lwr)

# Load predicted data on response scale into dataframe
preddata<-xx
preddata$fit<-fit2
preddata$lwr <- lwr2 
preddata$upr <- upr2 

### Calculate confidence intervals for LRP 
LRP_lwr <- LRP - critval*attr(Dose, "SE")[[1]]
LRP_upr <- LRP + critval*attr(Dose, "SE")[[1]]
# These seem much tighter than other method?



list.out<-list()
list.out$Logistic_Data <- ModDat
list.out$model <- Fit_Mod
list.out$Preds <- preddata
list.out$LRP<-data.frame(fit = LRP, lwr = LRP_lwr, upr = LRP_upr)


list.out


# # Equivalent code to sum escapaement across indicators within inlets, specific to WCVI for plotting
# 
# # Add inlets, excluding San Juan and Nitinat, as they contain only one stock each
# # NAs must be removed from the sums first
# BarkleyNA <- apply( cbind( WCVIEsc$Sarita, WCVIEsc$Somass, WCVIEsc$Nahmint), 1, sum, na.rm=F)
# ClayoquotNA <- apply( cbind( WCVIEsc$Bedwell.Ursus, WCVIEsc$Megin, WCVIEsc$Moyeha, WCVIEsc$Tranquil), 1, sum, na.rm=F)
# Nootka.EsperanzaNA <- apply( cbind( WCVIEsc$Burman , WCVIEsc$Conuma, WCVIEsc$Gold, WCVIEsc$Leiner, WCVIEsc$Tahsis, 
#                                     WCVIEsc$Zeballos), 1, sum, na.rm=F)
# KyuquotNA <- apply( cbind( WCVIEsc$Artlish, WCVIEsc$Kaouk, WCVIEsc$Tahsish), 1, sum, na.rm=F)
# QuatsinoNA <- apply( cbind( WCVIEsc$Cayeghle, WCVIEsc$Marble), 1, sum,  na.rm=F)
# 
# WCVIEsc <- WCVIEsc %>% mutate(Barkley = BarkleyNA, Clayoquot = ClayoquotNA, Nootka.Esperanza =Nootka.EsperanzaNA,
#                               Kyuquot = KyuquotNA, Quatsino = QuatsinoNA )
# 
# WCVIEsc

#-------------------------------------------------------------------------------------
# TMB version of code to estimate Sgen from SMSY and SREP
#-------------------------------------------------------------------------------------

# Do not need TMB code given I need to run this over bootstraps PRIOR to logististic regression. 
# I could input all boostrapped PI draws, and estimate LRP internally for each draw, however wrangling with data is difficult 
# in TMB and since there is only one estimation step: logistic regression, could simply implement in R

# To do: (1) Add a switch for 1/2 alpha
# (2) Check results against R code


SMSY <- wcviRPs %>% pull(SMSY)
SREP <- wcviRPs %>% pull(SREP)

# Calculate scale for each stock
digits <- count.dig(SMSY)
Scale <- 10^(digits)

SMSY <- SMSY/Scale
SREP <- SREP/Scale


data <- list()
data$SMSY <- SMSY
data$SREP <- SREP
data$Inlets <- read.csv("DataIn/WCVIStocks.csv") %>% filter (Stock != "Cypre") %>% pull(SA_ind)
data$N_inlets <- length(unique(read.csv("DataIn/WCVIStocks.csv") %>% filter (Stock != "Cypre") %>% pull(SA_ind)))

#data$Scale <- Scale


param <- list()
param$RicB <- 1/(data$SREP/2) #initialize SMAX at half SREP
param$logSgen <- log(data$SMSY/2)#initialize SMAX at half SMSY

# Compile model if changed:
#dyn.unload(dynlib("TMB_Files/WA_Sgen"))
#compile("TMB_Files/WA_Sgen.cpp")
dyn.load(dynlib("TMB_Files/WA_Sgen"))
obj <- MakeADFun(data, param, DLL="WA_Sgen", silent=TRUE)

# b is bounded between 1/3 of SREP and SREP
lower <- 1/data$SREP
upper <- 3/data$SREP
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))#, lower=rep(0,20), upper=log(SMSY))
pl <- obj$env$parList(opt$par) 
#summary(sdreport(obj), p.value=TRUE)


# THe summation of Sgens across inlets is not working because Sgen's are each scaled differntly
# Actually, do this in R as organizing data is easier there (and just almost as fast to run)

#exp(pl$logSgen)*Scale
#1/((1/pl$RicB)*Scale)
