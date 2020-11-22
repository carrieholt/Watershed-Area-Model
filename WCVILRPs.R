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

SGENcalcsv2 <- map2_dfr (wcviRPs$SMSY/Scale,wcviRPs$SREP/Scale, Sgen.fn, half.a = TRUE, const.SMAX = FALSE)
wcviRPs <- wcviRPs %>% mutate (SGENha.cSREP = SGENcalcsv2$SGEN) %>% mutate( SGENha.cSREP = round( SGENha.cSREP*Scale, 0 ) )
wcviRPs <- wcviRPs %>% mutate (SMSYha.cSREP = SGENcalcsv2$SMSY) %>% mutate( SMSYha.cSREP = round( SMSYha.cSREP*Scale, 0 ) )
wcviRPs <- wcviRPs %>% mutate (SREPha.cSREP = SGENcalcsv2$SREP) %>% mutate( SREPha.cSREP = round( SREPha.cSREP*Scale, 0 ) )
#wcviRPs <- wcviRPs %>% mutate (SMAXrev = 1/SGENcalcs$bpar) %>% mutate(SMAXrev=round(SMAXrev,0))

SGENcalcsv3 <- map2_dfr (wcviRPs$SMSY/Scale, wcviRPs$SREP/Scale, Sgen.fn, half.a = TRUE, const.SMAX = TRUE)
wcviRPs <- wcviRPs %>% mutate (SGENha.cSMAX = SGENcalcsv3$SGEN) %>% mutate( SGENha.cSMAX = round( SGENha.cSMAX*Scale, 0 ) )
wcviRPs <- wcviRPs %>% mutate (SMSYha.cSMAX = SGENcalcsv3$SMSY) %>% mutate( SMSYha.cSMAX = round( SMSYha.cSMAX*Scale, 0 ) )
wcviRPs <- wcviRPs %>% mutate (SREPha.cSMAX = SGENcalcsv3$SREP) %>% mutate( SREPha.cSMAX = round( SREPha.cSMAX*Scale, 0 ) )

wcviRPs 
#write.csv(wcviRPs, "DataOut/wcviRPs.csv")
#read.csv("DataOut/wcviRPs.csv")


# TMB:. add a switch for 1/2 alpha
# Check results against R code

# Do I even need TMB code for this given the I need to run this over bootstraps PRIOR to logististic regression. Or can I input all boostrapped, and estimate internally...

# TMB code to esimate Sgen

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
