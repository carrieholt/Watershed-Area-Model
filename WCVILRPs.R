# Code to estimate WCVI LRPs
# Libraries
library(tidyverse)
library(ggplot2)
library(gsl)

# Functions
sGenOptimum <- function ( S, theta ) {
  # Function called from sGenSolver 
  loga <- theta[1]
  b <- theta[2]
  prt <- S * exp( loga - b * S)
  sMSY <- ( 1 - gsl::lambert_W0 (exp ( 1 - loga) ) ) / b
  epsilon <- log(sMSY) - log(prt)
  nLogLike <- - sum( dnorm ( epsilon, 0, 1, log = T))
  return( nLogLike )
}


sGenSolver <- function (loga, b) {
  # Function to estimate Sgen from loga and b Ricker parameters
  theta <- c(loga, b)
  sMSY <- (1 - gsl::lambert_W0(exp(1 - loga))) / b
  fit <- optimize(f = sGenOptimum, interval = c(0, sMSY),
                  theta = theta)
  return(fit$minimum)
}


bEst <- function(b, SMSY, SREP){
  # Function to estimate Ricker b paramter from SMSY and SREP
  # (1 − b * SMSY) * exp( 1 − b * SMSY) = exp(1 − loga) (Explicit formula for SMSY; Eqn 10 Scheueurell 2014)
  # Rearranging:
  # log (1 − b * SMSY) + (1 − b * SMSY) = 1 - loga
  # log (1 − b * SMSY) - b * SMSY = - loga
  # loga <-  b*SMSY - log(1-b*SMSY)  
  
  # SREP <- loga / b # from Ricker 1975
  # loga <- SREP * b # Rearrange and substitue this into equation above
  # SREP * b <- b*SMSY - log(1-b*SMSY)  
  # SREP <- (b * SMSY - log( 1 - b * SMSY) ) / b
  # Use this equation to solve for b
  
  pred <- (b * SMSY - log( 1 - b * SMSY) ) / b
  epsilon <- log(SREP) - log(pred)
  nLogLike <- - sum( dnorm ( epsilon, 0, 1, log = T))
  return(nLogLike)
}
bSolver <- function (SMSY, SREP) {
  # Function to estimate b from SMSY and SREP, where SMAX is bounded between 1/3 of SREP and SREP
  fit <- optimize(f = bEst, interval = c(1/SREP, 3/SREP), SMSY=SMSY, SREP=SREP)
  return(fit$minimum)
}

Sgen.fn <- function ( SMSY, SREP, half.a = FALSE, const.SMAX =FALSE, explicit = TRUE , plot=FALSE) {
  # Function to convert SMSY and SREP from watershed-area model into Sgen
  # half.a = should we assume productivity has dropped by half?
  # explicit = should we use the explicit relationship between SMSY and Ricker parameters as in Scheuerell 2014?

  
  if (explicit){
    a.par <- exp( bSolver(SMSY, SREP) * SREP )
    b.par <- log( a.par ) / SREP
    if (half.a) a.par <- 0.5 * exp( bSolver(SMSY, SREP) * SREP )
    if (half.a) b.par <- log( a.par ) / SREP # b parameter adjusted based on half alpha, assuming SREP is unchanged
    # Alternaively, could assume SMAX (= 1/b) remains constant, but SREP changes
    if (half.a) if(const.SMAX) b.par <- bSolver(SMSY, SREP) #(same b as for constant prod case)
    if (half.a) if(const.SMAX) SREP <- log(a.par) / b.par
    if (half.a) SMSY <- (1 - gsl::lambert_W0(exp(1 - log(a.par) ))) / b.par
  }
  
  if( !explicit ){
    a.par <- exp( (0.5 - SMSY / SREP ) / 0.07 )                
    if (half.a) a.par <- (1/2) * exp( ( 0.5 - SMSY / SREP) / 0.07 )  ## Half alpha calculated here
    b.par <- log (a.par) / SREP  ## b parameter adjusted based on half alpha, assuming SREP is unchanged
    # Alternaively, could assume SMAX (= 1/b) remains constant, but SREP changes
    if (half.a) if(const.SMAX) b.par <- exp( (0.5 - SMSY / SREP) / 0.07) / SREP #(same b as for constant prod case)
    if (half.a) if(const.SMAX) SREP <- log(a.par) / b.par
    if (half.a) sMSY <-(1 - gsl::lambert_W0(exp(1 - log(a.par) ))) / b.par
  }

  if(plot){
    Rpred <- NA
    for (i in 1:1000){ Rpred[i]<- a.par * i * exp (- b.par * i)}
    plot(1:1000, Rpred, type="l", ylim = c (0, 1400), xlab = "Spawners",  ylab = "Recruits" )
    abline(a=0, b=1)
    abline(v=sgen.out, lty="dotted")
    abline(v=SMSY, lty="dashed")
    abline(v=(1/b.par), lty="dotdash")
    abline(v=SREP)
  }
  
  sgen.out <- sGenSolver( log(a.par), b.par )
  return( list( SGEN = sgen.out , SMSY = SMSY, SREP = SREP, apar = a.par, bpar = b.par) )
  
}

#SMSY <- 345
#SREP <- 971
Sgen.fn(SMSY, SREP, half.a = FALSE, const.SMAX = FALSE, plot=TRUE)
Sgen.fn(SMSY, SREP)
#Sgen.fn(SMSY, SREP, half.a = TRUE, const.SMAX = FALSE, plot=TRUE)


wcviRPs_long <- read.csv("DataOut/WCVI_SMSY.csv")

stock_SMSY <- wcviRPs_long %>% filter(Stock != "Cypre") %>% filter (Param == "SMSY") %>% rename(SMSY=Estimate, SMSYLL=LL, SMSYUL=UL) %>% select (-Param, -X, -CU)
stock_SREP <- wcviRPs_long %>% filter(Stock != "Cypre") %>% filter (Param == "SREP") %>% rename(SREP=Estimate, SREPLL=LL, SREPUL=UL) %>% select (-Param, -X)
wcviRPs <- stock_SMSY %>% left_join(stock_SREP, by="Stock")


# Not sure how to add a function to mutatte. I tried vecorizting Sgen.fn, but it doesn't work
# Sgen.fn_v <- Vectorize(Sgen.fn)
# wcviRPs %>% mutate( SGEN = Sgen.fn_v, SMSY, SREP )
# apply doesn work:
#SgenList <- apply (X = wcviRPs[,c('SMSY','SREP')], MARGIN = 1, FUN = Sgen.fn, SMSY=SMSY, SREP=SREP)
# Could try purrr package, and map function

# This works, but is clunky:
# Sgen <- unlist(t(mapply(FUN=Sgen.fn, wcviRPs$SMSY, wcviRPs$SREP))[,1])
# SMSY <- unlist(t(mapply(FUN=Sgen.fn, wcviRPs$SMSY, wcviRPs$SREP))[,2])
#SREP <- unlist(t(mapply(FUN=Sgen.fn, wcviRPs$SMSY, wcviRPs$SREP))[,3])
# wcviRPs <- wcviRPs %>% mutate (SGEN=Sgen) %>% mutate(SGEN=round(SGEN,0))

# This is better (using PURRR)
SGENcalcs <- map2_dfr (wcviRPs$SMSY,wcviRPs$SREP, Sgen.fn)
wcviRPs <- wcviRPs %>% mutate (SGEN = SGENcalcs$SGEN) %>% mutate(SGEN=round(SGEN,0))
wcviRPs <- wcviRPs %>% mutate (SMSYrev = SGENcalcs$SMSY) %>% mutate(SMSYrev=round(SMSYrev,0))
wcviRPs <- wcviRPs %>% mutate (SREPrev = SGENcalcs$SREP) %>% mutate(SREPrev=round(SREPrev,0))
wcviRPs <- wcviRPs %>% mutate (a.par = SGENcalcs$apar) %>% mutate(a.par=round(a.par,2))
wcviRPs <- wcviRPs %>% mutate (SMAXrev = 1/SGENcalcs$bpar) %>% mutate(SMAXrev=round(SMAXrev,0))


# Add Sgen, SMSY, SREP under half.alpha scnearios for constSREP and constSMAX. Rename columns, and re-arrange. Send to WCVI staff.
# print a plot with 3 panels and legend to show impact on half alpha on an example stocks to show staff
# Code a.par and b.par (using non-explicit calcs) and estimation of Sgen in TMB to get  SE for SGEN. Add a switch for 1/2 alpha
# Check results against R code


%>% group_by(CU)
