# Code to estimate WCVI LRPs
# Libraries
library(tidyverse)
library(ggplot2)
library(gsl)

# Functions


bEst <- function(b, SMSY, SREP){
  # Function to estimate Ricker b paramter from SMSY and SREP
  # (1 − b * SMSY) * exp( 1 − b * SMSY) = exp(1 − loga) (Explicit formula for SMSY; Eqn 10 Scheueurell 2014)
  # Rearranging:
  # log (1 − b * SMSY) + (1 − b * SMSY) = 1 - loga
  # log (1 − b * SMSY) - b * SMSY = - loga
  # loga <-  b*SMSY - log(1-b*SMSY)  
  
  # SREP <- loga / b # from Ricker 1975
  # loga <- SREP * b # Substitue this into equation above
  # SREP * b <- b*SMSY - log(1-b*SMSY)  
  # SREP <- (b * SMSY - log( 1 - b * SMSY) ) / b
  # Use this equation to solve for b
  
  pred <- (b * SMSY - log( 1 - b * SMSY) ) / b
  epsilon <- log(SREP) - log(pred)
  nLogLike <- - sum( dnorm ( epsilon, 0, 1, log = T))
  return(nLogLike)
}
bSolver <- function (SMSY, SREP) {
  # Function to estimate b from SMSY and SREP
  fit <- optimize(f = bEst, interval = c(1/SREP, 3/SREP), SMSY=SMSY, SREP=SREP)
  return(fit$minimum)
}

BMconverter.fn <- function ( SMSY, SREP, half.a = FALSE, explicit = TRUE ) {
  # Function to convert SMSY and SREP from watershed-area model into Sgen
  a.par <- exp( (0.5 - SMSY / SREP ) / 0.07 )                
  if (half.a) a.par <- (1/2) * exp( ( 0.5 - SMSY / SREP) / 0.07 )  ## Half alpha calculated here
  b.par <- log (a.par) / SREP  ## b parameter adjusted based on half alpha, assuming SREP is unchanged
  # Alternaively, could assume SMAX (= 1/b) remains constant, but SREP changes
  # b.par <- exp( (0.5 - SMSY / SREP) / 0.07) / SREP
  
  if (explicit){
    b.par <- bSolver(SMSY, SREP)
    a.par <- exp( b.par * SREP )
  }
  
  
  sgen.out <- sGenSolver( log(a.par), b.par )
  return( sgen.out )
}
SMSY <- 345
SREP <- 971

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


stockRPs <- read.csv("DataOut/WCVI_SMSY.csv")

stockRPs %>% filter(Stock != "Cypre") %>% group_by(CU)