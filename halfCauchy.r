  test <- seq(0.001,3,0.1)
  plot(x=test, y=abs(dcauchy(test,0,1)), type="l", ylab="Probability density", xlab="" , ylim=c(0,0.4))
  
  #Normal
  lines(x=test, y=abs(dnorm(test,0,1)), lty="dashed")
  #lines(x=test, y=(dnorm(test,0.44,0.1)), lty="dashed")
  
  
  # student T, reduced to cauchy with scale=1. Adapted from Kasper Kristensen
  # https://www.stockassessment.org/datadisk/stockassessment/userdirs/user80/WBSS_HAWG_2018/TMB/include/distributions_R.hpp
  
  dt<- function(x, df=1, give_log=FALSE)
  {
    logres <- lgamma((df+1)/2) - 1/2*log(df*pi) -lgamma(df/2) - (df+1)/2*log(1+x*x/df)
    #logres <- log(1/(pi*(1+x^2)))  #Equivalent ( see http://www.math.wm.edu/~leemis/chart/UDR/PDFs/TStandardcauchy.pdf)
    if(!give_log) return (exp(logres))
    else return (logres)
  }
  dt(test)
  
  lines(x=test, y=dt(test), col="red", lty="dotted")
  
  # Cauchy = normal/sqrt(chi^2), Gelman et al. 2006. Doesnt work
  prior.scale <- 1
  xi <- dnorm (test, 0, prior.scale)
  tau.eta <-  dgamma (test,0.5,0.5) # chi^2 with 1 d.f.
  sigma.theta <- abs(xi)/sqrt(tau.eta)
  plot(x=test, y=sigma.theta, type="l")
  
  
  