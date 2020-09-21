# Code to plot prior distributions for sigma

library(invgamma)
plot_inv_gamma_only <- FALSE

# Functions
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}

# Code
test <- seq(0.00001,9,len=10000)
  
plot(x=test, y=abs(dcauchy(test,0,1)), type="n", ylab="Probability density", xlab="Ricker Sigma or LogA sigma" , ylim=c(0,0.8), xlim=c(0,2.5))
cols<-viridis(4, alpha=0.9)
    
#Inverse gamma on sqrt(variance)=sigma
shape<-rate<-0.001
lines(sqrt(test), sqrt(dinvgamma(test,shape=shape, rate=rate)), lwd=1, col=t_col(color=cols[1], percent=90))
shape<-rate<-0.01
lines(sqrt(test), sqrt(dinvgamma(test,shape=shape, rate=rate)), lwd=1, col=t_col(color=cols[1], percent=70))
shape<-rate<-0.1
lines(sqrt(test), sqrt(dinvgamma(test,shape=shape, rate=rate)), lwd=1, col=t_col(color=cols[1], percent=50))
shape<-rate<-1
lines(sqrt(test), sqrt(dinvgamma(test,shape=shape, rate=rate)), lwd=1, col=t_col(color=cols[1], percent=30))

if (!plot_inv_gamma_only){
  #Half-Normal on sigma
  lines(x=test, y=abs(dnorm(test,0,1)), col=cols[2], lwd=2)
  #lines(x=test, y=(dnorm(test,0.44,0.1)), lty="dashed")
  
  # Half-cauchy on sigma, as implemented in TMB
  # student T, reduced to cauchy with scale=1. Adapted from Kasper Kristensen
  # https://www.stockassessment.org/datadisk/stockassessment/userdirs/user80/WBSS_HAWG_2018/TMB/include/distributions_R.hpp
  dt<- function(x, df=1, give_log=FALSE)
  {
    logres <- lgamma((df+1)/2) - 1/2*log(df*pi) -lgamma(df/2) - (df+1)/2*log(1+x*x/df)
    #logres <- log(1/(pi*(1+x^2)))  #Equivalent ( see http://www.math.wm.edu/~leemis/chart/UDR/PDFs/TStandardcauchy.pdf)
    if(!give_log) return (exp(logres))
    else return (logres)  
  }
  
  lines(x=test, y=dt(test), col=cols[3], lwd=2)
  
  #Altrenative, equivalent half cauchy using R function
  #lines(x=test, y=abs(dcauchy(test,0,1)), lwd=2, col="red")
  
  #Uniform 0-2
  lines(x=c(0,2,2,3), y=c(0.5,0.5,0,0), col=cols[4], lwd=2)
  
  legend(x=1.1, y=0.75, legend=c("Inverse gamma(1,1)", "Inverse gamma(0.1,0.1)", "Inverse gamma(0.01,0.01)", 
                                 "Inverse gamma(0.001,0.001)", "Half Normal (0,1)", "Half Cauchy (0,1)", "Uniform (0,2)"),
         col=c(t_col(cols[1], 30), t_col(cols[1], 50), t_col(cols[1], 70), t_col(cols[1], 90), cols[2:4]), bty="n", lwd=2)
  
}

# Cauchy = normal/sqrt(chi^2), Gelman et al. 2006. Doesnt work
prior.scale <- 1
xi <- dnorm (test, 0, prior.scale)
tau.eta <-  dgamma (test,0.5,0.5) # chi^2 with 1 d.f.
sigma.theta <- abs(xi)/sqrt(tau.eta)
#plot(x=test, y=sigma.theta, type="l")

  
if (plot_inv_gamma_only){
  legend(x=1.1, y=0.75, legend=c("Inverse gamma(1,1)", "Inverse gamma(0.1,0.1)", "Inverse gamma(0.01,0.01)", 
                                 "Inverse gamma(0.001,0.001)"),
         col=c(t_col(cols[1], 30), t_col(cols[1], 50), t_col(cols[1], 70), t_col(cols[1], 90)), bty="n", lwd=2)
  
}
  

# Add histogram of Ricker sigmas
  