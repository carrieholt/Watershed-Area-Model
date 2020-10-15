# Code to plot prior distributions for sigma

library(invgamma)
library(viridis)

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

plotPriors <- function (plot_inv_gamma_only, Delta){
  #plot_inv_gamma_only TRUE = then only plot inverse gamma distributions and not cauchy, normal or uniform
  #Delta TRUE = plot priors for the SD of the Watershed-area model
  #Delta FALSE = plot priors for the SD of the Ricker model 
  
  test <- seq(0.00001,9,len=10000)
  
  #plot(x=test, y=abs(dcauchy(test,0,1)), type="n", ylab="Probability density", xlab="Ricker Sigma (or LogA sigma)" , ylim=c(0,0.8), xlim=c(0,2.5))
  plot(x=test, y=abs(dcauchy(test,0,1)), type="n", ylab="Probability density", xlab="Sigma for watershed-area regression" , ylim=c(0,0.8), xlim=c(0,2.5))
  cols<-viridis(4, alpha=0.9)
  
  
  if(!Delta){
    # Add histogram of Ricker sigmas
    
    # Ricker sigma from Parken data, excluding Upper Columbia and Siluetz, where they included AR(1) term in model
    # RicSig <- read.csv("DataIn/ParkenRicSig.csv") 
    
    # Ricker sigma from PSE Chinook SR data extracted 15 Oct 2020. See sigR_metaanalysis.R
    RicSig <- exp(read.table("DataOut/PSE_sigma.csv")$Estimate) 
    hist_out <- hist(RicSig, plot=FALSE)
    lowlimit <- hist_out$breaks[1] # Indicates where bars should start on the plot, assuming units of 0.1
    barplot(height =  c( rep(0,lowlimit*10),(hist_out$density)*0.2), width=0.1, col=grey(0.8), border=grey(0.7), space=0, add=TRUE)
    RicSig_Thorson <- 0.69 #Thorson et al. 2014 marginal sigma from predictive distribution from hierarchical Ricker model of 20 salmonid stocks from Myers et al 1995
    RicSigSD_Thorson <- 0.27
    abline(v=RicSig_Thorson, lwd=2, lty="dotted", col=grey(0.8))
    #polgyon(x=c(RicSig_Thorson-RigSig_Thorson))
  }
  
  # Better to use external data to form priors, 
  # E.g., NCEAS State of Alaska Salmon and People project:	
  # Brendan Connors: Chinook (n = 75) most from US West Coast and of questionable quality, remaining ~20 stocks are from AK and of higher quality
  
  
  
  #Inverse gamma on sqrt(variance)=sigma
  if(!Delta){
    shape<-rate<-0.001
    lines(sqrt(test), sqrt(dinvgamma(test,shape=shape, rate=rate)), lwd=2, col=t_col(color=cols[1], percent=90))
  }
  
  shape<-rate<-0.01
  lines(sqrt(test), sqrt(dinvgamma(test,shape=shape, rate=rate)), lwd=4, col=t_col(color=cols[1], percent=70))
  shape<-rate<-0.1
  lines(sqrt(test), sqrt(dinvgamma(test,shape=shape, rate=rate)), lwd=4, col=t_col(color=cols[1], percent=50))
  shape<-rate<-1
  lines(sqrt(test), sqrt(dinvgamma(test,shape=shape, rate=rate)), lwd=4, col=t_col(color=cols[1], percent=30))
  
  if (!plot_inv_gamma_only){
    
    if(!Delta){
      #Half-Normal on sigma
      lines(x=test, y=abs(dnorm(test,0,1)), col=cols[2], lwd=4)
    
    
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
    
      lines(x=test, y=dt(test), col=cols[3], lwd=4)
    
      #Altrenative, equivalent half cauchy using R function
      #lines(x=test, y=abs(dcauchy(test,0,1)), lwd=2, col="red")
    
      #Uniform 0-1
      lines(x=c(0,2,2,3), y=c(0.5,0.5,0,0), col=cols[4], lwd=2)
    }
    
    if(Delta){
      #Uniform 0-1
      lower <- 0.21 #See KFrun.R
      upper <- sd(log(read.csv("DataIn/ParkenSMSY.csv")$SMSY))
      lines(x=c(0,lower,lower, upper,upper,3), y=c(0,0,0.8,0.8,0,0), col=cols[4], lwd=4)
      #Normal bounded
      norm_scalar <- 1.8
      testb <- seq(lower, upper,len=100)
      lines(x=testb, y=abs(dnorm(testb,0.8,0.28))/norm_scalar, col=cols[2], lwd=4)#See KRrun.R for N(0.8,0.28)
      lines(x=c(0,lower, lower), y=c(0, 0, dnorm(lower,0.8,0.28))/norm_scalar, col=cols[2], lwd=4)
      lines(x=c(upper,upper, 3), y=c(dnorm(upper,0.8,0.28)/norm_scalar, 0, 0), col=cols[2], lwd=4)
      
    }
    
    if(Delta){
      # For Delta sigma, vertical line SD of SMSY among 25 stocks
      abline(v= sd(log(read.csv("DataIn/ParkenSMSY.csv")$SMSY)), col=grey(0.5))
      abline(v= sqrt((0.293+0.146)/2), col=grey(0.5), lty="dashed")#See Parken et al. (2006)
      abline(v= 0.21, col=grey(0.5), lty="dotted")#See KFrun.R
    }
    
    if(!Delta){
      legend(x=1.1, y=0.78, legend=c("Inverse gamma(1,1)", "Inverse gamma(0.1,0.1)", "Inverse gamma(0.01,0.01)", 
                                     "Inverse gamma(0.001,0.001)", "Half Normal (0,1)", "Half Cauchy (0,1)", "Uniform (0,2)"),
             col=c(t_col(cols[1], 30), t_col(cols[1], 50), t_col(cols[1], 70), t_col(cols[1], 90), cols[2:4]), bty="n", lwd=2) 
    }
    if(Delta){
      # legend(x=1.4, y=0.78, legend=c("Inverse gamma(1,1)", "Inverse gamma(0.1,0.1)", "Inverse gamma(0.01,0.01)",
      #                               "Half Normal (0,1)", "Half Cauchy (0,1)", "SD of ln(SMSY) Parken et al.",
      #                               "sigma WA regression Parken et al.", "Med. SD of time-varying SMSYs"),
      #        col=c(t_col(cols[1], 30), t_col(cols[1], 50), t_col(cols[1], 70), cols[2:3], rep(grey(0.5),3)),
      #        lty=c(rep("solid", 6), "dashed", "dotted"), bty="n", lwd=2, cex=0.8)
      # legend(x=1.4, y=0.78, legend=c("Inverse gamma(1,1)", "Inverse gamma(0.1,0.1)", "Inverse gamma(0.01,0.01)", 
      #                                "Normal bounded", "Uniform bounded", "SD of ln(SMSY) Parken et al.", 
      #                                "sigma WA regression Parken et al.", "Med. SD of time-varying SMSYs"),
      #        col=c(t_col(cols[1], 30), t_col(cols[1], 50), t_col(cols[1], 70), cols[2], cols[4], rep(grey(0.5),3)),  
      #        lty=c(rep("solid", 6), "dashed", "dotted"), bty="n", lwd=2, cex=0.8) 
      
    }
  }
  
  if (plot_inv_gamma_only){
    if(!Delta){
      legend(x=1.1, y=0.78, legend=c("Inverse gamma(1,1)", "Inverse gamma(0.1,0.1)", "Inverse gamma(0.01,0.01)", 
                                     "Inverse gamma(0.001,0.001)"),
             col=c(t_col(cols[1], 30), t_col(cols[1], 50), t_col(cols[1], 70), t_col(cols[1], 90)), bty="n", lwd=2)
    }
    if(Delta){
      # For Delta sigma, vertical line SD of SMSY among 25 stocks
      abline(v= sd(log(read.csv("DataIn/ParkenSMSY.csv")$SMSY)), col=grey(0.5))
      abline(v= sqrt((0.293+0.146)/2), col=grey(0.5), lty="dashed")#See Parken et al. (2006)
      abline(v= 0.21, col=grey(0.5), lty="dotted")#See KFrun.R
      
      
      legend(x=1.4, y=0.78, legend=c("Inverse gamma(1,1)", "Inverse gamma(0.1,0.1)", "Inverse gamma(0.01,0.01)", 
             "SD of ln(SMSY) Parken et al.", "sigma WA regression Parken et al.", "Med. SD of time-varying SMSYs"),
             col=c(t_col(cols[1], 30), t_col(cols[1], 50), t_col(cols[1], 70), rep(grey(0.5),3)), 
             lty=c(rep("solid", 4), "dashed", "dotted"),
             bty="n", lwd=2, cex=0.8)
    }
    
  }
  
}

# Cauchy = normal/sqrt(chi^2), Gelman et al. 2006. Doesnt work
# prior.scale <- 1
# xi <- dnorm (test, 0, prior.scale)
# tau.eta <-  dgamma (test,0.5,0.5) # chi^2 with 1 d.f.
# sigma.theta <- abs(xi)/sqrt(tau.eta)
# plot(x=test, y=sigma.theta, type="l")


par(mfrow=c(1,1))  

png(paste("DataOut/RicPriors_InvGamma.png", sep=""), width=7, height=7, units="in", res=500)
plotPriors(plot_inv_gamma_only=TRUE, Delta=FALSE)
dev.off()

#png(paste("DataOut/DeltaPriors_InvGamma.png", sep=""), width=7, height=7, units="in", res=500)
#plotPriors(plot_inv_gamma_only=TRUE, Delta=TRUE)
#dev.off()

#png(paste("DataOut/RicPriors_sm.png", sep=""), width=7, height=7, units="in", res=500)
#plotPriors(plot_inv_gamma_only=FALSE, Delta=FALSE)
#dev.off()


#png(paste("DataOut/DeltaPriors.png", sep=""), width=7, height=7, units="in", res=500)
#plotPriors(plot_inv_gamma_only=FALSE, Delta=TRUE)
#dev.off()
