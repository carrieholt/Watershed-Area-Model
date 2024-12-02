#------------------------------------------------------------------------------
# Options for likelihood penalties IWAM
# 1 Dec 2024
#------------------------------------------------------------------------------

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


#------------------------------------------------------------------------------
# Sigma in watershed-area regression
# Option 1: Gamma penalty on precision
sigma <- seq(0.00001,9,len=10000)
cols<-viridis(4, alpha=0.9)

plot(x=sigma, y=sigma, type="n", ylab="Density", 
     xlab="Sigma for watershed-area regression" ,
     main="Gamma penalty on precision",
     ylim=c(0,0.6), xlim=c(0,2.5))
abline(v= sd(log(read.csv("DataIn/ParkenSMSY.csv")$SMSY)), col=grey(0.5), 
       lty="dashed")
abline(v= sd(log(read.csv("DataIn/ParkenSREP.csv")$SREP)), col=grey(0.5), 
       lty="dotted")
#abline(v= sqrt(0.293), col=grey(0.5), lty="dashed")#See Parken et al. (2006)
#abline(v= sqrt(0.146), col=grey(0.5), lty="dashed")#See Parken et al. (2006)

rate <- rev(c(0.5, 1, 2, 3, 4))
per <- seq(90, 10, length.out=length(rate))
mean <- 3
for (i in 1:length(rate)){
  lines(sigma, dgamma(1/(sigma^2),shape=rate[i]*mean, scale=1/rate[i]), lwd=4, 
        col=t_col(color=cols[1], percent= per[i]))
}
legend(x=0, y=0.6, legend=c("sigma ln(SMSY) Parken et al. (upper bound)", 
                             "sigma ln(SREP) Parken et al. (upper bound)"),
       col=grey(0.5),
       lty=c("dashed", "dotted"), bty="n", lwd=2, cex=0.5)


lines(sigma, dgamma(1/(sigma^2), shape=3, scale=1), col="red")
# Rate = 1, Scale=1, Shape= 3 seems would work
# It has low densities above the upper bound, and has very lowdensity at zero


# Option 2: Inverse gamma penalty on variance
sigma <- seq(0, 2.5, by = 0.01)
rate <- rev(c(0.5, 1, 2, 3, 5, 7, 10, 13))
per <- seq(90, 10, length.out=length(rate))
mean <- 0.75
plot(sigma, sigma, ylim = c(0,1.4), type = "n", 
     xlab = "Sigma for watershed-area regression", 
     ylab = "Density", main  = "Inverse gamma penalty on variance")
abline(v= sd(log(read.csv("DataIn/ParkenSMSY.csv")$SMSY)), col=grey(0.5), 
       lty="dashed")
abline(v= sd(log(read.csv("DataIn/ParkenSREP.csv")$SREP)), col=grey(0.5), 
       lty="dotted")
# abline(v= sqrt(0.293), col=grey(0.5), lty="dashed")#See Parken et al. (2006)
# abline(v= sqrt(0.146), col=grey(0.5), lty="dashed")#See Parken et al. (2006)

for(i in 1:length(rate)){
  lines(sigma, dgamma(sigma^2, shape = rate[i]*mean, 
                      scale = 1/rate[i])*(sigma^2), 
        lwd=4, col=t_col(color=cols[1], percent= per[i]))
}

legend(x=0, y=1.4, legend=c("sigma ln(SMSY) Parken et al. (upper bound)", 
                             "sigma ln(SREP) Parken et al.(upper bound)"),
       col=grey(0.5),
       lty=c("dashed", "dotted"), bty="n", lwd=2, cex=0.5)

lines(sigma, dgamma(sigma^2, shape=7.5, scale=1/10)*(sigma^2), col="red")
# lines(sigma, dgamma(sigma^2, shape=3.75, scale=1/5)*(sigma^2), col="red")

# Rate = 10, Scale = 1/10, Shape=7.5 works
# More diffuse penalties result in singular results for linear regression 
# (sig goes to 0)

# # Option 3: Gamma penalty on sigma 
# # (doesn't work because it allows zero sigmas)
# sigma <- seq(0, 5, by = 0.1)
# pen <- c(0.5, 1, 3, 5, 7, 10, 13)
# plot(sigma, sigma, ylim = c(0,1.5), type = "n", xlab = "sigma", 
# ylab = "density", main  = "")
# for(peni in pen){
#   lines(sigma, dgamma(sigma, shape = peni*mean, scale = 1/peni), type = "l")
# }

#------------------------------------------------------------------------------

# Updated Ricker sigma penalties
# Option 3: Gamma penalty on sigma
sigma <- seq(0, 2.5, by = 0.01)
rate <- c(13, 10, 8, 5, 3)
per <- seq(90, 10, length.out=length(rate))
mean <- 0.75

plot(x=sigma, y=sigma, type="n", ylab="Density", 
     xlab="Sigma for Ricker model" , 
     main="Gamma penalty on sigma", ylim=c(0,2), xlim=c(0,2.5))
RicSigPark <- read.csv("DataIn/ParkenRicSig.csv")$RicSig 
# Ricker sigma from PSE Chinook SR data extracted 15 Oct 2020. 
# See sigR_metaanalysis.R
RicSig <- exp(read.table("DataOut/PSE_sigma.csv")$Estimate) 
hist_out <- hist(RicSig, plot=FALSE)
hist_outPark <- hist(RicSigPark, plot=FALSE)
lowlimit <- hist_out$breaks[1] # Indicates where bars should start on the plot, 
# assuming units of 0.1
barplot(height =  c( rep(0,lowlimit*10),(hist_out$density)*0.2), width=0.1, 
        col="light blue", border=grey(0.7), space=0, add=TRUE)
#barplot(height= hist_outPark$density*0.2, width=0.1, col=grey(0.8, alpha=0.5), 
# border=grey(0.7, alpha=0.5), space=0, add=TRUE)
RicSig_Thorson <- 0.69 #Thorson et al. 2014 marginal sigma from predictive 
# distribution from hierarchical Ricker model of 20 salmonid stocks from Myers 
# et al 1995
RicSigSD_Thorson <- 0.27
abline(v=RicSig_Thorson, lwd=2, lty="dotted", col=grey(0.8))

legend(x=1, y=1.5, legend=c( #"Distribution of Ricker Sigma Parken", 
                             "Distribution of Ricker Sigma PSE",
                             "Thorson et al. 2014 Sigma"),
       col=c( "light blue"), 
       bty="n", pch=c(#22,
                      22, NA), 
       pt.bg= c(#grey(0.7, alpha=0.5), 
                "light blue", "light grey"), 
       lty=c(#NA, 
             NA, "dashed"), cex=0.5 )


for(i in 1:length(rate)){
  lines(sigma, dgamma(sigma, shape = rate[i]*mean, 
                      scale = 1/rate[i]), 
        lwd=4, col=t_col(color=cols[1], percent= per[i]))
}
# lines(sigma, dgamma(sigma, shape=2.25, scale=0.333), col="red")
# lines(sigma, dgamma(sigma, shape=9.75, scale=0.0769), col="red")
lines(sigma, dgamma(sigma, shape=7.5, scale=0.1), col="red")
# lines(sigma, dgamma(sigma, shape=6, scale=0.125), col="red")
# lines(sigma, dgamma(sigma, shape=3.75, scale=0.2), col="red")


# Suggested shape=7.5 , scale = 1/10, rate = 10, when combined with 
# option 1 for WA regression