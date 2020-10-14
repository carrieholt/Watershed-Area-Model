# Code to compare methods for estimating prediction intervals using a simple linear regression
library(TMB)

# Create Data
y=c(2,1,4,2,8,5,8,10,7,9)
x=1:10

# Look at data
plot(x,y)
abline(lm(y~x)$coef[1], lm(y~x)$coef[2])


# Read in data to TMB
data <- list()
data$x <- x
data$y <- y
param <- list()
param$a <- 1
param$b <- 1
param$sd <- 1


# Run TMB 
#compile("linreg.cpp")
dyn.load(dynlib("linreg"))
obj <- MakeADFun(data, param, DLL="linreg")
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))#lower=as.vector(lower), upper=as.vector(upper))
pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1
summary(sdreport(obj), p.value=TRUE)

# Simulate 1000 predicted data sets (at x fixed values)
Y_Preds <- matrix(nrow = 1000, ncol = length(x))
for(i in 1:1000){
  Y_Preds[i, ] <- obj$simulate()$y
}

# Find the lower and upper 2.5% intervals
Y_Pred_Summ <- apply(Y_Preds, 2, quantile, probs = c(0.025, 0.5, 0.975))

Y_Pred_Summ

# simulate from multivariate normal
sims <- rmvnorm(1000, unlist(pl), solve(obj$he()))
# now for each simulate estimates
Preds <- matrix(nrow = 1000, ncol = length(x))
for(i in 1:dim(sims)[1]){
  #Preds[i, ] <- rnorm(length(x), sims[i, 1] + sims[i, 2]*x, sims[i, 3])
  #Preds[i, ] <- sims[i, 1] + sims[i, 2]*x
  Preds[i, ] <- rnorm(length(x), sims[i, 1] + sims[i, 2]*x, pl$sd)
}
apply(Preds, 2, quantile, probs = c(0.025, 0.5, 0.975))


# EStiamte prediction intervals using R function predict.lm()
lmtest <- lm(y~x)
predict.lm(lmtest, interval="predict")



# These prediction intervals are larger than those using obj$simulate()

# Carrie code to calculate prediction intervals that match predict.lm:

# STEYX from https://support.microsoft.com/en-us/office/steyx-function-6ce74b2c-449d-4a6e-b9ac-f9cef5ba48ab
sumXErr <- sum( (x-mean(x))^2 )
sumYErr <- sum( (y-mean(y))^2 )
sumXYErr <- sum( (y - mean(y)) * (x - mean(x)) ) ^2
STEYX <- sqrt( (1/(length(x)-2)) * (sumYErr - sumXYErr/sumXErr) )

# SE of the prediction from http://www.real-statistics.com/regression/confidence-and-prediction-intervals/
SE.pred <- STEYX * sqrt( (1 + 1/length(x) + ((x-mean(x))^2)/sumXErr) ) 
t.crit <- qt(0.975,df=length(x)-2) #95% intervals

upr <- predict(lmtest) + SE.pred*t.crit
upr

lwr <- predict(lmtest) - SE.pred*t.crit
lwr

# Compare to R function. They match
predict.lm(lmtest, interval="predict")


#----------------------------------------------------------------------------
# A function to do calculate prediction intervals
# x = independent variable
# y = dependenet variable
# Newx = x variables for which you'd like the prediction intervals
# Predy = Predicted y variable at those Newx's


PredInt <- function(x,y,Newx=x, Predy){
  sumXErr <- sum( (x-mean(x))^2 )
  sumYErr <- sum( (y-mean(y))^2 )
  sumXYErr <- sum( (y - mean(y)) * (x - mean(x)) ) ^2
  STEYX <- sqrt( (1/(length(x)-2)) * (sumYErr - sumXYErr/sumXErr) )
  
  # SE of the prediction from http://www.real-statistics.com/regression/confidence-and-prediction-intervals/
  SE.pred <- STEYX * sqrt( (1 + 1/length(x) + ((Newx - mean(x))^2)/sumXErr) ) 
  t.crit <- qt(0.975,df=length(x)-2) #95% intervals
  
  upr <- Predy + SE.pred*t.crit
  lwr <- Predy - SE.pred*t.crit
  PI <- list()
  PI$upr <- upr
  PI$lwr <- lwr
  return(PI)
  
}

PredInt(x=x,y=y,Predy=predict(lmtest))



# try with using loga, logb, logsd -- gets even wonkier results
#compile("linreg_log.cpp")
dyn.load(dynlib("linreg_log"))

data$JA = 0
param <- list()
param$loga <- 0
param$logb <- 0
param$logsd <- 0


obj_log <- MakeADFun(data, param, DLL="linreg_log", hessian = T)
opt_log <- nlminb(obj_log$par, obj_log$fn, obj_log$gr, 
                  control = list(eval.max = 1e5, iter.max = 1e5))    #lower=as.vector(lower), upper=as.vector(upper))
pl_log <- obj_log$env$parList(opt_log$par) # Parameter estimate after phase 1
summary(sdreport(obj_log), p.value=TRUE)

Y_Preds_log <- matrix(nrow = 1000, ncol = length(x))
for(i in 1:1000){
  Y_Preds_log[i, ] <- obj_log$simulate()$y
}
apply(Y_Preds_log, 2, quantile, probs = c(0.025, 0.5, 0.975))

# what if simulate combinations using mvnorm 
library(mvtnorm)
# getting negative values when do this -- should have logA, logB, logsd as params
sims <- rmvnorm(1000, exp(unlist(pl_log)), solve(obj_log$he()))
# now for each simulate estimates
Preds <- matrix(nrow = 1000, ncol = length(x))
for(i in 1:dim(sims)[1]){
  #Preds[i, ] <- exp(sims[i, 1]) + exp(sims[i, 2])*x
  #Preds[i, ] <- rnorm(length(x), exp(sims[i, 1]) + exp(sims[i, 2])*x, exp(sims[i, 3]))
  Preds[i, ] <- sims[i, 1] + sims[i, 2]*x
}
apply(Preds, 2, quantile, probs = c(0.025, 0.5, 0.975))

# nope that didn't work!