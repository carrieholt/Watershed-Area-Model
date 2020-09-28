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
#compile("TMB_Files/linreg.cpp")
dyn.load(dynlib("TMB_Files/linreg"))
obj <- MakeADFun(data, param, DLL="linreg")
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))#lower=as.vector(lower), upper=as.vector(upper))
pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1
summary(sdreport(obj), p.value=TRUE)

# Look at example simulated output
obj$simulate()$y


# Simulate 1000 predicted data sets (at x fixed values)
Y_Preds <- matrix(nrow = 1000, ncol = length(x))
for(i in 1:1000){
  Y_Preds[i, ] <- obj$simulate()$y
}

# Find the lower and upper 2.5% intervals
Y_Pred_Summ <- apply(Y_Preds, 2, quantile, probs = c(0.025, 0.5, 0.975))

Y_Pred_Summ


# EStiamte prediction intervals using R function predict.lm()
lmtest <- lm(y~x)
predict.lm(lmtest, interval="predict")



# These prediction intervals are larger than those using obj$simulate()

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