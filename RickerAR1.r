
# For  RIcker AR1
SRDat <- read.csv("DataIn/SRDat.csv")
SRDat <- SRDat %>% filter (CU_ID==2)

TMB_Inputs <- list(Scale = 1000, logA_Start = 1, rho_Start = 0.1)


#Set up data and parameter lists for input into TMB model
Scale <- TMB_Inputs$Scale
data <- list()
data$S <- SRDat$Spawners/Scale 
data$logR <- log(SRDat$Recruits/Scale)
data$stk <- as.numeric(SRDat$CU_ID)

data$yr <- SRDat$yr_num


param <- list()
param$logA <- rep(TMB_Inputs$logA_Start, N_Stocks)
param$logB <- log(1/( (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/Scale) )
param$rho <- rep(TMB_Inputs$rho_Start, N_Stocks)
param$logSigma <- rep(-2, N_Stocks)


# Compile model if changed:
dyn.unload(dynlib("TMB_Files/Ricker_ar1"))
compile("TMB_Files/Ricker_ar1.cpp")

dyn.load(dynlib("TMB_Files/Ricker_ar1"))

obj <- MakeADFun(data, param, DLL="Ricker_ar1", silent=TRUE)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) 
summary(sdreport(obj))


#R functions for Ricker AR1
RickerAR1.model <- function(theta,R,S){
  a <- exp(as.numeric((theta[1])))
  b <- exp(as.numeric(theta[2]))
  rho <- as.numeric(theta[3])
  sig <- exp(as.numeric(theta[4]))
  n <- length(S)
  logPR <- NA
  dev <- NA
  for (i in 1:n){ 
    if (i==1) {
      logPR[i] <- log(a) + log(S[i]) - b * S[i]
      dev[i] <- log(R[i]) - logPR[i]
      }
    if (i>=2) {
      logPR[i] <- log(a) + log(S[i]) -b * S[i] + rho * dev[i-1]
      dev[i] <- log(R[i]) - logPR[i]
      }
  }
  epsilon.wna <- log(R) - logPR	#residuals
  epsilon <- as.numeric(na.omit(epsilon.wna))
  nloglike <- -sum(dnorm(epsilon,0,sig, log=T))
  #return(list(PR=PR, epsilon=epsilon, nloglike=nloglike)) 
  return(nloglike=nloglike) 
}

RickerAR1.solver <- function(R,S){
  init.vals<-c(1,-2.4,0.1,-2)
  SRfit=optim(par=init.vals,fn=RickerAR1.model,R=R,S=S, method="BFGS", hessian=T)	#CH: hessian=2nd derivative, optim good for up to 40 parameters
  V=solve(SRfit$hessian) #covariance matrix
  std=sqrt(abs(diag(V)))
  X=V/(std%o%std)	#outer product of standard dev matrix (CH comment)
  return(list(SRfit=SRfit, etheta=SRfit$par, V=V, X=X))
}

ch <- RickerAR1.solver(exp(data$logR), data$S)$SRfit$par
ch


# For  Ricker, gives same results at IFRoutTo2013rec.rds for Lower Thompson (CU_ID=2)
SRDat <- read.csv("DataIn/SRDat.csv")
SRDat <- SRDat %>% filter (CU_ID==2)

TMB_Inputs <- list(Scale = 1000, logA_Start = 1)


#Set up data and parameter lists for input into TMB model
Scale <- TMB_Inputs$Scale
data <- list()
data$S <- SRDat$Spawners/Scale 
data$logR <- log(SRDat$Recruits/Scale)
data$stk <- as.numeric(SRDat$CU_ID)

data$yr <- SRDat$yr_num


param <- list()
param$logA <- rep(TMB_Inputs$logA_Start, N_Stocks)
param$logB <- log(1/( (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/Scale) )
param$logSigma <- rep(-2, N_Stocks)


# Compile model if changed:
dyn.unload(dynlib("TMB_Files/Ricker"))
compile("TMB_Files/Ricker.cpp")

dyn.load(dynlib("TMB_Files/Ricker"))

obj <- MakeADFun(data, param, DLL="Ricker", silent=TRUE)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1
summary(sdreport(obj))



# R code for Ricker
Ricker.model <- function(theta,R,S){
  a <- exp(as.numeric((theta[1])))
  b <- exp(as.numeric(theta[2]))
  sig <- exp(as.numeric(theta[3]))
  PR=a*S*exp(-b*S)
  epsilon.wna=log(R)-log(PR)	#residuals
  epsilon=as.numeric(na.omit(epsilon.wna))
  nloglike = -sum(dnorm(epsilon,0,sig, log=T))
  #return(list(PR=PR, epsilon=epsilon, nloglike=nloglike)) #actually returns postive loglikelihood
  return(nloglike=nloglike)
}


Ricker.solver <- function(R,S){
  init.vals<-c(1,-2.4,-2)
  SRfit=optim(par=init.vals,fn=Ricker.model,R=R,S=S, method="BFGS", hessian=T)	#CH: hessian=2nd derivative, optim good for up to 40 parameters
  V=solve(SRfit$hessian) #covariance matrix
  std=sqrt(abs(diag(V)))
  X=V/(std%o%std)	#outer product of standard dev matrix (CH comment)
  return(list(SRfit=SRfit, etheta=SRfit$par, V=V, X=X))
}

ch <- exp(Ricker.solver(exp(data$logR), data$S)$SRfit$par[1:3])
ch
