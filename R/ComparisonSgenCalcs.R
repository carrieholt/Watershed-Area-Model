# Comparison of Sgen calculations
# Provided by Paul Van dam-Bates 8 Oct 2024. New versions not implemented in 
# code but will be implemented in IWAM repo and package


## Objective function for Carrie's code.
sGenOptimum <- function ( S, theta ) {
  # Function called from sGenSolver 
  loga <- theta[1]
  b <- theta[2]
  prt <- S * exp( loga - b * S)
  sMSY <- ( 1 - gsl::lambert_W0 (exp ( 1 - loga) ) ) / b
  epsilon <- log(sMSY) - log(prt)
  nlogLike <- - sum( dnorm ( epsilon, 0, 1, log = T))
  return( nlogLike )
}

## Carrie's code for finding SGen
sGenSolverGaussian <- function (loga, b) {
  # Function to estimate Sgen from loga and b Ricker parameters
  theta <- c(loga, b)
  sMSY <- (1 - gsl::lambert_W0(exp(1 - loga))) / b
  fit <- optimize(f = sGenOptimum, interval = c(0, sMSY),
                  theta = theta)
  return(fit$minimum)
}


## Carrie's code for finding SGen updated with nested func and reduced operations.
sGenSolverGaussian2 <- function (loga, b) {
  # Function to estimate Sgen from loga and b Ricker parameters
  sMSY <- (1 - gsl::lambert_W0(exp(1 - loga))) / b
  fn <- function(S){ -sum( dnorm ( log(sMSY) - log(S) - loga + b*S, 0, 1, log = T))}
  fit <- optimize(f = fn, interval = c(0, sMSY))
  return(fit$minimum)
}


## Alternative to Carrie's optimization using an absolute value instead of sum of squares.
sGenSolverAbs = function(loga, b){
  sMSY <- (1-gsl::lambert_W0(exp(1-loga)))/b
  fn <- function(S) { abs(log(sMSY) - log(S) - loga + b*S) }
  fit <- optimize(fn, interval = c(0, sMSY))
  return(fit$minimum)
}

## Using LambertW a second time.
sGenDirect <- function(loga, b){
  sMSY <- ( 1 - gsl::lambert_W0 (exp ( 1 - loga) ) ) / b
  a <- exp(loga)
  -1/b*gsl::lambert_W0(-b*sMSY/a)
}

microbenchmark::microbenchmark(
  new = sGenSolverAbs(loga = 2, b = 0.0001),
  minor_updated = sGenSolverGaussian2(loga = 2, b = 0.0001),
  old = sGenSolverGaussian(loga = 2, b = 0.0001),
  lamw = sGenDirect(loga = 2, b = 0.0001))
