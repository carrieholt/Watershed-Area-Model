  test <- seq(0.001,3,0.1)
  plot(x=test, y=abs(dcauchy(test,0,1)), type="l", ylab="Probability density", xlab="" , ylim=c(0,0.4))
  
  # Cauchy = normal/sqrt(chi^2), Gelman et al. 2006
  #theta <- xi*eta
  
  prior.scale <- 1
  #tau.xi <- prior.scale^(-2)
  #xi <- dnorm (test, 0, tau.xi)
  xi <- dnorm (test, 0, prior.scale)
  tau.eta <-  dgamma (test,0.5,0.5) # chi^2 with 1 d.f.
  plot(x=test, y=tau.eta, type="l")
  sigma.theta <- abs(xi)/sqrt(tau.eta)
  plot(x=test, y=sigma.theta, type="l")
