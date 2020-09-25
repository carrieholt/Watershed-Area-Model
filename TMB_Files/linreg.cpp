// code from "https://kaskr.github.io/adcomp/_book/Simulation.html"
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  DATA_VECTOR(x);
  PARAMETER(a);
  PARAMETER(b);
  PARAMETER(sd);
  vector<Type> mu = a + b * x;
  Type nll = -sum(dnorm(y, mu, sd, true));
  SIMULATE {
    y = rnorm(mu, sd);  // Simulate response
    REPORT(y);          // Report the simulation
  }
  return nll;
}


