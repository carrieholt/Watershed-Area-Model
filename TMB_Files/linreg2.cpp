// code from "https://kaskr.github.io/adcomp/_book/Simulation.html"
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  DATA_VECTOR(x);
  PARAMETER(a);
  PARAMETER(b);
  PARAMETER(log_sd);
  vector<Type> mu = a + b * x;
  Type nll=0.0;
  nll -= sum(dnorm(y, mu, exp(log_sd), true));
  nll -= log_sd;
  return nll;
}


