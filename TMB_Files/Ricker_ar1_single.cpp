#include <TMB.hpp >

template<class Type>
Type objective_function<Type>:: operator() ()
{
  DATA_VECTOR(S);
  DATA_VECTOR(logR);
  DATA_IVECTOR(stk);
  DATA_IVECTOR(yr);

  PARAMETER(logA);
  PARAMETER(logB);
  PARAMETER(rho);
  PARAMETER(logSigma);
  
  
  Type ans=0.0;
  int N_Obs = S.size(); 
  vector <Type> LogR_Pred(N_Obs);
  Type sigma = exp(logSigma);
  Type A = exp(logA);
  vector <Type> err(N_Obs);
  
  
  // Ricker likelihood
  for (int i = 0; i<N_Obs; i++){
    if (i == 0) {LogR_Pred(i) = logA + log(S(i)) - exp(logB) * S(i);
      err(i) = logR(i) - LogR_Pred(i);
    } else if (i >= 1) { LogR_Pred(i) = logA + log(S(i)) - exp(logB) * S(i) + rho*err(i-1);
      err(i) = logR(i) - LogR_Pred(i);
      }
    
    //ans += -dnorm(LogR_Pred(i), logR(i),  sigma, true);
    ans += -dnorm(logR(i), LogR_Pred(i),  sigma, true);
  }
  
  
  ADREPORT(A);
  ADREPORT(LogR_Pred);
  
  return ans;
  
}
