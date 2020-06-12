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
  PARAMETER(logSigma);
  PARAMETER(rho);
  
  
  Type ans=0.0;
  int N_Obs = S.size(); 
  vector <Type> LogR_Pred(N_Obs);
  Type sigma = exp(logSigma);
  Type A = exp(logA);
  vector <Type> err(N_Obs);
  
  
  // Ricker likelihood
  for (int i = 0; i<N_Obs; i++){
    if (i == 0) { err(i) = 0;
    } else if (i >= 1) { err(i) = err(i-1) * rho;
      }
    LogR_Pred(i) = logA + log(S(i)) - exp(logB) * S(i) + err(i);
    ans += -dnorm(LogR_Pred(i), logR(i),  sigma, true);
  }
  
  
  ADREPORT(A);
  ADREPORT(LogR_Pred);
  
  return ans;
  
}
