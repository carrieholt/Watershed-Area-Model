#include <TMB.hpp >

template<class Type>
Type objective_function<Type>:: operator() ()
{
  DATA_VECTOR(S);
  DATA_VECTOR(logR);
  DATA_IVECTOR(stk);
  DATA_IVECTOR(yr);

  PARAMETER_VECTOR(logA);
  PARAMETER_VECTOR(logB);
  PARAMETER_VECTOR(rho);
  PARAMETER_VECTOR(logSigma);

  
  Type ans=0.0;
  int N_Obs = S.size(); 
  vector <Type> LogR_Pred(N_Obs);
  vector <Type> sigma = exp(logSigma);
  vector <Type> A = exp(logA);
  vector <Type> err(N_Obs);
  
  
  // Ricker likelihood
  for (int i = 0; i<N_Obs; i++){
    if (yr(i) == 0) {LogR_Pred(i) = logA(stk(i)) + log(S(i)) - exp(logB(stk(i))) * S(i);
      err(i) = logR(i) - LogR_Pred(i);
    } else if (yr(i) >= 1) { LogR_Pred(i) = logA(stk(i)) + log(S(i)) - exp(logB(stk(i))) * S(i) + rho(stk(i))*err(i-1);
      err(i) = logR(i) - LogR_Pred(i);
      }
    
    ans += -dnorm(LogR_Pred(i), logR(i),  sigma(stk(i)), true);
  }
  
  
  ADREPORT(A);
  ADREPORT(LogR_Pred);
  
  return ans;
  
}
