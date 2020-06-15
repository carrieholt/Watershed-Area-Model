#include <TMB.hpp >

template<class Type>
Type objective_function<Type>:: operator() ()
{
  DATA_VECTOR(S);
  DATA_VECTOR(logR);
  DATA_IVECTOR(stk);
  DATA_IVECTOR(yr);


  PARAMETER_VECTOR(logA_ar);
  PARAMETER_VECTOR(logB_ar);
  PARAMETER_VECTOR(rho);
  PARAMETER_VECTOR(logSigma_ar);
  PARAMETER_VECTOR(logA_std);
  PARAMETER_VECTOR(logB_std);
  PARAMETER_VECTOR(logSigma_std);
  
  
  Type ans=0.0;
  int N_Obs = S.size(); 
  vector <Type> LogR_Pred_ar(N_Obs);
  vector <Type> sigma_ar = exp(logSigma_ar);
  vector <Type> A_ar = exp(logA_ar);
  vector <Type> err(N_Obs);
  
  
  // Ricker AR1 likelihood
  for (int i = 0; i<N_Obs; i++){
    if (yr(i) == 0) {LogR_Pred_ar(i) = logA_ar(stk(i)) + log(S(i)) - exp(logB_ar(stk(i))) * S(i);
        err(i) = logR(i) - LogR_Pred_ar(i);
        
      } else if (yr(i) >= 1) { LogR_Pred_ar(i) = logA_ar(stk(i)) + log(S(i)) - exp(logB_ar(stk(i))) * S(i) + rho(stk(i))*err(i-1);
        err(i) = logR(i) - LogR_Pred_ar(i);
        }
    
    ans += -dnorm(LogR_Pred_ar(i), logR(i),  sigma_ar(stk(i)), true);
  }
  
  vector <Type> LogR_Pred_std(N_Obs);
  vector <Type> sigma_std = exp(logSigma_std);
  vector <Type> A_std = exp(logA_std);
  
  // Standard Ricker model
  for (int i = 0; i<N_Obs; i++){
    
    LogR_Pred_std(i) = logA_std(stk(i)) + log(S(i)) - exp(logB_std(stk(i))) * S(i);
    ans += -dnorm(LogR_Pred_std(i), logR(i),  sigma_std(stk(i)), true);
    
  }

    
  ADREPORT(A_ar);
  ADREPORT(A_std);
  ADREPORT(LogR_Pred_ar);
  ADREPORT(LogR_Pred_std);
  
  return ans;
  
}
