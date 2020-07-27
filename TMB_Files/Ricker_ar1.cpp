#include <TMB.hpp >

// Set up Lambert's W function to use to calculate SMSY
// Code taken from https://kaskr.github.io/adcomp/lambert_8cpp_source.html
// Step 1: Code up a plain C version
// Double version of Lambert W function
double LambertW(double x) {
  double logx = log(x);
  double y = (logx > 0 ? logx : 0);
  int niter = 100, i=0;
  for (; i < niter; i++) {
    if ( fabs( logx - log(y) - y) < 1e-9) break;
    y -= (y - exp(logx - y)) / (1 + y);
  }
  if (i == niter) Rf_warning("W: failed convergence");
  return y;
}

TMB_ATOMIC_VECTOR_FUNCTION(
  // ATOMIC_NAME
  LambertW
  ,
  // OUTPUT_DIM
  1,
  // ATOMIC_DOUBLE
  ty[0] = LambertW(tx[0]); // Call the 'double' version
,
// ATOMIC_REVERSE
Type W  = ty[0];                    // Function value from forward pass
Type DW = 1. / (exp(W) * (1. + W)); // Derivative
px[0] = DW * py[0];                 // Reverse mode chain rule
)
  
// Scalar version
template<class Type>
Type LambertW(Type x)
  {
    CppAD::vector<Type> tx(1);
    tx[0] = x;
    return LambertW(tx)[0];
  }
  
  
template <class Type> 
vector <Type> minus_one_to_one(vector <Type> x) 
  { 
    return Type(2) * invlogit(x) - Type(1); 
  } 
  
template<class Type>
Type objective_function<Type>:: operator() ()
{
  DATA_VECTOR(S_ar);
  DATA_VECTOR(logRS_ar);
  DATA_IVECTOR(stkInd_ar);//This is an sequential index
  DATA_IVECTOR(yr_ar);
  //DATA_SCALAR(N_Stks);


  PARAMETER_VECTOR(logA_ar);
  PARAMETER_VECTOR(logB_ar);
  PARAMETER_VECTOR(rho);
  PARAMETER_VECTOR(logSigma_ar);

  
  Type ans=0.0;
  int N_Obs_ar = S_ar.size(); 

  //vector <Type> LogR_Pred_ar(N_Obs_ar);
  vector <Type> LogRS_Pred_ar(N_Obs_ar);
  vector <Type> sigma_ar = exp(logSigma_ar);
  vector <Type> err(N_Obs_ar);
  vector <Type> rho_bounded = minus_one_to_one(rho);
  vector <Type> nLL_ar(N_Obs_ar);

  
  // Ricker AR1 likelihood
  for (int i = 0; i<N_Obs_ar; i++){
    //if (yr_ar(i) == 0) {LogR_Pred_ar(i) = logA_ar(stk_ar(i)) + log(S_ar(i)) - exp(logB_ar(stk_ar(i))) * S_ar(i);
    if (yr_ar(i) == 0) {LogRS_Pred_ar(i) = logA_ar(stkInd_ar(i)) - exp(logB_ar(stkInd_ar(i))) * S_ar(i);
      
        //err(i) = logR_ar(i) - LogR_Pred_ar(i);
        err(i) = logRS_ar(i) - LogRS_Pred_ar(i);
        
      //} else if (yr_ar(i) >= 1) { LogR_Pred_ar(i) = logA_ar(stk_ar(i)) + log(S_ar(i)) - exp(logB_ar(stk_ar(i))) * S_ar(i) + rho(stk_ar(i))*err(i-1);
    } else if (yr_ar(i) >= 1) { LogRS_Pred_ar(i) = logA_ar(stkInd_ar(i)) - exp(logB_ar(stkInd_ar(i))) * S_ar(i) + rho_bounded(stkInd_ar(i))*err(i-1);

        //err(i) = logR_ar(i) - LogR_Pred_ar(i);
        err(i) = logRS_ar(i) - LogRS_Pred_ar(i);
        
        }
    
    //ans += -dnorm(LogR_Pred_ar(i), logR_ar(i),  sigma_ar(stk_ar(i)), true);
    ans += -dnorm(LogRS_Pred_ar(i), logRS_ar(i),  sigma_ar(stkInd_ar(i)), true);
    nLL_ar(i) = -dnorm(LogRS_Pred_ar(i), logRS_ar(i),  sigma_ar(stkInd_ar(i)), true);
  }
  

  // Calculate SMSY and SREP
  int N_stks_ar = logA_ar.size(); 
  vector <Type> SMSY_ar(N_stks_ar);  
  vector <Type> B_ar = exp(logB_ar);
  vector <Type> SREP_ar(N_stks_ar);
  for(int i=0; i<N_stks_ar; i++){
    SMSY_ar(i) =  (1 - LambertW(exp(1-logA_ar(i))) ) / B_ar(i) ;
  }
  SREP_ar = logA_ar / B_ar;

  REPORT(rho_bounded);
  ADREPORT(SMSY_ar);
  ADREPORT(SREP_ar);
  ADREPORT(LogRS_Pred_ar);
  REPORT(nLL_ar);
  return ans;
  
}
