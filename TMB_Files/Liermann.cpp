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
  Type LambertW(Type x){
    CppAD::vector<Type> tx(1);
    tx[0] = x;
    return LambertW(tx)[0];
  }
  
  
template<class Type>
Type objective_function<Type>:: operator() ()
{
  DATA_VECTOR(S_std);
  DATA_VECTOR(logRS_std);
  DATA_IVECTOR(stk_std);
  DATA_IVECTOR(yr_std);
  DATA_SCALAR(Tau_dist);
  DATA_SCALAR(logMuA_mean);
  DATA_SCALAR(logMuA_sig);
  DATA_SCALAR(Tau_A_dist);
  DATA_VECTOR(WA); 
  DATA_VECTOR(Scale);
  DATA_IVECTOR(Stream);
  DATA_VECTOR(PredlnWA);

  PARAMETER_VECTOR(logA_std);
  PARAMETER_VECTOR(logB_std);
  PARAMETER_VECTOR(logSigma_std);
  PARAMETER(logMuA);
  PARAMETER(logSigmaA);
  PARAMETER(logDelta1);
  PARAMETER(logDelta1ocean);
  PARAMETER(logDelta2);
  PARAMETER(Delta2ocean);
  PARAMETER(logDeltaSigma);
  

  
  Type ans=0.0;
  int N_Obs = S_std.size(); 
  int N_stks_std = Scale.size();
  
  vector <Type> LogRS_Pred_std(N_Obs);
  vector <Type> sigma_std = exp(logSigma_std);
  Type sigmaA = exp(logSigmaA);
  vector <Type> nLL_std(N_Obs);
  
  // Standard Ricker model
  for (int i = 0; i<N_Obs; i++){
    LogRS_Pred_std(i) = logA_std(stk_std(i)) - exp(logB_std(stk_std(i))) * S_std(i);
    ans += -dnorm(LogRS_Pred_std(i), logRS_std(i),  sigma_std(stk_std(i)), true);    
    nLL_std(i) = -dnorm(LogRS_Pred_std(i), logRS_std(i),  sigma_std(stk_std(i)), true);
    
  }

    // Add hierarchical structure to A ==============
  for(int i=0; i<N_stks_std; i++){
    // add prior on logA
    ans += -dnorm(logA_std(i), logMuA, sigmaA, true );
    // add prior on sigma 
    ans += -dgamma(pow(sigma_std(i),-2), Tau_dist, 1/Tau_dist, true);
  }
  
  // Add priors for hyperpars ====================
  // MuA prior
  ans += -dnorm(logMuA, logMuA_mean, logMuA_sig, true);
  // sigmaA prior
  ans += -dgamma(pow(sigmaA,-2), Tau_A_dist, 1/Tau_A_dist, true);

  //Calculate SMSY and SREP
  vector <Type> SMSY_std(N_stks_std);  
  vector <Type> SREP_std(N_stks_std);  
  
  for(int i=0; i<N_stks_std; i++){
    SMSY_std(i) =  (1 - LambertW(exp(1-logA_std(i))) ) / exp(logB_std(i)) ;
  }
  SREP_std = logA_std / exp(logB_std);
  
  
  //Liermann's model with both stream and ocean type=================
  vector <Type> PredlnSMSY(N_stks_std);
  Type sigma_delta = exp(logDeltaSigma);
  
  for (int i=0; i<N_stks_std; i++){
    PredlnSMSY(i) = logDelta1 + logDelta1ocean * Stream(i) + ( exp(logDelta2) + Delta2ocean * Stream(i) ) * log(WA(i)) ;
    ans += -dnorm( PredlnSMSY(i), log(SMSY_std(i) * Scale(i) ),  sigma_delta, true);
  }
  
  // Get predicted values for plotting  WA regresssion with CIs
  int N_pred = PredlnWA.size();
  vector <Type> PredlnSMSYs_CI(N_pred);
  vector <Type> PredlnSMSYo_CI(N_pred);

  for (int i=0; i<N_pred; i++){
    PredlnSMSYs_CI(i) = logDelta1 + exp(logDelta2) * PredlnWA(i);
    PredlnSMSYo_CI(i) = logDelta1 + logDelta1ocean + (exp(logDelta2) + Delta2ocean) * PredlnWA(i);
  }
  
  ADREPORT(SMSY_std);
  ADREPORT(SREP_std);
  ADREPORT(LogRS_Pred_std);
  ADREPORT(PredlnSMSYs_CI);
  ADREPORT(PredlnSMSYo_CI);
  REPORT(nLL_std);
  return ans;
  
}

