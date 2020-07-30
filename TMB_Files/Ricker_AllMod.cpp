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
  DATA_VECTOR(S_std);
  //DATA_VECTOR(logR_std);
  DATA_VECTOR(logRS_std);
  
  // I need to include both ind_std and stk_std, and use ind_stk in the estimation step for loops, but stk_std in the model compilation step at bottom.
  DATA_IVECTOR(stk_std);
  DATA_IVECTOR(yr_std);
  DATA_VECTOR(S_ar);
  //DATA_VECTOR(logR_ar);
  DATA_VECTOR(logRS_ar);
  DATA_IVECTOR(stk_ar);
  DATA_IVECTOR(yr_ar);
  DATA_VECTOR(S_surv);
  //DATA_VECTOR(logR_surv);
  DATA_VECTOR(logRS_surv);
  DATA_IVECTOR(stk_surv);
  DATA_IVECTOR(yr_surv);
  DATA_VECTOR(Surv_surv);
  DATA_VECTOR(MeanLogSurv_surv);

  DATA_VECTOR(WA);
  DATA_VECTOR(Scale);
  //DATA_SCALAR(Tau_dist);
  DATA_VECTOR(Stream);
  DATA_INTEGER(N_stream);
  DATA_INTEGER(N_ocean);

  
  PARAMETER_VECTOR(logA_std);
  PARAMETER_VECTOR(logB_std);
  PARAMETER_VECTOR(logSigma_std);
  PARAMETER_VECTOR(logA_ar);
  PARAMETER_VECTOR(logB_ar);
  PARAMETER_VECTOR(rho);
  PARAMETER_VECTOR(logSigma_ar);
  PARAMETER_VECTOR(logA_surv);
  PARAMETER_VECTOR(logB_surv);
  PARAMETER_VECTOR(logSigma_surv);
  PARAMETER_VECTOR(gamma);
  PARAMETER(logDelta1);
  PARAMETER(Delta2);
  PARAMETER(logDeltaSigma);
  //PARAMETER(slogDelta1);
  //PARAMETER(sDelta2);
  //PARAMETER(slogDeltaSigma);
  //PARAMETER_VECTOR(logSgen);
  
  
  Type ans=0.0;
  int N_Obs_std = S_std.size(); 
  int N_Obs_ar = S_ar.size(); 
  int N_Obs_surv = S_surv.size(); 
  
  //vector <Type> LogR_Pred_ar(N_Obs_ar);
  vector <Type> LogRS_Pred_ar(N_Obs_ar);
  vector <Type> sigma_ar = exp(logSigma_ar);
  vector <Type> err(N_Obs_ar);
  vector <Type> rho_bounded = minus_one_to_one(rho);
  vector <Type> nLL_ar(N_Obs_ar);

  
  // Ricker AR1 likelihood
  for (int i = 0; i<N_Obs_ar; i++){
    //if (yr_ar(i) == 0) {LogR_Pred_ar(i) = logA_ar(stk_ar(i)) + log(S_ar(i)) - exp(logB_ar(stk_ar(i))) * S_ar(i);
    if (yr_ar(i) == 0) {LogRS_Pred_ar(i) = logA_ar(stk_ar(i)) - exp(logB_ar(stk_ar(i))) * S_ar(i);
      
        //err(i) = logR_ar(i) - LogR_Pred_ar(i);
        err(i) = logRS_ar(i) - LogRS_Pred_ar(i);
        
      //} else if (yr_ar(i) >= 1) { LogR_Pred_ar(i) = logA_ar(stk_ar(i)) + log(S_ar(i)) - exp(logB_ar(stk_ar(i))) * S_ar(i) + rho(stk_ar(i))*err(i-1);
    } else if (yr_ar(i) >= 1) { LogRS_Pred_ar(i) = logA_ar(stk_ar(i)) - exp(logB_ar(stk_ar(i))) * S_ar(i) + rho_bounded(stk_ar(i))*err(i-1);

        //err(i) = logR_ar(i) - LogR_Pred_ar(i);
        err(i) = logRS_ar(i) - LogRS_Pred_ar(i);
        
        }
    
    //ans += -dnorm(LogR_Pred_ar(i), logR_ar(i),  sigma_ar(stk_ar(i)), true);
    ans += -dnorm(LogRS_Pred_ar(i), logRS_ar(i),  sigma_ar(stk_ar(i)), true);
    nLL_ar(i) = -dnorm(LogRS_Pred_ar(i), logRS_ar(i),  sigma_ar(stk_ar(i)), true);
  }
  
  //vector <Type> LogR_Pred_std(N_Obs_std);
  vector <Type> LogRS_Pred_std(N_Obs_std);
  vector <Type> sigma_std = exp(logSigma_std);
  vector <Type> nLL_std(N_Obs_std);
  
  // Standard Ricker model
  for (int i = 0; i<N_Obs_std; i++){
    
    //LogR_Pred_std(i) = logA_std(stk_std(i)) + log(S_std(i)) - exp(logB_std(stk_std(i))) * S_std(i);
    LogRS_Pred_std(i) = logA_std(stk_std(i)) - exp(logB_std(stk_std(i))) * S_std(i);
    
    //ans += -dnorm(LogR_Pred_std(i), logR_std(i),  sigma_std(stk_std(i)), true);
    ans += -dnorm(LogRS_Pred_std(i), logRS_std(i),  sigma_std(stk_std(i)), true);
    nLL_std(i) = -dnorm(LogRS_Pred_std(i), logRS_std(i),  sigma_std(stk_std(i)), true);
    
  }

  //vector <Type> LogR_Pred_surv(N_Obs_surv);
  vector <Type> LogRS_Pred_surv(N_Obs_surv);
  vector <Type> sigma_surv = exp(logSigma_surv);
  //Type sigma_surv = exp(logSigma_surv);
  vector <Type> nLL_surv(N_Obs_surv);
  
  // Ricker likelihood with survival covariate
  for (int i = 0; i<N_Obs_surv; i++){
    //LogR_Pred_surv(i) = logA_surv + log(S_surv(i)) - exp(logB_surv) * S_surv(i) + gamma * Surv_surv(i);
    //LogRS_Pred_surv(i) = logA_surv - exp(logB_surv) * S_surv(i) + gamma * Surv_surv(i);
    LogRS_Pred_surv(i) = logA_surv(stk_surv(i)) - exp(logB_surv(stk_surv(i))) * S_surv(i) + gamma(stk_surv(i)) * Surv_surv(i);
    
    //ans += -dnorm(LogR_Pred_surv(i), logR_surv(i),  sigma_surv, true);
    ans += -dnorm(LogRS_Pred_surv(i), logRS_surv(i),  sigma_surv(stk_surv(i)), true);
    nLL_surv(i) = -dnorm(LogRS_Pred_surv(i), logRS_surv(i),  sigma_surv(stk_surv(i)), true);
  }
  
  
  //Compile parameters
  //int N_stks = logA_std.size(); 
  //vector <Type> logA_(N_stks);
  //vector <Type> logAadj_(N_stks);
  //vector <Type> logB_(N_stks);
  //vector <Type> logSigma_(N_stks);
  //vector <Type> SMSYadj(N_stks);  
  //vector <Type> SMSY(N_stks);  
  
  //for (int i = 0; i<N_stks; i++){
  //  if(model(i) == 0) {
  //    logA_(i) = logA_std(i);
  //    logB_(i) = logB_std(i);
  //    logSigma_(i) = logSigma_std(i);
  //    logAadj_(i) = logA_std(i) + pow( exp( logSigma_std(i) ), 2) /2;
      
   // }    
   // if(model(i) == 1) {
   //   logA_(i) = logA_ar(i);
   //   logB_(i) = logB_ar(i);
   //   logSigma_(i) = logSigma_ar(i);
   //   logAadj_(i) = logA_ar(i) + pow( exp( logSigma_ar(i) ), 2) /2;
   // }    
    
 // }
  

  // Calculate SMSY using Lambert's W function
  // Approach from Scheurell 2016
  //vector <Type> B = exp(logB_);
  //for(int i=0; i<N_stks; i++){
  //  SMSYadj[i] =  (1 - LambertW(exp(1-logAadj_[i])) ) / B[i] ;
  //  SMSY[i] =  (1 - LambertW(exp(1-logA_[i])) ) / B[i] ;
  //}
  
  // Calculate SREP
  //vector <Type> SREPadj(N_stks);
  //SREPadj = logAadj_ / B;
  
  // Now estimate Sgen
  //vector <Type> LogSMSY(N_stks);
  //vector <Type> Sgen = exp(logSgen);
 
  //LogSMSY = logA_ + logSgen - B * Sgen;
  //vector <Type> Diff = exp(LogSMSY)-SMSY;
  //ans += -sum(dnorm(Diff, 0, Sgen_sig, true ));
  
  // Code for SMSY SREP without AR model
  int N_stks_std = logA_std.size(); 
  int N_stks_ar = logA_ar.size(); 
  int N_stks_surv = logA_surv.size(); 
  //int N_stks_surv = logA_surv.size(); 
  //vector <Type> logAadj_(N_stks_std);
  //vector <Type> SMSYadj(N_stks_std);  
  vector <Type> SMSY_std(N_stks_std);  
  vector <Type> B_std = exp(logB_std);
  vector <Type> SREP_std(N_stks_std);
  vector <Type> SMSY_ar(N_stks_ar);  
  vector <Type> B_ar = exp(logB_ar);
  vector <Type> SREP_ar(N_stks_ar);
  //Type SMSY_surv;  
  //Type B_surv = exp(logB_surv);
  //Type SREP_surv;
  vector <Type> SMSY_surv(N_stks_surv);  
  vector <Type> B_surv = exp(logB_surv);
  vector <Type> SREP_surv(N_stks_surv);
  
  for(int i=0; i<N_stks_std; i++){
    SMSY_std(i) =  (1 - LambertW(exp(1-logA_std(i))) ) / B_std(i) ;
  }
  SREP_std = logA_std / B_std;
  
  for(int i=0; i<N_stks_ar; i++){
    SMSY_ar(i) =  (1 - LambertW(exp(1-logA_ar(i))) ) / B_ar(i) ;
  }
  SREP_ar = logA_ar / B_ar;
  
  for(int i=0; i<N_stks_surv; i++){
    SMSY_surv(i) =  (1 - LambertW(exp(1 - (logA_surv(i) + gamma(i) * MeanLogSurv_surv(i)))) )/ B_surv(i) ;
  }
  SREP_surv = logA_surv / B_surv;
  
  
  // Estimate watershed-area regression.
  // First aggregate SMSY and SREP in model order
  int N_stks = N_stks_std + N_stks_ar + N_stks_surv;
  vector <Type> SMSY(N_stks);
  vector <Type> SREP(N_stks);
  
  for(int i=0; i < N_stks_std; i++){
    SMSY[i] = SMSY_std[i];
    SREP[i] = SREP_std[i];
  }
  for(int i=0; i < N_stks_ar; i++){
    SMSY[N_stks_std + i] = SMSY_ar[i];
    SREP[N_stks_std + i] = SREP_ar[i];
  }
  for(int i=0; i < N_stks_surv; i++){
    SMSY[N_stks_std + N_stks_ar + i] = SMSY_surv[i];
    SREP[N_stks_std + N_stks_ar + i] = SREP_surv[i];
  }
  
//Separate into stream and ocean type stocks

  vector <Type> SMSY_stream(N_stream);
  vector <Type> SREP_stream(N_stream);
  vector <Type> WA_stream(N_stream);
  vector <Type> Scale_stream(N_stream);
  vector <Type> SMSY_ocean(N_ocean);
  vector <Type> SREP_ocean(N_ocean);
  vector <Type> WA_ocean(N_ocean);
  vector <Type> Scale_ocean(N_ocean);
  

  int j = 0;
  int k = 0;
  
  for(int ii=0; ii < N_stks; ii++){
    if(Stream[ii]==1){
      SMSY_stream[j] = SMSY[ii];
      SREP_stream[j] = SREP[ii];
      WA_stream[j] = WA[ii];
      Scale_stream[j] = Scale[ii];
      j += 1;
    }
    if(Stream[ii]==2){
      SMSY_ocean[k] = SMSY[ii];
      SREP_ocean[k] = SREP[ii];
      WA_ocean[k] = WA[ii];
      Scale_ocean[k] = Scale[ii];
      k += 1;
    }
    
  }


  vector <Type> PredlnSMSY(N_stks);
  Type Delta2_bounded = invlogit(Delta2);
  Type sigma_delta = exp(logDeltaSigma);
  
  for (int i=0; i<N_stks; i++){
    //PredlnSMSY(i) = logDelta1 + exp(logDelta2) * log(WA(i));
    PredlnSMSY(i) = logDelta1 + Delta2_bounded * log(WA(i));
    ans += -dnorm(PredlnSMSY(i), log(SMSY(i)*Scale(i)),  sigma_delta, true);
  }
  // Add Inverse gamma prior on sigma_delta^2
  //ans += -dgamma(pow(sigma_delta,-2), Tau_dist, 1/Tau_dist, true);

  //vector <Type> sPredlnSMSY(N_stream);
  //Type sDelta2_bounded = invlogit(sDelta2);
  //Type ssigma_delta = exp(slogDeltaSigma);
  
  //for (int i=0; i<N_stream; i++){
    //sPredlnSMSY(i) = slogDelta1 + sDelta2_bounded * log(WA_stream(i));
    //ans += -dnorm(sPredlnSMSY(i), log(SMSY_stream(i)*Scale_stream(i)),  ssigma_delta, true);
  //}
  
  
  //ADREPORT(A_ar);
  //ADREPORT(A_std);
  //ADREPORT(logA_);
  //ADREPORT(logAadj_);
  //ADREPORT(SMSYadj);
  REPORT(rho_bounded);
  ADREPORT(SMSY_std);
  ADREPORT(SREP_std);
  ADREPORT(SMSY_ar);
  ADREPORT(SREP_ar);
  ADREPORT(SMSY_surv);
  ADREPORT(SREP_surv);
  ADREPORT(SMSY);
  ADREPORT(SMSY_stream)
  ADREPORT(SREP_stream)
  ADREPORT(WA_stream)
  ADREPORT(logDelta1)
  //ADREPORT(logDelta2)
  ADREPORT(Delta2_bounded)
  ADREPORT(sigma_delta)
  //ADREPORT(slogDelta1)
  //ADREPORT(sDelta2_bounded)
  //ADREPORT(ssigma_delta)
  //ADREPORT(gamma);
  //ADREPORT(LogR_Pred_ar);
  //ADREPORT(LogR_Pred_std);
  //ADREPORT(LogR_Pred_surv);
  ADREPORT(LogRS_Pred_ar);
  ADREPORT(LogRS_Pred_std);
  ADREPORT(LogRS_Pred_surv);
  REPORT(nLL_std);
  REPORT(nLL_ar);
  REPORT(nLL_surv);
  REPORT(ans);
  //ADREPORT(Sgen);
  return ans;
  
}
