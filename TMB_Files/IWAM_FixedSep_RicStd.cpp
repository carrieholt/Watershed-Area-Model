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
  DATA_VECTOR(logRS_std);
  
  // I need to include both ind_std and stk_std, and use ind_stk in the estimation step for loops, but stk_std in the model compilation step at bottom.
  DATA_IVECTOR(stk_std);
  DATA_IVECTOR(yr_std);
  //DATA_SCALAR(Tau_sig);
  
  DATA_VECTOR(WA);
  DATA_VECTOR(Scale);
  DATA_IVECTOR(Stream);
  DATA_INTEGER(N_stream);
  DATA_INTEGER(N_ocean);
  //DATA_IVECTOR(order_noChick);

  //Hierarchical hyper pars
  //DATA_SCALAR(logMuDelta1_mean);
  //DATA_SCALAR(logMuDelta1_sig);
  //DATA_SCALAR(logMuDelta2_mean);
  //DATA_SCALAR(logMuDelta2_sig);
  //DATA_SCALAR(Tau_Delta1_dist);
  //DATA_SCALAR(Tau_Delta2_dist);
  
  DATA_VECTOR(PredlnWA);
  
  PARAMETER_VECTOR(logA_std);
  PARAMETER_VECTOR(logB_std);
  PARAMETER_VECTOR(logSigma_std);
  PARAMETER(logDelta1);
  PARAMETER(logDelta1ocean);
  PARAMETER(logDelta2);
  PARAMETER(logDelta2ocean);
  PARAMETER(logDeltaSigma);

  // Separate stream and ocean type Deltas- fixed effects
  //PARAMETER(slogDelta1);
  //PARAMETER(sDelta2);
  //PARAMETER(slogDeltaSigma);
  //PARAMETER(ologDelta1);
  //PARAMETER(ologDelta2);
  //PARAMETER(ologDeltaSigma);
  //PARAMETER_VECTOR(logSgen);
  ////Lierman WA model pars (gives same resutls as above)
  //PARAMETER(logDelta1);
  //PARAMETER(logDelta1ocean);
  //PARAMETER(logDelta2);
  //PARAMETER(logDelta2ocean);
  //PARAMETER(logDeltaSigma);
  
  ////Hierarchical pars
  //PARAMETER_VECTOR(logDelta1);
  //PARAMETER(logDelta1);
  //PARAMETER_VECTOR(logDelta2);
  //PARAMETER(logDeltaSigma);
  //PARAMETER(logMuDelta1);
  //PARAMETER(SigmaDelta1);
  //PARAMETER(logMuDelta2);
  //PARAMETER(SigmaDelta2);
  
  Type ans=0.0;
  int N_Obs_std = S_std.size(); 

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
  //// Add inverse gamma penalty/prior on sigma_std
  int N_stks_std = logA_std.size(); 
  //for (int i = 0; i<N_stks_std; i++){
  //  ans += -dgamma(pow(sigma_std(i),-2), Tau_sig, 1/Tau_sig, true);
  //}
  
  // Code for SMSY SREP 
  vector <Type> SMSY_std(N_stks_std);  
  vector <Type> B = exp(logB_std);
  vector <Type> SREP_std(N_stks_std);

  for(int i=0; i<N_stks_std; i++){
    SMSY_std(i) =  (1 - LambertW(exp(1-logA_std(i))) ) / B(i) ;
  }
  SREP_std = logA_std / B;
  

  
  
//Separate into stream and ocean type stocks

  vector <Type> SMSY_stream(N_stream);
  vector <Type> SREP_stream(N_stream);
  vector <Type> Scale_stream(N_stream);
  vector <Type> SMSY_ocean(N_ocean);
  vector <Type> SREP_ocean(N_ocean);
  vector <Type> Scale_ocean(N_ocean);
  

  int j = 0;
  int k = 0;
  
  for(int ii=0; ii < N_stks_std; ii++){
    if(Stream[ii]==0){
      SMSY_stream[j] = SMSY_std[ii];
      SREP_stream[j] = SREP_std[ii];
      //WA_stream[j] = WA[ii];
      Scale_stream[j] = Scale[ii];
      j += 1;
    }
    if(Stream[ii]==1){
      SMSY_ocean[k] = SMSY_std[ii];
      SREP_ocean[k] = SREP_std[ii];
      //WA_ocean[k] = WA[ii];
      Scale_ocean[k] = Scale[ii];
      k += 1;
    }
    
  }

  //// WA model with all data =============
  //vector <Type> PredlnSMSY(N_stks);
  //Type Delta2_bounded = invlogit(Delta2);
  //Type sigma_delta = exp(logDeltaSigma);
  
//  for (int i=0; i<N_stks; i++){
  //  //PredlnSMSY(i) = logDelta1 + exp(logDelta2) * log(WA(i));
  //  PredlnSMSY(i) = logDelta1 + Delta2_bounded * log(WA(i));
  //  ans += -dnorm(PredlnSMSY(i), log(SMSY(i)*Scale(i)),  sigma_delta, true);
//  }

  //// Add Inverse gamma prior on sigma_delta^2
  ////ans += -dgamma(pow(sigma_delta,-2), Tau_dist, 1/Tau_dist, true);

  // WA model with only stream-type data =============
  //vector <Type> sPredlnSMSY(N_stream);
  //Type sDelta2_bounded = invlogit(sDelta2);
  //Type ssigma_delta = exp(slogDeltaSigma);
  
  //for (int i=0; i<N_stream; i++){
  //  sPredlnSMSY(i) = slogDelta1 + sDelta2_bounded * log(WA_stream(i));
  //  ans += -dnorm(sPredlnSMSY(i), log(SMSY_stream(i)*Scale_stream(i)),  ssigma_delta, true);
  //}
  
  // WA model with only ocean-type data =============
  //vector <Type> oPredlnSMSY(N_ocean);
  ////Type oDelta2_bounded = invlogit(oDelta2);
  //Type osigma_delta = exp(ologDeltaSigma);
  
  //for (int i=0; i<N_ocean; i++){
  //  oPredlnSMSY(i) = ologDelta1 + exp(ologDelta2) * log(WA_ocean(i));
  //  ans += -dnorm(oPredlnSMSY(i), log(SMSY_ocean(i)*Scale_ocean(i)),  osigma_delta, true);
  //}
  
  //Liermann's model with both stream and ocean type=================
  //int N_stks_short = order_noChick.size();
  vector <Type> PredlnSMSY(N_stks_std);
  Type sigma_delta = exp(logDeltaSigma);
  
  for (int i=0; i<N_stks_std; i++){
  //for (int i=0; i<N_stks_short; i++){
    PredlnSMSY(i) = logDelta1 + logDelta1ocean * Stream(i) + ( exp(logDelta2) + exp(logDelta2ocean) * Stream(i) ) * log(WA(i)) ;
    ans += -dnorm( PredlnSMSY(i), log(SMSY_std(i) * Scale(i) ),  sigma_delta, true);
  }
  
  ////Model with stream and ocean-types add random slope (logDelta1) and yi-intercept (Delta2) ==============
  //vector <Type> PredlnSMSY(N_stks);
  //Type sigma_delta = exp(logDeltaSigma);
  
  //for (int i=0; i<N_stks; i++){
    //PredlnSMSY(i) = logDelta1(Stream(i)) + exp(logDelta2(Stream(i))) * log(WA(i));
    //PredlnSMSY(i) = logDelta1 + exp(logDelta2(Stream(i))) * log(WA(i));
    //ans += -dnorm(PredlnSMSY(i), log(SMSY(i)*Scale(i)),  sigma_delta, true);//sigma_delta(Stream(i)), true);
  //}
  
  //// Add hierarchical structure to logDelta1, logDelta2, sigmaDelta==============
  //int n_lh = 2;
  //for(int i=0; i< n_lh; i++){
    //// add prior on logDelta1
    //ans += -dnorm(logDelta1(i), logMuDelta1, SigmaDelta1, true );
    // add prior on logDelta2
    //ans += -dnorm(logDelta2(i), logMuDelta2, SigmaDelta2, true );
    //// add prior on sigma
    ////ans += -dgamma(pow(sigmaDelta1(i),-2), Tau_dist, 1/Tau_dist, true);
  //}
  
  // Add priors for hyperpars ====================
  //// MuDelta1 prior
  //ans += -dnorm(logMuDelta1, logMuDelta1_mean, logMuDelta1_sig, true);
  //// SigmaDelta1 prior
  //ans += -dgamma(pow(SigmaDelta1,-2), Tau_Delta1_dist, 1/Tau_Delta1_dist, true);
  //// MuDelta2 prior
  //ans += -dnorm(logMuDelta2, logMuDelta2_mean, logMuDelta2_sig, true);
  // SigmaDelta2 prior
  //ans += -dgamma(pow(SigmaDelta2,-2), Tau_Delta2_dist, 1/Tau_Delta2_dist, true);
  
  // Get predicted values for plotting  WA regresssion with CIs
  int N_pred = PredlnWA.size();
  vector <Type> PredlnSMSYs_CI(N_pred);
  vector <Type> PredlnSMSYo_CI(N_pred);
  //vector <Type> PredlnSMSY_S(N_pred);
  //vector <Type> PredlnSMSY_O(N_pred);
  
  for (int i=0; i<N_pred; i++){
    PredlnSMSYs_CI(i) = logDelta1 + exp(logDelta2) * PredlnWA(i);
    PredlnSMSYo_CI(i) = logDelta1 + logDelta1ocean + (exp(logDelta2) + exp(logDelta2ocean)) * PredlnWA(i);
    //PredlnSMSY_S(i) = logDelta1(0) + exp(logDelta2(0)) * PredlnWA(i);
    //PredlnSMSY_O(i) = logDelta1(1) + exp(logDelta2(1)) * PredlnWA(i);
    //PredlnSMSY_S(i) = logDelta1 + exp(logDelta2(0)) * PredlnWA(i);
    //PredlnSMSY_O(i) = logDelta1 + exp(logDelta2(1)) * PredlnWA(i);
  }
  
  //ADREPORT(A_ar);
  //ADREPORT(A_std);
  //ADREPORT(logA_);
  //ADREPORT(logAadj_);
  //ADREPORT(SMSYadj);
  ADREPORT(SMSY_std);
  ADREPORT(SMSY_stream)
  //ADREPORT(SREP_stream)
  //ADREPORT(WA_stream)
  ADREPORT(SMSY_ocean*Scale_ocean)
  //ADREPORT(WA_ocean)
  //ADREPORT(logDelta1)
  ////ADREPORT(exp(logDelta2))
  //ADREPORT(Delta2_bounded)
  //ADREPORT(sigma_delta)
  //ADREPORT(slogDelta1)
  //ADREPORT(sDelta2_bounded)
  //ADREPORT(ssigma_delta)
  //ADREPORT(ologDelta1)
  //ADREPORT(ologDelta2)
  //ADREPORT(osigma_delta)
  //ADREPORT(gamma);
  //ADREPORT(LogR_Pred_ar);
  //ADREPORT(LogR_Pred_std);
  //ADREPORT(LogR_Pred_surv);
  ////Lierman pars
  //ADREPORT(logDelta1);
  //ADREPORT(logDelta1ocean);
  //ADREPORT(exp(logDelta2));
  //ADREPORT(exp(logDelta1ocean));
  ////Hierachical Pars
  //ADREPORT(logDelta1);
  
  ADREPORT(LogRS_Pred_std);
  ADREPORT(PredlnSMSYs_CI);
  ADREPORT(PredlnSMSYo_CI);
  //ADREPORT(PredlnSMSY_O);
  //ADREPORT(PredlnSMSY_S);
  REPORT(nLL_std);
  REPORT(ans);
  //ADREPORT(Sgen);
  return ans;
  
}
