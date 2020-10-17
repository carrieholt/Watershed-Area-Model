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
  DATA_SCALAR(logMuAs_mean);
  DATA_SCALAR(logMuAs_sig);
  DATA_SCALAR(logMuAo_mean);
  DATA_SCALAR(logMuAo_sig);
  DATA_SCALAR(HalfNormMean);
  DATA_SCALAR(HalfNormSig);
  DATA_SCALAR(HalfNormMeanA);
  DATA_SCALAR(HalfNormSigA);
  DATA_VECTOR(WA); 
  DATA_VECTOR(Scale);
  DATA_IVECTOR(Stream);
  DATA_INTEGER(SigRicPriorNorm);
  DATA_INTEGER(SigRicPriorGamma);
  DATA_INTEGER(SigRicPriorCauchy);
  DATA_INTEGER(SigDeltaPriorNorm);
  DATA_INTEGER(SigDeltaPriorGamma);
  DATA_INTEGER(SigDeltaPriorCauchy);
  DATA_SCALAR(Tau_dist);
  DATA_SCALAR(Tau_D_dist);
  //DATA_SCALAR(logDeltaSigma);
  //DATA_SCALAR(logNuSigma);
  
  DATA_SCALAR(sigDelta_mean);
  DATA_SCALAR(sigDelta_sig);
  DATA_SCALAR(sigNu_mean);
  DATA_SCALAR(sigNu_sig);
  
  //DATA_INTEGER(N_stream);
  //DATA_INTEGER(N_ocean);
  DATA_VECTOR(PredlnWA);
  DATA_VECTOR(TestlnWAo);
  DATA_VECTOR(TestlnWAs);
  

  PARAMETER_VECTOR(logA_std);
  PARAMETER_VECTOR(logB_std);
  PARAMETER_VECTOR(logSigma_std);
  //PARAMETER_VECTOR(logA_s);
  //PARAMETER_VECTOR(logB_s);
  //PARAMETER_VECTOR(logSigma_s);
  //PARAMETER_VECTOR(logA_o);
  //PARAMETER_VECTOR(logB_o);
  //PARAMETER_VECTOR(logSigma_o);
  PARAMETER(logMuAs);
  //PARAMETER(logSigmaAs);
  PARAMETER(logSigmaA);
  PARAMETER(logMuAo);
  //PARAMETER(logSigmaAo);
  PARAMETER(logDelta1);
  PARAMETER(logDelta1ocean);
  PARAMETER(logDelta2);
  PARAMETER(Delta2ocean);
  PARAMETER(logDeltaSigma);
  
  PARAMETER(logNu1);
  PARAMETER(logNu1ocean);
  PARAMETER(logNu2);
  PARAMETER(Nu2ocean);
  PARAMETER(logNuSigma);
  

  
  Type ans=0.0;
  int N_Obs = S_std.size(); 
  int N_stks_std = Scale.size();
  
  
  // Next: separate out S, logRS, stk (new 0:nstks by lh), yr, WA, Scale for stream and ocean type,and create new stream variable than is all 0 for stream and 1 for ocean
  // (the new WA, Scale, and Stream data series are needed so that they match up with the new order of stocks)
//  int j = 0;
//  int k = 0;
//  vector <Type> S_s(?);
//  vector <Type> logRS_s(?);
//  vector <int> stkNumber_s(?);
//  vector <int> stk_s(?);
//  vector <int> yr_s(?);
//  vector <Type> WA_s(N_stream);
//  vector <Type> Scale_s(N_stream);
//  vector <int> Stream_s(N_stream);
  
//  vector <Type> S_o(N_ocean);
//  vector <Type> logRS_o(N_ocean);
//  vector <int> stkNumber_o(N_ocean);
//  vector <int> stk_o(N_ocean);
//  vector <int> yr_o(N_ocean);
//  vector <Type> WA_o(N_ocean);
//  vector <Type> Scale_o(N_ocean);
//  vector <int> Stream_o(N_ocean);
//  int StreamNumber = 0;
//  int OceanNumber = 1;
//  
//  for(int ii=0; ii < N_Obs; ii++){
//    if(Stream[stk_std(ii)]==0){
//      S_s[j] = S_std[ii];
//      logRS_s[j] = logRS_std[ii];
//      stkNumber_s[j] = stk_std[ii];//original stock number
//      yr_s[j] = yr_std[ii];
      
 //     stk_s[j] = j;//stock index// This is wrong since j is sequenced by observation not stock
//      WA_s[j] = WA[stk_std(ii)];// This is wrong since this loops is over obs, not stocks (like WA[]), same for below
//      Scale_s[j] = Scale[stk_std(ii)];
//      Stream_s[j] = StreamNumber; 
      
//      j += 1;
//    }
//    if(Stream[stk_std(ii)]==1){
//      S_o[k] = S_std[ii];
//      logRS_o[k] = logRS_std[ii];
//      stkNumber_o[k] = stk_std[ii];
//      stk_o[k] = k;
//      yr_o[k] = yr_std[ii];
//      WA_o[k] = WA[ii];
//      Scale_o[k] = Scale[ii];
//      Stream_o[k] = OceanNumber; 
//      k += 1;
 //   }
//  }

//  int N_Obs_s = S_s.size();
//  int N_Obs_o = S_o.size();
  vector <Type> LogRS_Pred_std(N_Obs);
//  vector <Type> LogRS_Pred_s(N_Obs_s);
//  vector <Type> LogRS_Pred_o(N_Obs_o);
  vector <Type> sigma_std = exp(logSigma_std);
//  vector <Type> sigma_s = exp(logSigma_s);
//  vector <Type> sigma_o = exp(logSigma_o);
  Type sigmaA = exp(logSigmaA);
//  Type sigmaAs = exp(logSigmaAs);
//  Type sigmaAo = exp(logSigmaAo);
  vector <Type> nLL_std(N_Obs);
//  vector <Type> nLL_s(N_Obs_s);
//  vector <Type> nLL_o(N_Obs_o);
  
  // Standard Ricker model: stream
  
  //for (int i = 0; i<N_Obs_s; i++){
  for (int i = 0; i<N_Obs; i++){
    LogRS_Pred_std(i) = logA_std(stk_std(i)) - exp(logB_std(stk_std(i))) * S_std(i);
    //LogRS_Pred_s(i) = logA_s(stk_s(i)) - exp(logB_s(stk_s(i))) * S_s(i);
    ans += -dnorm(LogRS_Pred_std(i), logRS_std(i),  sigma_std(stk_std(i)), true);    
    //ans += -dnorm(LogRS_Pred_s(i), logRS_s(i),  sigma_s(stk_s(i)), true);    
    nLL_std(i) = -dnorm(LogRS_Pred_std(i), logRS_std(i),  sigma_std(stk_std(i)), true);
    //nLL_s(i) = -dnorm(LogRS_Pred_s(i), logRS_s(i),  sigma_s(stk_s(i)), true);
  }
  //// Standard Ricker model: ocean
  //for (int i = 0; i<N_Obs_o; i++){
  //  //for (int i = 0; i<N_Obs; i++){
  //  //LogRS_Pred_std(i) = logA_std(stk_std(i)) - exp(logB_std(stk_std(i))) * S_std(i);
  //  LogRS_Pred_o(i) = logA_o(stk_o(i)) - exp(logB_o(stk_o(i))) * S_o(i);
  //  //ans += -dnorm(LogRS_Pred_std(i), logRS_std(i),  sigma_std(stk_std(i)), true);    
  //  ans += -dnorm(LogRS_Pred_o(i), logRS_o(i),  sigma_o(stk_o(i)), true);    
  //  //nLL_std(i) = -dnorm(LogRS_Pred_std(i), logRS_std(i),  sigma_std(stk_std(i)), true);
  //  nLL_o(i) = -dnorm(LogRS_Pred_o(i), logRS_o(i),  sigma_o(stk_o(i)), true);
  //}
  
  
  // Add hierarchical structure to A ==============
  for(int i=0; i<N_stks_std; i++){
    // add prior on logA_std, 
    ans += -dnorm(logA_std(i), logMuAs + logMuAo * Stream(i), sigmaA, true );
     // add prior on sigma 
    if (SigRicPriorGamma == 1) {
       ans += -dgamma(pow(sigma_std(i),-2), Tau_dist, 1/Tau_dist, true);
    }
    if (SigRicPriorNorm == 1) {
      ans += -abs( dnorm( sigma_std(i), HalfNormMean, HalfNormSig, true) );
    }
    if (SigRicPriorCauchy == 1) {
    ans += - abs( dt( sigma_std(i), Type(1), true ));
    }
  }
  
  //// Add hierarchical structure to A: stream ==============
  //for(int i=0; i<N_stream; i++){
  //  // add prior on logA_std, 
  //  ans += -dnorm(logA_s(i), logMuAs, sigmaAs, true );
  //  // add prior on sigma 
  //  ans += -dgamma(pow(sigma_s(i),-2), Tau_dist, 1/Tau_dist, true);
  //}
  //// Add hierarchical structure to A: ocean ==============
  //for(int i=0; i<N_ocean; i++){
  //  // add prior on logA_std, 
  //  ans += -dnorm(logA_o(i), logMuAo, sigmaAo, true );
  //  // add prior on sigma 
  //  ans += -dgamma(pow(sigma_o(i),-2), Tau_dist, 1/Tau_dist, true);
  //}
  
    // Add priors for hyperpars ====================
  // MuA prior for stream type
  ans += -dnorm(logMuAs, logMuAs_mean, logMuAs_sig, true);
  // MuA prior for ocean type
  ans += -dnorm(logMuAo, logMuAo_mean, logMuAo_sig, true);
  // sigmaA prior
  if (SigRicPriorGamma == 1) {
    ans += -dgamma(pow(sigmaA,-2), Tau_dist, 1/Tau_dist, true);
  }
  if (SigRicPriorNorm == 1) {
    ans += -abs( dnorm( sigmaA, HalfNormMeanA, HalfNormSigA, true) );
  }
  if (SigRicPriorCauchy == 1) {
    ans += - abs(dt( sigmaA, Type(1), true));
  }
  
  
  //// First aggregate Ricker parameters into a single vector, and other co-variates in life-history order
  //vector <Type> logA_std(N_stks_std);
  //vector <Type> logB_std(N_stks_std);
  //vector <Type> WA_lhOrder(N_stks_std);
  //vector <Type> Scale_lhOrder(N_stks_std);
  //vector <int> stkNumber_lhOrder(N_stks_std);
  //vector <int> Stream_lhOrder(N_stks_std);

  //for(int i=0; i < N_stream; i++){
  //  logA_std[i] = logA_s[i];
  //  logB_std[i] = logB_s[i];
  //  WA_lhOrder[i] = WA_s[i];
  //  Scale_lhOrder[i] = Scale_s[i];
  //  stkNumber_lhOrder[i] = stkNumber_s[i];
  //  Stream_lhOrder[i] = Stream_s[i];
  //}
  //for(int i=0; i < N_ocean; i++){
  //  logA_std[N_stream + i] = logA_o[i];
  //  logB_std[N_stream + i] = logB_o[i];
  //  WA_lhOrder[N_stream + i] = WA_o[i];
  //  Scale_lhOrder[N_stream + i] = Scale_o[i];
  //  stkNumber_lhOrder[N_stream + i] = stkNumber_o[i];
  //  Stream_lhOrder[N_stream + i] = Stream_o[i];
  //}

  
  //Calculate SMSY and SREP
  vector <Type> SMSY_std(N_stks_std);  
  vector <Type> SREP_std(N_stks_std);  
  
  //For SMSY calculation, 
  for(int i=0; i<N_stks_std; i++){
    SMSY_std(i) =  (1 - LambertW( exp (1- logA_std(i)) ) ) / exp(logB_std(i)) ;
  }
  SREP_std = logA_std / exp(logB_std);
  
  
  //Liermann's model with both stream and ocean type=================
  vector <Type> PredlnSMSY(N_stks_std);
  vector <Type> PredlnSREP(N_stks_std);
  Type sigma_delta = exp(logDeltaSigma);
  Type sigma_nu = exp(logNuSigma);
  
  for (int i=0; i<N_stks_std; i++){
    PredlnSMSY(i) = logDelta1 + logDelta1ocean * Stream(i) + ( exp(logDelta2) + Delta2ocean * Stream(i) ) * log(WA(i)) ;
    ans += -dnorm( PredlnSMSY(i), log(SMSY_std(i) * Scale(i) ),  sigma_delta, true);
    PredlnSREP(i) = logNu1 + logNu1ocean * Stream(i) + ( exp(logNu2) + Nu2ocean * Stream(i) ) * log(WA(i)) ;
    ans += -dnorm( PredlnSREP(i), log(SREP_std(i) * Scale(i) ),  sigma_nu, true);
    //PredlnSMSY(i) = logDelta1 + logDelta1ocean * Stream(i) + ( exp(logDelta2) + Delta2ocean * Stream(i) ) * log(WA(i)) ;
    //ans += -dnorm( PredlnSMSY(i), log(SMSY_std(i) * Scale(i) ),  sigma_delta, true);
  }
  
  // Normal prior on sigma_delta and sigma_nu
  if (SigDeltaPriorNorm == 1) {
    ans += -dnorm(sigma_delta, sigDelta_mean, sigDelta_sig, true);
    ans += -dnorm(sigma_nu, sigNu_mean, sigNu_sig, true);
  }
  
  // Inverse gamma prior on sigma_delta and sigma_nu
  if (SigDeltaPriorGamma == 1) {
    ans += -dgamma(pow(sigma_delta,-2), Tau_D_dist, 1/Tau_D_dist, true);
    ans += -dgamma(pow(sigma_nu,-2), Tau_D_dist, 1/Tau_D_dist, true);
  }
  
  // Half cauchy prior on sigma_delta and sigma_nu
  if (SigDeltaPriorCauchy == 1) {
    ans += -abs( dt( sigma_delta, Type(1), true));
    ans += - abs( dt( sigma_nu, Type(1), true ));
  }
  
    // Get predicted values for plotting  WA regresssion with CIs
  int N_pred = PredlnWA.size();
  vector <Type> PredlnSMSYs_CI(N_pred);
  vector <Type> PredlnSMSYo_CI(N_pred);
  vector <Type> PredlnSREPs_CI(N_pred);
  vector <Type> PredlnSREPo_CI(N_pred);
  
  for (int i=0; i<N_pred; i++){
    PredlnSMSYs_CI(i) = logDelta1 + exp(logDelta2) * PredlnWA(i);
    PredlnSMSYo_CI(i) = logDelta1 + logDelta1ocean + (exp(logDelta2) + Delta2ocean) * PredlnWA(i);
    PredlnSREPs_CI(i) = logNu1 + exp(logNu2) * PredlnWA(i);
    PredlnSREPo_CI(i) = logNu1 + logNu1ocean + (exp(logNu2) + Nu2ocean) * PredlnWA(i);
  }
  
  //// Get predicted values for stream-type Test stocks with CIs
  int N_tests = TestlnWAs.size();
  vector <Type> TestlnSMSYs(N_tests);
  vector <Type> TestlnSREPs(N_tests);
  
  for (int i=0; i<N_tests; i++){
    TestlnSMSYs(i) = logDelta1 + exp(logDelta2) * TestlnWAs(i);
    TestlnSREPs(i) = logNu1 + exp(logNu2) * TestlnWAs(i);
  }
  
  
  ///Get predicted values for ocean-type Test stocks with CIs
  int N_testo = TestlnWAo.size();
  vector <Type> TestlnSMSYo(N_testo);
  vector <Type> TestlnSREPo(N_testo);
  
  for (int i=0; i<N_testo; i++){
    TestlnSMSYo(i) = logDelta1 + + logDelta1ocean + (exp(logDelta2) + Delta2ocean) * TestlnWAo(i);
    TestlnSREPo(i) = logNu1 + logNu1ocean + (exp(logNu2) + Nu2ocean) * TestlnWAo(i);
  }
  
  vector <Type> lnSMSY = log(SMSY_std*Scale);
  vector <Type> lnSREP = log(SREP_std*Scale);
  
  
  ADREPORT(SMSY_std);
  ADREPORT(SREP_std);
  ADREPORT(LogRS_Pred_std);
  ADREPORT(PredlnSMSY);
  ADREPORT(PredlnSREP);
  ADREPORT(lnSMSY);
  ADREPORT(lnSREP);
  ADREPORT(PredlnSMSYs_CI);
  ADREPORT(PredlnSMSYo_CI);
  ADREPORT(PredlnSREPs_CI);
  ADREPORT(PredlnSREPo_CI);
  //REPORT(stkNumber_lhOrder);//This is the order of stocks in the output (different from order in input!)
  ADREPORT(TestlnSMSYs);
  ADREPORT(TestlnSREPs);
  ADREPORT(TestlnSMSYo);
  ADREPORT(TestlnSREPo);
  REPORT(nLL_std);
  //REPORT(nLL_o);
  return ans;
  
}
