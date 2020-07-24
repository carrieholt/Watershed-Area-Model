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
  DATA_VECTOR(S);
  //DATA_VECTOR(logR);
  DATA_VECTOR(logRS);
  DATA_IVECTOR(stk);
  DATA_IVECTOR(yr);
  DATA_IVECTOR(model);
  //DATA_SCALAR(Sgen_sig);
  

  //PARAMETER_VECTOR(logA_ar);
  //PARAMETER_VECTOR(logB_ar);
  //PARAMETER_VECTOR(rho);
  //PARAMETER_VECTOR(logSigma_ar);
  PARAMETER_VECTOR(logA_std);
  PARAMETER_VECTOR(logB_std);
  PARAMETER_VECTOR(logSigma_std);
  //PARAMETER_VECTOR(logSgen);
  
  
  Type ans=0.0;
  int N_Obs = S.size(); 
  //vector <Type> LogR_Pred_ar(N_Obs);
  //vector <Type> sigma_ar = exp(logSigma_ar);
  //vector <Type> A_ar = exp(logA_ar);
  //vector <Type> err(N_Obs);

  
  // Ricker AR1 likelihood
  //for (int i = 0; i<N_Obs; i++){
  //  if (yr(i) == 0) {LogR_Pred_ar(i) = logA_ar(stk(i)) + log(S(i)) - exp(logB_ar(stk(i))) * S(i);
  //      err(i) = logR(i) - LogR_Pred_ar(i);
        
  //    } else if (yr(i) >= 1) { LogR_Pred_ar(i) = logA_ar(stk(i)) + log(S(i)) - exp(logB_ar(stk(i))) * S(i) + rho(stk(i))*err(i-1);
  //      err(i) = logR(i) - LogR_Pred_ar(i);
   //     }
    
  //  ans += -dnorm(LogR_Pred_ar(i), logR(i),  sigma_ar(stk(i)), true);
 // }
  
  //vector <Type> LogR_Pred_std(N_Obs);
  vector <Type> LogRS_Pred_std(N_Obs);
  vector <Type> sigma_std = exp(logSigma_std);
  vector <Type> A_std = exp(logA_std);
  
  // Standard Ricker model
  for (int i = 0; i<N_Obs; i++){
    
    //LogR_Pred_std(i) = logA_std(stk(i)) + log(S(i)) - exp(logB_std(stk(i))) * S(i);
    LogRS_Pred_std(i) = logA_std(stk(i)) - exp(logB_std(stk(i))) * S(i);
    
    //ans += -dnorm(LogR_Pred_std(i), logR(i),  sigma_std(stk(i)), true);
    ans += -dnorm(LogRS_Pred_std(i), logRS(i),  sigma_std(stk(i)), true);    
  }

  //Compile parameters
 // int N_stks = model.size(); 
 // vector <Type> logA_(N_stks);
//  vector <Type> logAadj_(N_stks);
  //vector <Type> logB_(N_stks);
  //vector <Type> logSigma_(N_stks);
  //vector <Type> SMSYadj(N_stks);  
 // vector <Type> SMSY(N_stks);  
  
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
  int N_stks = model.size(); 
  vector <Type> logAadj_(N_stks);
  vector <Type> SMSYadj(N_stks);  
  vector <Type> SMSY(N_stks);  
  vector <Type> B = exp(logB_std);
  vector <Type> SREPadj(N_stks);
  
  for(int i=0; i<N_stks; i++){
    logAadj_(i) = logA_std(i) + pow( exp( logSigma_std(i) ), 2) /2;
    SMSYadj(i) =  (1 - LambertW(exp(1-logAadj_(i))) ) / B(i) ;
    SMSY(i) =  (1 - LambertW(exp(1-logA_std(i))) ) / B(i) ;
  }
  SREPadj = logAadj_ / B;
  
  //ADREPORT(A_ar);
  ADREPORT(A_std);
  //ADREPORT(logA_);
  ADREPORT(logAadj_);
  ADREPORT(SMSYadj);
  ADREPORT(SMSY);
  ADREPORT(SREPadj);
  
  //ADREPORT(LogR_Pred_ar);
  //ADREPORT(LogR_Pred_std);
  ADREPORT(LogRS_Pred_std);
  //ADREPORT(Sgen);
  return ans;
  
}
