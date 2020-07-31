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
  DATA_VECTOR(logRS);
  DATA_IVECTOR(stk);
  DATA_IVECTOR(yr);
  DATA_VECTOR(Scale);
  
  PARAMETER_VECTOR(logA);
  PARAMETER_VECTOR(logB);
  PARAMETER_VECTOR(logSigma);

  
  Type ans=0.0;
  int N_Obs = S.size(); 
  int N_stks = Scale.size();

  vector <Type> LogRS_Pred(N_Obs);
  vector <Type> sigma = exp(logSigma);
  vector <Type> nLL(N_Obs);
  
  // Standard Ricker model
  for (int i = 0; i<N_Obs; i++){
    LogRS_Pred(i) = logA(stk(i)) - exp(logB(stk(i))) * S(i);
    ans += -dnorm(LogRS_Pred(i), logRS(i),  sigma(stk(i)), true);    
    nLL(i) = -dnorm(LogRS_Pred(i), logRS(i),  sigma(stk(i)), true);
    
  }

  vector <Type> SMSY(N_stks);  
  vector <Type> SREP(N_stks);  
  
  for(int i=0; i<N_stks; i++){
    SMSY(i) =  Scale(i) * (1 - LambertW(exp(1-logA(i))) ) / exp(logB(i)) ;
  }
  SREP = Scale * logA / exp(logB);
  
  ADREPORT(logA);
  ADREPORT(SMSY);
  ADREPORT(SREP);
  ADREPORT(LogRS_Pred);
  REPORT(nLL);
  return ans;
  
}
