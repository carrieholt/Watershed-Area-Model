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
  DATA_VECTOR(SMSY);
  DATA_VECTOR(WA);
  DATA_VECTOR(Scale);
  DATA_SCALAR(Tau_dist);

  PARAMETER(logDelta1);
  PARAMETER(logDelta2);
  PARAMETER(logDeltaSigma);

  
  Type ans=0.0;
  int N_stks = SMSY.size();
  vector <Type> PredlnSMSY(N_stks);
  //Type Delta2_bounded = invlogit(Delta2);
  Type sigma_delta = exp(logDeltaSigma);

  for (int i=0; i<N_stks; i++){
    PredlnSMSY(i) = logDelta1 + exp(logDelta2) * log(WA(i));
    ans += -dnorm(PredlnSMSY(i), log(SMSY(i)*Scale(i)),  sigma_delta, true);
  }
  // Add Inverse gamma prior on sigma_delta^2
  ans += -dgamma(pow(sigma_delta,-2), Tau_dist, 1/Tau_dist, true);

  ADREPORT(logDelta1);
  ADREPORT(logDelta2);
  ADREPORT(sigma_delta);
  return ans;
  
}
