#include <TMB.hpp>

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
  DATA_VECTOR(SMSY);
  DATA_VECTOR(SREP);
  //DATA_IVECTOR(Stocks);
  DATA_IVECTOR(Inlets);
  //DATA_IVECTOR(CUs);
  //DATA_IVECTOR(MCTrials);
  //DATA_INTEGER(N_stks);
  DATA_INTEGER(N_inlets);
  //DATA_INTEGER(N_cus);

  
  PARAMETER_VECTOR(RicB);
  PARAMETER_VECTOR(logSgen);
  //PARAMETER(B_0);
  //PARAMETER(B_1);
  
  Type ans=0.0;
  int N_obs = SMSY.size(); 
  
  vector <Type> PredSREP(N_obs);
  Type sig = 1;
  
  // Estimate RicB using explicit formulation derived from Scheuerell 2014
  PredSREP = ( RicB * SMSY - log( 1 - RicB * SMSY ) ) / RicB;
  vector <Type> Diff1 = PredSREP - SREP;
  ans = -sum(dnorm(Diff1, 0, sig, true ));
  
  // Calculate RicA
  vector <Type> RicA(N_obs);
  RicA =  exp( RicB * SREP );
  
  // //Approximation for SMSY from Hilborn and Walters 1992
  // vector <Type> RicA = exp( (0.5 - SMSY / SREP ) / 0.07 ); 
  // vector <Type> RicB = log (RicA) / SREP;
  
  // Now estimate Sgen
  vector <Type> Sgen = exp(logSgen);
  vector <Type> PredlogSMSY(N_obs);
  Type Sgen_sig = 1;
  
  PredlogSMSY = log(RicA) + logSgen - RicB * Sgen;
  vector <Type> Diff2 = exp(PredlogSMSY) - SMSY;
  ans += -sum(dnorm(Diff2, 0, Sgen_sig, true ));
  
  vector <Type> Sgen_Inlet(N_inlets);//When I add MC trials, this should be N_inlets x N_MCtrials long
  Sgen_Inlet.setZero();
  //Wnen I add MC trials, I need to add another index Inlets_MC, which is N_inlets x N_MCtrials long, and repeats rep(0:6, N_MCtrials)
  
  // I need to scale Sgens so that I can add them. Then see if they are correct.
  
  //for (int mc = 0; mc < N_MCtrials; mc++){
      for (int i = 0; i < N_inlets; i++){
        for (int n = 0; n< N_obs; n++){
          //if(MCTrials(n)==mc){
            if(Inlets(n)==i) {//when I add MC trials, the i index should be changed to Inlets_MC(i)
              Sgen_Inlet(i) = Sgen_Inlet(i) + Sgen(n); 
            }
          //}
        }
      }
  //}
     
  ADREPORT(Sgen_Inlet)
  return ans;
  
}
