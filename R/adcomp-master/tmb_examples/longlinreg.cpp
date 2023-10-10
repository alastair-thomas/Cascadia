// Linear regression - 10^6 observations.
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Optional: demonstrate some TMBad specific optimizations
#ifdef TMBAD_FRAMEWORK
  using vectorize::vector;
#endif
  DATA_VECTOR(Y);
  DATA_VECTOR(x);
  PARAMETER(a);
  PARAMETER(b);
  PARAMETER(logSigma);
  Type sigmasq=exp(2*logSigma);
  ADREPORT(sigmasq);
  Type nll = -sum(dnorm(Y, a+b*x, exp(logSigma), true));
  return nll;
}
