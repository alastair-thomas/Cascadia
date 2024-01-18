
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  
  DATA_VECTOR(depth); // depths of centroids of subfaults
  DATA_VECTOR(subsidence); // subsidence data
  DATA_VECTOR(V); // standard deviation of errors
  DATA_MATRIX(okada); // Okada matrix
  DATA_IVECTOR(spde_idx); // spde index's. Saved as an index vector
  DATA_STRUCT(spde, spde_t); // Mesh for basis functions, spde_t is the class
  
  // parameters with domain over entire real down
  PARAMETER_VECTOR(x); // random effect spatial field
  PARAMETER(log_lambda);
  PARAMETER(mu);
  PARAMETER(log_kappa);
  PARAMETER(log_tau);
  
  // why are these on a log scale when passed?
  Type tau = exp(log_tau);
  Type kappa = exp(log_kappa);
  Type lambda = exp(log_lambda);

  // the precision matrix should be multiplied by tau?
  // following equation 3.6 from my thesis:
  // x | kappa, tau ~ MVN(0, 1/tau * Q^(-1))
  SparseMatrix<Type> Q = tau * Q_spde(spde,kappa); // Q_spde is an inbuilt function
  
  // spde model for the spatial field
  // already negative and on the log scale
  Type nll1 = GMRF(Q)(x);
  
  vector<Type> untaperedSlips(spde_idx.size());
  untaperedSlips.setZero();
  
  // Here I extract the parts of x which give slips
  for(int i=0; i<spde_idx.size(); i++){
    untaperedSlips(i) = exp(mu + x(spde_idx(i)));
  }
  
  vector<Type> taper = exp(-lambda*depth);
  vector<Type> taperedSlips = taper*untaperedSlips;
  
  // apply the okada matrix to the taperedSlips
  vector<Type> okadaSubsidence = okada * taperedSlips;
  
  // use n normal densities and sum up.
  Type ll2 = sum(dnorm(subsidence, okadaSubsidence, V, true));
  
  // final nll is the product of p(x|theta) p(y|x, theta).
  // calculate negative log likelihood
  Type nll = nll1 - ll2;
  
  
  double nu = 1.0;            // nu = alpha-d/2 = 2-1 by eqn (2) in Lindgren 
  Type rho = sqrt(8.0*nu)/kappa;  // Distance at which correlation has dropped to 0.1 (p.  4 in Lindgren)
  
  Type a = 4.0 + kappa*kappa;
  Type sigmaSquared = 1.0 / (4.0*3.14*nu*(a - 4.0));
    
  // report values needed for more analysis
  REPORT(Q);
  // For testing purposes
  // report values for testing
  REPORT(nll1);            // the negative log likehood from the SPDE approximation
  REPORT(untaperedSlips);  // untapered slips = exp(mu + Ax)
  REPORT(taperedSlips);    // tapered slips   = t * untaperedslips
  REPORT(okadaSubsidence); // okada subsidence = Gs
  REPORT(ll2);             // The log likelihood from the sum of normals
  REPORT(nll);             // nll1 - ll2
  REPORT(x);               // the unscaled distribution.
  
  ADREPORT(rho);
  ADREPORT(sigmaSquared);
  
  return nll;
  
}
