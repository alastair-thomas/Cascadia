
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  
  DATA_VECTOR(depth); // depths of centroids of subfaults
  DATA_VECTOR(subsidence); // subsidence data
  DATA_VECTOR(sigma); // standard deviation of errors
  DATA_MATRIX(okada);
  DATA_IVECTOR(spde_idx); // be careful that it is set up correctly, sparse vs dense
  DATA_STRUCT(spde, spde_t); // Mesh for basis functions, spde_t is the class
  
  // parameters with domain over entire real down
  PARAMETER_VECTOR(x);
  PARAMETER(log_lambda);
  //PARAMETER(log_dmax);
  PARAMETER(mu);
  PARAMETER(log_kappa);
  PARAMETER(log_tau);
  
  // why are these on a log scale when passed?
  Type tau = exp(log_tau);
  Type kappa = exp(log_kappa);
  Type lambda = exp(log_lambda);
  //Type dmax = exp(log_dmax);
  //Type lambda = 1.0/2000.0;
  //Type dmax = 30000.0;

  SparseMatrix<Type> Q = Q_spde(spde,kappa); // Q_spde is an inbuilt function
  
  // spde model for the spatial field
  // already negative and on the log scale
  Type nll1 = GMRF(Q)(x);

  
  vector<Type> untaperedSlips(spde_idx.size());
  untaperedSlips.setZero();
  
  // I just extract the parts of x which give slips
  // I don't understand the tau here
  for(int i=0; i<spde_idx.size(); i++){
    untaperedSlips(i) = exp(mu + (x(spde_idx(i)) / tau));
  }
  
  vector<Type> taper = 1.0 - exp(lambda*(depth-30000.0)); // 1.0 to make it a float
  vector<Type> taperedSlips = taper*untaperedSlips;
  //vector<Type> taperedSlips = untaperedSlips;
  
  // Then take the Okada model of the slips.
  
  // apply the okada matrix to the taperedSlips
  vector<Type> okadaSubsidence = okada * taperedSlips;
  
  // final nll is the product of p(x|theta) p(y|x, theta).
  // the MVNORM function assumes zero mean scale by mean works?
  // use n normal densities and sum up.
  
  Type ll2 = sum(dnorm(subsidence, okadaSubsidence, sigma, true));
  
  // calculate negative log likelihood
  Type nll = nll1 - ll2;
  
  // Then I need to sort out the tau
  // I don't know how to do this.
  
  //double nu = 1.0;            // nu = alpha-d/2 = 2-1 by eqn (2) in Lindgren 
  //Type rho = sqrt(8.0*nu)/kappa;  // Distance at which correlation has dropped to 0.1 (p.  4 in Lindgren)
  //ADREPORT(rho);
  // Add a report for the Q matrix so I can use it in the R code.
    
  return nll;
}
