// I need to implement the negative log likelihood here

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  
  DATA_VECTOR(depth); // depths of centroids of subfaults
  DATA_VECTOR(subsidence);
  DATA_MATRIX(sigma); // diagonal covariance matrix
  DATA_MATRIX(okada);
  DATA_MATRIX(spde_proj); // be careful that it is set up correctly, sparse vs dense
  
  DATA_STRUCT(spde, spde_t); // Mesh for basis functions, spde_t is the class
  
  // parameters with domain over entire real down
  PARAMETER_VECTOR(x);
  PARAMETER(log_lambda);
  PARAMETER(log_dMax);
  PARAMETER(mu);
  PARAMETER(log_kappa);
  PARAMETER(log_tau);
  
  Type tau = exp(log_tau);
  Type kappa = exp(log_kappa);

  Type nll = 0.0;

  SparseMatrix<Type> Q = Q_spde(spde,kappa); // Q_spde is an inbuilt function
  
  nll = GMRF(Q)(x);                              // Negative log likelihood

  // scale by tau afterwards
  
  vector<Type> untaperedSlips = exp(mu + x / sqrt(tau));
  
  vector<Type> taper = 1.0 - exp(-lambda*depth); // 1.0 to make it a float
  vector<Type> taperedSlips = taper*untaperedSlips;
  // Then take the Okada model of the slips.
  
  // apply the okada matrix to the taperedSlips
  vector<Type> okadaSlips = okada %*% taperedSlips;
  
  // final nll is the product of p(x|theta) p(y|x, theta).
  // the MVNORM function assumes zero mean scale by mean works?
  nll = nll * (MVNORM(sigma)(subsidence) + okadaSlips);
  
  // Then I need to sort out the tau
  // I don't know how to do this.
  
  double nu = 1.0;            // nu = alpha-d/2 = 2-1 by eqn (2) in Lindgren 
  Type rho = sqrt(8.0*nu)/kappa;  // Distance at which correlation has dropped to 0.1 (p.  4 in Lindgren)
  ADREPORT(rho);
    
  return nll;
}
