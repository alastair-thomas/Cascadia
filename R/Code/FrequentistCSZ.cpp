
#include <TMB.hpp>
using namespace R_inla;
using namespace density;
using namespace Eigen;


// helper function to make sparse SPDE precision matrix
// Inputs:
//    logkappa: log(kappa) parameter value
//    logtau: log(tau) parameter value
//    M0, M1, M2: these sparse matrices are output from:
//     R::INLA::inla.spde2.matern()$param.inla$M*
// Cite: https://github.com/umut-altay/GeoAdjust/blob/main/simulations.cpp
// Modified to use tau and kappa
template<class Type>
SparseMatrix<Type> spde_Q(Type kappa, Type tau, SparseMatrix<Type> M0,
                          SparseMatrix<Type> M1, SparseMatrix<Type> M2) {
  SparseMatrix<Type> Q;
  Type kappa2 = pow(kappa, 2.);
  Type kappa4 = pow(kappa, 4.);
  Q = pow(tau, 2.)  * (kappa4*M0 + Type(2.0)*kappa2*M1 + M2);
  return Q;
}
//------------------------------------------------------------------------------

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // The data passed from R
  //
  DATA_VECTOR(depth);        // Depths of centroids of subfaults
  DATA_VECTOR(subsidence);   // Subsidence data
  DATA_VECTOR(V);            // Standard deviation of errors
  DATA_MATRIX(okada);        // Okada matrix
  DATA_IVECTOR(spde_idx);    // SPDE index's
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);
  //DATA_STRUCT(spde, spde_t); // Matrices to calculate the GMRF representation
  //----------------------------------------------------------------------------
  
  // Parameters with domain over entire real domain
  //
  PARAMETER_VECTOR(x);   // Random effect spatial field
  PARAMETER(logLambda); // Controls the taper effect
  PARAMETER(mu);         // The mean untapered slip
  PARAMETER(logKappa);  // The scale parameter
  PARAMETER(logTau);   // The variance parameter
  //----------------------------------------------------------------------------
  
  // Transform parameters
  Type tau = exp(logTau);
  Type kappa = exp(logKappa);
  Type lambda = exp(logLambda);
  //----------------------------------------------------------------------------

  // Calculate the contribution of the spatial effect via SPDE approximation:
  //
  // Build the sparse precision matrix requires for the GMRF approximation
  // The function Q_spde builds the precision matrix with unit tau
  // The formula for alpha=2 is:
  // Q = tau^2 (kappa^4 C + 2 kappa^2 G + G C^-1 G) from 
  //SparseMatrix<Type> Q = (tau2)*Q_spde(spde,kappa);
  SparseMatrix<Type> Q = spde_Q(kappa, tau, M0, M1, M2);
  
  // Get the contribution of the spatial field to the likelihood
  // Using the SPDE approximation
  Type nll1 = GMRF(Q)(x);
  //----------------------------------------------------------------------------
  
  // Now calculate the non spatial contribution form the data:
  //
  // Create the untapered slips
  vector<Type> untaperedSlips(spde_idx.size());
  untaperedSlips.setZero();

  // Calculate the untapered slips based on the spatial field
  for(int i=0; i<spde_idx.size(); i++){
    untaperedSlips(i) = exp(mu)*exp(x(spde_idx(i)));
  }
  
  // Create the taper based on depth and the lambda parameter
  vector<Type> taper = exp(-lambda*depth);
  
  // Taper the slips
  vector<Type> taperedSlips = taper*untaperedSlips;
  
  // Apply the okada matrix to the taperedSlips to get subsidence
  vector<Type> okadaSubsidence = okada * taperedSlips;
  
  // Subsidence is normally distributed
  // Calculate the contribution of each one on the log scale
  Type ll2 = sum(dnorm(subsidence, okadaSubsidence, V, true));
  Type nll2 = -ll2;
  //----------------------------------------------------------------------------
  
  // Final nll is the product of p(x|theta) p(y|x, theta).
  // Therefore the sum of log likelihoods
  Type nll = nll1 + nll2;
  //----------------------------------------------------------------------------
  
  // Calculate useful re-parameterisations
  //
  // Calculate an understandable range parameter
  double nu = 1.0;            // nu = alpha-d/2 = 2-2/2 by eqn (2) in Lindgren 
  Type rho = sqrt(8.0*nu)/kappa;  // Distance where correlation is ~0.1
  
  // Calculate the marginal variance
  Type sigma2 = 1.0 / (4.0*3.14159265359*pow(kappa, 2.)*pow(tau, 2.));
  //----------------------------------------------------------------------------
  
  // Report values for testing
  REPORT(nll1);            // -log(p(x|theta))
  REPORT(untaperedSlips);  // Untapered slips  = exp(mu + Ax)
  REPORT(taperedSlips);    // Tapered slips    = taper*untaperedslips
  REPORT(okadaSubsidence); // Okada subsidence = Gs
  REPORT(ll2);             // -log(p(y|x, theta))
  REPORT(nll);             // Total negative log likelihood
  REPORT(x);               // Spatial random effects
  //----------------------------------------------------------------------------
  
  // Report the important parameters so that SDs can be found
  ADREPORT(rho);
  ADREPORT(sigma2);
  //----------------------------------------------------------------------------
  
  // Return the total negative log likelihood
  return nll;
  
}
