
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
// AT - Modified to use tau and kappa
template<class Type>
SparseMatrix<Type> spde_Q(Type kappa, Type tau, SparseMatrix<Type> M0,
                          SparseMatrix<Type> M1, SparseMatrix<Type> M2) {
  SparseMatrix<Type> Q;
  Type kappa2 = pow(kappa, 2.0);
  Type kappa4 = pow(kappa, 4.0);
  Q = pow(tau, 2.)  * (kappa4*M0 + Type(2.0)*kappa2*M1 + M2);
  return Q;
}
//------------------------------------------------------------------------------


// helper function to use the same penalized complexity prior on
// matern params that is used in INLA
// Cite: https://github.com/umut-altay/GeoAdjust/blob/main/simulations.cpp
// AT - Modified to use tau and kappa
//    - Also now returns the negative log likelihood
template<class Type>
Type dPCPriSPDE(Type tau, Type kappa, vector<Type> maternPri)
{
  // The PC prior in matern gives:
  // P(range < a) = b
  // P(sigma > c) = d
  // First extract the components
  Type a = maternPri[0];
  Type b = maternPri[1];
  Type c = maternPri[2];
  Type d = maternPri[3];
  
  Type penalty; // prior contribution to nll
  
  double alpha = 2.0; // alpha
  double dimension = 2.0;  // dimension
  double nu = 1.0; // nu = alpha - dimension/2 by Lindgren
  Type lambda1 = -log(b) * pow(a, dimension/2.0);
  Type lambda2 = -log(d) / c;
  Type range   = sqrt(8.0*nu) / kappa;
  Type sigma   = 1.0 / sqrt(4.0 * 3.14159265359 * pow(tau, 2.0) * pow(kappa, 2.0));
  
  penalty = -alpha * log(range) - lambda1 * pow(range, -dimension/2.0) - lambda2 * sigma;
  // Note: (rho, sigma) --> (x=log kappa, y=log tau) -->
  // transforms: rho = sqrt(8)/e^x & sigma = 1/(sqrt(4pi)*e^x*e^y)
  // --> Jacobian: |J| propto e^(-y -2x)
  Type jacobian = -log(tau) - 2.0*log(kappa);
  penalty += jacobian;
  
  return -penalty;
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
  DATA_SPARSE_MATRIX(A);     // SPDE Projection Matrix
  DATA_SPARSE_MATRIX(M0);    // Matrix to calculate GMRF precision matrix
  DATA_SPARSE_MATRIX(M1);    // `` ''
  DATA_SPARSE_MATRIX(M2);    // `` ''
  DATA_VECTOR(maternPri);    // The hyper-parameters for the PC prior on matern
  DATA_VECTOR(taperPri);     // The shape and scale for gamma prior on the taper
  //----------------------------------------------------------------------------
  
  // Parameters with domain over entire real domain
  //
  PARAMETER_VECTOR(x);  // Random effect spatial field
  PARAMETER(logLambda); // Controls the taper effect
  PARAMETER(mu);        // The mean untapered slip
  PARAMETER(logKappa);  // The scale parameter
  PARAMETER(logTau);    // The variance parameter
  //----------------------------------------------------------------------------
  
  // Transform parameters
  Type tau    = exp(logTau);
  Type kappa  = exp(logKappa);
  Type lambda = exp(logLambda);
  //----------------------------------------------------------------------------
  
  // Calculate the contribution of the prior on the spatial effects
  //
  Type nll1 = dPCPriSPDE(tau, kappa, maternPri);
  //----------------------------------------------------------------------------
  
  // Calculate the contribution of the prior on for lambda
  //
  Type shape = taperPri[0]; // extract the shape hyper-parameter
  Type scale = taperPri[1]; // extract the scale hyper-parameter
  Type ll2 = dgamma(lambda, shape, scale, true); // get contribution
  ll2 += logLambda; // Add in jacobian factor
  Type nll2 = -ll2; // negate
  //----------------------------------------------------------------------------
  
  // Calculate the contribution of the spatial effect via SPDE approximation:
  //
  // Build the sparse precision matrix requires for the GMRF approximation
  SparseMatrix<Type> Q = spde_Q(kappa, tau, M0, M1, M2);
  
  // Get the contribution of the spatial field to the likelihood
  // Using the SPDE approximation
  Type nll3 = GMRF(Q)(x);
  //----------------------------------------------------------------------------
  
  // Now calculate the non spatial contribution form the data:
  //
  // Create the untapered slips
  vector<Type> untaperedSlips = exp(mu)*exp(A*x);
  
  // Create the taper based on depth and the lambda parameter
  vector<Type> taper = exp(-lambda*depth);
  
  // Taper the slips
  vector<Type> taperedSlips = taper*untaperedSlips;
  
  // Apply the okada matrix to the taperedSlips to get subsidence
  vector<Type> okadaSubsidence = okada * taperedSlips;
  
  // Subsidence is normally distributed
  // Calculate the contribution of each one on the log scale
  Type ll4 = sum(dnorm(subsidence, okadaSubsidence, V, true));
  Type nll4 = -ll4;
  //----------------------------------------------------------------------------
  
  // Final likelihood is the product of p(theta) p(x|theta) p(y|x, theta).
  // Therefore the sum of log likelihoods
  Type nll = nll1 + nll2 + nll3 + nll4;
  //----------------------------------------------------------------------------
  
  // // Calculate useful re-parameterisations
  // //
  // // Calculate the effective range parameter
  // double nu = 1.0;                // nu = alpha-d/2 = 2-2/2 by eqn (2) in Lindgren 
  // Type rho = sqrt(8.0*nu)/kappa;  // Distance where correlation is ~0.1
  // 
  // // Calculate the marginal variance
  // Type sigma2 = 1.0 / (4.0*3.14159265359*pow(kappa, 2.)*pow(tau, 2.));
  // //----------------------------------------------------------------------------
  // 
  // // Report values for testing
  // REPORT(nll1);            // -log(p(x|theta))
  // REPORT(untaperedSlips);  // Untapered slips  = exp(mu + Ax)
  // REPORT(taperedSlips);    // Tapered slips    = taper*untaperedslips
  // REPORT(okadaSubsidence); // Okada subsidence = Gs
  // REPORT(ll2);             // -log(p(y|x, theta))
  // REPORT(nll);             // Total negative log likelihood
  // REPORT(x);               // Spatial random effects
  // //----------------------------------------------------------------------------
  // 
  // // Report the important parameters so that SDs can be found
  // ADREPORT(rho);
  // ADREPORT(sigma2);
  // //----------------------------------------------------------------------------
  // 
  // Return the total negative log likelihood
  return nll;
  
}
