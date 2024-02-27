
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
  Q = pow(tau, Type(2.0))  * (kappa4*M0 + Type(2.0)*kappa2*M1 + M2);
  return Q;
}
//------------------------------------------------------------------------------


// https://kaskr.github.io/adcomp/R__inla_8hpp_source.html#l00071

template<class Type>
SparseMatrix<Type> spde_Q_Aniso(spde_aniso_t<Type> spde,
                                Type kappa, Type tau,
                                matrix<Type> H){
  
  int i;
  Type kappa_pow2 = kappa*kappa;
  Type kappa_pow4 = kappa_pow2*kappa_pow2;
  
  int n_s = spde.n_s;
  int n_tri = spde.n_tri;
  vector<Type> Tri_Area = spde.Tri_Area;
  matrix<Type> E0 = spde.E0;
  matrix<Type> E1 = spde.E1;
  matrix<Type> E2 = spde.E2;
  matrix<int> TV = spde.TV;
  SparseMatrix<Type> G0 = spde.G0;
  SparseMatrix<Type> G0_inv = spde.G0_inv;
  SparseMatrix<Type> G1_aniso(n_s,n_s); 
  SparseMatrix<Type> G2_aniso(n_s,n_s);
  
  // Calculate adjugate of H
  matrix<Type> adj_H(2,2);
  adj_H(0,0) = H(1,1);
  adj_H(0,1) = -1 * H(0,1);
  adj_H(1,0) = -1 * H(1,0);
  adj_H(1,1) = H(0,0);
  
  
  // Calculate new SPDE matrices
  //
  // Calculate G1 - pt. 1
  array<Type> Gtmp(n_tri,3,3);
  for(i=0; i<n_tri; i++){    
    // 1st line: E0(i,) %*% adjH %*% t(E0(i,)), etc.    
    Gtmp(i,0,0) = (E0(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E0(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));  
    Gtmp(i,0,1) = (E1(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E1(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));  
    Gtmp(i,0,2) = (E2(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E2(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
    Gtmp(i,1,1) = (E1(i,0)*(E1(i,0)*adj_H(0,0)+E1(i,1)*adj_H(1,0)) + E1(i,1)*(E1(i,0)*adj_H(0,1)+E1(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
    Gtmp(i,1,2) = (E2(i,0)*(E1(i,0)*adj_H(0,0)+E1(i,1)*adj_H(1,0)) + E2(i,1)*(E1(i,0)*adj_H(0,1)+E1(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
    Gtmp(i,2,2) = (E2(i,0)*(E2(i,0)*adj_H(0,0)+E2(i,1)*adj_H(1,0)) + E2(i,1)*(E2(i,0)*adj_H(0,1)+E2(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
  }
  // Calculate G1 - pt. 2
  for(i=0; i<n_tri; i++){
    G1_aniso.coeffRef(TV(i,1),TV(i,0)) = G1_aniso.coeffRef(TV(i,1),TV(i,0)) + (Gtmp(i,0,1));  
    G1_aniso.coeffRef(TV(i,0),TV(i,1)) = G1_aniso.coeffRef(TV(i,0),TV(i,1)) + (Gtmp(i,0,1));  
    G1_aniso.coeffRef(TV(i,2),TV(i,1)) = G1_aniso.coeffRef(TV(i,2),TV(i,1)) + (Gtmp(i,1,2));  
    G1_aniso.coeffRef(TV(i,1),TV(i,2)) = G1_aniso.coeffRef(TV(i,1),TV(i,2)) + (Gtmp(i,1,2));  
    G1_aniso.coeffRef(TV(i,2),TV(i,0)) = G1_aniso.coeffRef(TV(i,2),TV(i,0)) + (Gtmp(i,0,2));  
    G1_aniso.coeffRef(TV(i,0),TV(i,2)) = G1_aniso.coeffRef(TV(i,0),TV(i,2)) + (Gtmp(i,0,2));  
    G1_aniso.coeffRef(TV(i,0),TV(i,0)) = G1_aniso.coeffRef(TV(i,0),TV(i,0)) + (Gtmp(i,0,0));  
    G1_aniso.coeffRef(TV(i,1),TV(i,1)) = G1_aniso.coeffRef(TV(i,1),TV(i,1)) + (Gtmp(i,1,1));  
    G1_aniso.coeffRef(TV(i,2),TV(i,2)) = G1_aniso.coeffRef(TV(i,2),TV(i,2)) + (Gtmp(i,2,2));  
  }
  G2_aniso = G1_aniso * G0_inv * G1_aniso; 
  
  return pow(tau, Type(2.0)) * (kappa_pow4*G0 + Type(2.0)*kappa_pow2*G1_aniso + G2_aniso);
}


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
  DATA_IVECTOR(spde_idx);    // SPDE index's
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
  Type tau = exp(logTau);
  Type kappa = exp(logKappa);
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
  
  // I might need to add a constraint which means the x sums to zero
  //
  
  // Get the contribution of the spatial field to the likelihood
  // Using the SPDE approximation
  Type nll3 = GMRF(Q)(x);
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
  Type ll4 = sum(dnorm(subsidence, okadaSubsidence, V, true));
  Type nll4 = -ll4;
  //----------------------------------------------------------------------------
  
  // Final likelihood is the product of p(theta) p(x|theta) p(y|x, theta).
  // Therefore the sum of log likelihoods
  Type nll = nll1 + nll2 + nll3 + nll4;
  //----------------------------------------------------------------------------
  
  // Calculate useful re-parameterisations
  //
  // Calculate the effective range parameter
  double nu = 1.0;                // nu = alpha-d/2 = 2-2/2 by eqn (2) in Lindgren 
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
