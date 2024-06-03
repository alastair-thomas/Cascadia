
#include <TMB.hpp>
using namespace R_inla;
using namespace density;
using namespace Eigen;
using namespace tmbutils;

// List of matrices
// Cite: https://github.com/kaskr/adcomp/issues/96
// AT - Changed from sparse matrix to matrix
//    - Also doesn't check anything, it assumed inputs are correct
//    - Used for the Okada matrices
template<class Type>
struct LOM_t : vector<matrix<Type> > {
  LOM_t(SEXP x){  // x is the list passed from R
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP m = VECTOR_ELT(x, i);
      (*this)(i) = asMatrix<Type>(m);
    }
  }
};

// List of vectors
// Cite: https://github.com/kaskr/adcomp/issues/96
// AT - Changed from sparse matrix to vectors
//    - Also doesn't check anything, it assumed inputs are correct
//    - Used for subsidence data
template<class Type>
struct LOV_t : vector<vector<Type> > {
  LOV_t(SEXP x){  // x is the list passed from R
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP v = VECTOR_ELT(x, i);
      (*this)(i) = asVector<Type>(v);
    }
  }
};

// Helper function to make sparse SPDE precision matrix
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
  Q = pow(tau, 2.0)  * (kappa4*M0 + Type(2.0)*kappa2*M1 + M2);
  return Q;
}
//------------------------------------------------------------------------------


// Helper function to use the same penalized complexity prior on
// Matern params that are used in INLA
// Cite: https://github.com/umut-altay/GeoAdjust/blob/main/simulations.cpp
// AT - Modified to use tau and kappa
//    - Also now returns the negative log likelihood
template<class Type>
Type dPCPriSPDE(Type tau, Type kappa,
                vector<Type> maternPri)
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
  //  transforms: rho = sqrt(8)/e^x & sigma = 1/(sqrt(4pi)*e^x*e^y)
  //  --> Jacobian: |J| propto e^(-y -2x)
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
  DATA_VECTOR(depth);             // Depths of centroids of subfaults
  DATA_STRUCT(subsidence, LOV_t); // List of subsidence data vectors
  DATA_STRUCT(v, LOV_t);          // List of subsidence error vectors
  DATA_STRUCT(okada, LOM_t);      // List of Okada matrices
  DATA_SPARSE_MATRIX(A);          // SPDE projection matrix
  DATA_SPARSE_MATRIX(M0);         // Matrix to calculate GMRF precision matrix
  DATA_SPARSE_MATRIX(M1);         // `` ''
  DATA_SPARSE_MATRIX(M2);         // `` ''
  DATA_VECTOR(maternPri);         // PC prior on matern
  DATA_VECTOR(taperPri);          // Prior on the taper
  //----------------------------------------------------------------------------
  
  // Parameters with domain over entire real domain
  //
  PARAMETER_MATRIX(X);        // Random effect spatial field
  PARAMETER_VECTOR(w);        // The shared spatial effect
  PARAMETER(logLambda);       // Controls the taper effect
  PARAMETER(mu);              // The mean untapered slip
  PARAMETER(logKappaX);       // The individual scale parameter
  PARAMETER(logTauX);         // The individual variance parameter
  PARAMETER(logKappaW);       // The shared scale parameter
  PARAMETER(logTauW);         // The shared variance parameter
  //----------------------------------------------------------------------------
  
  // Transform parameters
  //
  Type tauX   = exp(logTauX);
  Type kappaX = exp(logKappaX);
  Type tauW   = exp(logTauW);
  Type kappaW = exp(logKappaW);
  Type lambda = exp(logLambda);
  
  //----------------------------------------------------------------------------
  
  // Calculate the contribution of the prior on the individual spatial effects
  //
  Type nll1 = dPCPriSPDE(tauX, kappaX, maternPri);
  //----------------------------------------------------------------------------
  
  // Calculate the contribution of the prior on the shared spatial effects
  //
  Type nll2 = dPCPriSPDE(tauW, kappaW, maternPri);
  //----------------------------------------------------------------------------
  
  // Calculate the contribution of the prior on for lambda
  //
  Type shape = taperPri[0]; // extract the shape hyper-parameter
  Type scale = taperPri[1]; // extract the scale hyper-parameter
  Type ll3 = dgamma(lambda, shape, scale, true); // get contribution
  Type jacobian = logLambda; // Jacobian factor
  Type nll3 = -ll3 - jacobian; // negate
  //----------------------------------------------------------------------------
  
  // Calculate the contribution of the spatial effects via SPDE approximation:
  // This is done for each earthquake
  // And then for the shared spatial field

  // Build the sparse precision matrix requires for the GMRF approximation
  SparseMatrix<Type> Qx = spde_Q(kappaX, tauX, M0, M1, M2);
  
  Type M = X.cols();   // The number of earthquakes
  Type nll4 = 0.0;     // Initialise the spatial negative log likelihood
  
  // Loop over each earthquake
  for(int j=0; j<M; j++){
    nll4 += GMRF(Qx)(X.col(j)); // contribution of jth field
  }
  
  // Add the contribution for the shared spatial effect
  SparseMatrix<Type> Qw = spde_Q(kappaW, tauW, M0, M1, M2);
  nll4 += GMRF(Qw)(w);
  //----------------------------------------------------------------------------
  
  // Calculate the non spatial contribution from the data:
  // This is done for each earthquake

  // Set variables that are consistent across earthquakes
  //
  Type nll5 = 0.0;                         // Total contribution from data
  vector<Type> taper = exp(-lambda*depth); // Taper function
  vector<Type> spatialEffectW = A*w;       // Shared spatial effect

  for(int j=0; j<M; j++){

    // Calculate the correct spatial effects
    vector<Type> spatialEffectX = A*X.col(j);
    vector<Type> spatialEffect = spatialEffectX + spatialEffectW;
    
    // Make the untapered slips
    vector<Type> untaperedSlips = exp(mu)*exp(spatialEffect);

    // Taper the slips
    vector<Type> taperedSlips = taper*untaperedSlips;

    // Apply the okada matrix to the taperedSlips to get subsidence
    vector<Type> okadaSubsidence = -(okada(j) * taperedSlips);

    // Subsidence is normally distributed
    // Calculate the contribution of each one on the log scale
    nll5 -= sum(dnorm(subsidence(j), okadaSubsidence, v(j), true));
  }
  
  //----------------------------------------------------------------------------
  
  // Final likelihood is the product of p(theta) p(x|theta) p(y|x, theta).
  // Therefore the sum of log likelihoods
  Type nll = nll1 + nll2 + nll3 + nll4 + nll5;
  //----------------------------------------------------------------------------
  
  // Return the total negative log likelihood
  return nll;
  
}
