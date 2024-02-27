
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
  Q = pow(tau, 2.)  * (kappa4*M0 + Type(2.0)*kappa2*M1 + M2);
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
  DATA_IVECTOR(spdeIDX);          // SPDE index's
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
  //PARAMETER(logitPhi);        // Contribution of shared spatial effect
  PARAMETER(logLambda);       // Controls the taper effect
  PARAMETER(mu);              // The mean untapered slip
  PARAMETER(logKappaX);       // The individual scale parameter
  PARAMETER(logTauX);         // The individual variance parameter
  PARAMETER(logKappaW);       // The shared scale parameter
  PARAMETER(logTauW);         // The shared variance parameter
  //----------------------------------------------------------------------------
  
  // Transform parameters
  //
  Type tauX = exp(logTauX);
  Type kappaX = exp(logKappaX);
  Type tauW = exp(logTauW);
  Type kappaW = exp(logKappaW);
  Type lambda = exp(logLambda);
  //Type phi = invlogit(logitPhi);
  
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
  
  // Calculate the contribution of the prior on phi
  // phi ~ Unif[0,1]
  // y = logit(phi) = log(phi / 1-phi)
  // phi = invlogit(y) = 1 / (1 + exp(-y))
  // P(phi) = 1
  // P(logit(phi)) = P(phi) * d/dy (invlogit(y))
  // d/dy (invlogit(y)) = -(-exp(-y))(1+exp(-y))^(-2) = exp(-y)/ (1 + exp(-y))^2
  // So this derivative is the Jacobian factor
  //
  //Type nll4 = 0.0;
  //nll4 += logitPhi + Type(2.0)*log(Type(1.0) + exp(-logitPhi));
  //----------------------------------------------------------------------------
  
  // Calculate the contribution of the prior on mu.
  // Currently set to be a flat prior.
  //
  // 
  //----------------------------------------------------------------------------
  
  // Calculate the contribution of the spatial effects via SPDE approximation:
  // This is done for each earthquake
  // And then for the shared spatial field

  // Build the sparse precision matrix requires for the GMRF approximation
  SparseMatrix<Type> Qx = spde_Q(kappaX, tauX, M0, M1, M2);
  
  Type M = X.cols();   // The number of earthquakes
  Type nll5 = 0.0;     // Initialise the spatial negative log likelihood
  
  // Loop over each earthquake
  for(int j=0; j<M; j++){
    vector<Type> thisX = X.col(j); // extract the jth spatial field
    nll5 += GMRF(Qx)(thisX);
  }
  
  // Add the contribution for the shared spatial effect
  SparseMatrix<Type> Qw = spde_Q(kappaW, tauW, M0, M1, M2);
  nll5 += GMRF(Qw)(w);
  //----------------------------------------------------------------------------
  
  // Calculate the non spatial contribution from the data:
  // This is done for each earthquake

  // Set variables that are consistent across earthquakes
  //
  Type nll6 = 0.0;                         // Total contribution from data
  vector<Type> taper = exp(-lambda*depth); // Taper function

  // Add contribution for each earthquake
  for(int j=0; j<M; j++){

    // Extract all the required information about this earthquake
    //
    vector<Type> thisX = X.col(j);               // Spatial effect
    matrix<Type> thisOkada = okada(j);           // Okada matrix
    vector<Type> thisSubsidence = subsidence(j); // Subsidence data
    vector<Type> thisV = v(j);                   // Variance data

    // Create the untapered slips
    vector<Type> untaperedSlips(spdeIDX.size());
    untaperedSlips.setZero();

    // Calculate the untapered slips based on the spatial fields
    for(int i=0; i<spdeIDX.size(); i++){
      Type x_ji = thisX(spdeIDX(i)); // jth earthquake, ith subfault
      Type w_i  = w(spdeIDX(i));     // ith subfault, shared component
      untaperedSlips(i) = exp(mu + w_i + x_ji);
    }

    // Taper the slips
    vector<Type> taperedSlips = taper*untaperedSlips;

    // Apply the okada matrix to the taperedSlips to get subsidence
    vector<Type> okadaSubsidence = thisOkada * taperedSlips;

    // Subsidence is normally distributed
    // Calculate the contribution of each one on the log scale
    nll6 -= sum(dnorm(thisSubsidence, okadaSubsidence, thisV, true));
  }
  //----------------------------------------------------------------------------
  
  // Final likelihood is the product of p(theta) p(x|theta) p(y|x, theta).
  // Therefore the sum of log likelihoods
  Type nll = nll1 + nll2 + nll3 + nll5 + nll6;
  //----------------------------------------------------------------------------
  
  // Calculate useful re-parameterisations
  //
  // Calculate an understandable range parameter
  double nu = 1.0;            // nu = alpha-d/2 = 2-2/2 by eqn (2) in Lindgren 
  Type rhoX = sqrt(8.0*nu)/kappaX;  // Distance where correlation is ~0.1
  
  // Calculate the marginal variance
  Type sigma2X = Type(1.0) / (4.0*3.14159265359*pow(kappaX, 2.)*pow(tauX, 2.));
  //----------------------------------------------------------------------------
  
  // Report values for testing
  REPORT(nll1); // spde prior X
  REPORT(nll2); // spde prior W
  REPORT(nll3); // taper prior 
  //REPORT(nll4); // shared spatial effect contribution prior
  REPORT(nll5); // spatial field
  REPORT(nll6); // data contribution
  //----------------------------------------------------------------------------
  
  ADREPORT(rhoX);
  ADREPORT(sigma2X);
  //----------------------------------------------------------------------------
  
  // Return the total negative log likelihood
  return nll;
  
}
