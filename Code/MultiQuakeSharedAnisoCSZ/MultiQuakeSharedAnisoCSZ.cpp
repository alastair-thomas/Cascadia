
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


// Function to create the precision matrix for the SPDE model
// Deals with anisotropy
// Cite: https://kaskr.github.io/adcomp/R__inla_8hpp_source.html#l00071
// AT - modified to take tau as well.
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
  DATA_SPARSE_MATRIX(A);          // SPDE transformation matrix
  DATA_VECTOR(maternPri);         // PC prior on matern
  DATA_VECTOR(taperPri);          // Prior on the taper
  DATA_SCALAR(anisoPri);          // Prior on the anisotropy stretch
  DATA_STRUCT(spde,spde_aniso_t); // SPDE object
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
  PARAMETER_VECTOR(logh);     // Anisotropy parameters
  //----------------------------------------------------------------------------
  
  // Transform/ extract parameters
  //
  Type tauX   = exp(logTauX);
  Type kappaX = exp(logKappaX);
  Type tauW   = exp(logTauW);
  Type kappaW = exp(logKappaW);
  Type lambda = exp(logLambda);
  Type psi    = exp(logh[0]); // anisotropic stretch
  Type thetaS = logh[1];      // anisotropic rotation
  //----------------------------------------------------------------------------
  
  // Set constants used throughout:
  //
  Type M = X.cols(); // The number of earthquakes (SHOULD BE E)
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
  
  // Calculate the contribution of the prior on anisotropic vector
  // Normal which makes psi 95% of the time in [1/2, 2]
  //
  Type ll4 = dnorm(logh[0], Type(0.0), anisoPri, true);
  Type nll4 = -ll4;
  //----------------------------------------------------------------------------
  
  // Calculate the contribution of the spatial effects via SPDE approximation:
  // This is done for each earthquake
  
  // Need to parameterize H matrix such that det(H)=1 (preserving volume) 
  // Note that H appears in (20) in Lindgren et al 2011
  
  
  // Method of eigenvalues
  matrix<Type> H1(2,2);
  H1(0,0) = sin(thetaS);
  H1(1,0) = cos(thetaS);
  H1(0,1) = cos(thetaS);
  H1(1,1) = -sin(thetaS);
  
  matrix<Type> H3 = H1.transpose();
  
  matrix<Type> H2(2,2);
  H2(0,0) = psi;
  H2(1,0) = Type(0.0);
  H2(0,1) = Type(0.0);
  H2(1,1) = Type(1.0)/psi;
  
  matrix<Type> H = H1*H2*H3;
  
  
  // Create GMRF precision matrix
  SparseMatrix<Type> Qx = spde_Q_Aniso(spde,kappaX,tauX,H);
  
  Type nll5 = 0.0;   // Initialise the spatial negative log likelihood
  
  // Loop over each earthquake
  for(int j=0; j<M; j++){
    nll5 += GMRF(Qx)(X.col(j)); // effect of jth spatial field
  }
  
  // Add the contribution for the shared spatial effect
  SparseMatrix<Type> Qw = spde_Q_Aniso(spde,kappaW,tauW,H);
  nll5 += GMRF(Qw)(w);
  //----------------------------------------------------------------------------
  
  // Calculate the non spatial contribution from the data:
  // This is done for each earthquake

  Type nll6 = 0.0; // Initialise contribution of data
  
  vector<Type> taper = exp(-lambda*depth); // Taper function
  vector<Type> spatialEffectW = A*w; // Shared spatial component
  
  // Add contribution for each earthquake
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
    nll6 -= sum(dnorm(subsidence(j), okadaSubsidence, v(j), true));
  }
  //----------------------------------------------------------------------------
  
  // Final likelihood is the product of p(theta) p(x|theta) p(y|x, theta).
  // Therefore the sum of log likelihoods
  Type nll = nll1 + nll2 + nll3 + nll4 + nll5 + nll6;
  //----------------------------------------------------------------------------
  
  // Return the total negative log likelihood
  return nll;
  
}
