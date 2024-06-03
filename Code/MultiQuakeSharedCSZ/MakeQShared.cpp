
#include <TMB.hpp>
using namespace R_inla;
using namespace density;
using namespace Eigen;
using namespace tmbutils;

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

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // The data passed from R
  //
  DATA_SPARSE_MATRIX(M0);         // Matrix to calculate GMRF precision matrix
  DATA_SPARSE_MATRIX(M1);         // `` ''
  DATA_SPARSE_MATRIX(M2);         // `` ''
  //----------------------------------------------------------------------------
  
  // Parameters with domain over entire real domain
  //
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
  
  //----------------------------------------------------------------------------
  
  // Build the sparse precision matrix requires for the GMRF approximation
  SparseMatrix<Type> Qx = spde_Q(kappaX, tauX, M0, M1, M2);
  SparseMatrix<Type> Qw = spde_Q(kappaW, tauW, M0, M1, M2);
  
  REPORT(Qx);
  REPORT(Qw);
  //----------------------------------------------------------------------------
  
  return Type(0.0);
  
}
