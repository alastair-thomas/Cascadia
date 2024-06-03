
#include <TMB.hpp>
using namespace R_inla;
using namespace density;
using namespace Eigen;
using namespace tmbutils;

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
  DATA_STRUCT(spde,spde_aniso_t); // SPDE object
  //----------------------------------------------------------------------------
  
  // Parameters with domain over entire real domain
  //
  PARAMETER(logKappaX);   // The scale parameter
  PARAMETER(logTauX);     // The variance parameter
  PARAMETER(logKappaW);   // The scale parameter
  PARAMETER(logTauW);     // The variance parameter
  PARAMETER_VECTOR(logh); // Ansisotropic inputs
  //----------------------------------------------------------------------------
  
  // Extract and transform parameters
  //
  Type tauX    = exp(logTauX);
  Type kappaX  = exp(logKappaX);
  Type tauW    = exp(logTauW);
  Type kappaW  = exp(logKappaW);
  Type psi    = exp(logh[0]); // anisotropic stretch
  Type thetaS = logh[1];      // anisotropic rotation
  
  //----------------------------------------------------------------------------
  
  // Calculate the Precision Matrix

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

  // Create GMRF precision matrices
  SparseMatrix<Type> Qx = spde_Q_Aniso(spde, kappaX, tauX, H);
  SparseMatrix<Type> Qw = spde_Q_Aniso(spde, kappaW, tauW, H);
  
  REPORT(Qx);
  REPORT(Qw);
  //----------------------------------------------------------------------------
  return Type(0.0);
  
}
