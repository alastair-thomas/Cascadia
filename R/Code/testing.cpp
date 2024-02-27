#include <TMB.hpp>
using namespace R_inla;
using namespace density;
using namespace Eigen;
using namespace tmbutils;

// List of matrices
// Cite: https://github.com/kaskr/adcomp/issues/96
// AT - Changed from sparse matrix to matrix
//    - Also doesn't check anything, it assumed inputs are correct
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

//------------------------------------------------------------------------------

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_STRUCT(G, LOM_t); // All the okada matrices
  DATA_STRUCT(v, LOV_t); // list of variance vectors
  
  PARAMETER(a);
  
  Type nll = 0;
  
  matrix<Type> thisG = G(1);
  Type thisTrace = thisG.trace();
  nll += thisTrace;
  
  for(int i=0; i<2; i++){
    vector<Type> thisV = v(i);
    Type thisSum = thisV.sum();
    nll += thisSum;
  }
  
  nll += a;
  
  // Return the total negative log likelihood
  return -log(nll);
  
}
