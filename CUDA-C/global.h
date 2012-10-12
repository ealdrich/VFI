//////////////////////////////////////////////////////////////////////////////
///
/// @file global.h
///
/// @brief Global header file.
///
//////////////////////////////////////////////////////////////////////////////

#ifndef __FILE_GLOBALVARS_H_SEEN__
#define __FILE_GLOBALVARS_H_SEEN__

typedef double REAL;

// Class for storing economic and computational parameters of the model
class parameters{
 public:
  REAL eta, beta, alpha, delta, mu, rho, sigma, lambda, tol;
  int nk, nz, howard;
  char maxtype;
  void load(const char*);
};

// Function declarations
REAL curr_second (void);
void ar1(const parameters& param, REAL* Z, REAL* P);
void kGrid(const parameters& param, const REAL* Z, REAL* K);
void vfInit(const parameters& param, const REAL* Z, REAL* V);

// to determine whether single or double precision is being used
extern const float singletype;
extern const double doubletype;
extern const REAL realtype;

#endif
