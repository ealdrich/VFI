//////////////////////////////////////////////////////////////////////////////
///
/// @file global.h
///
/// @brief Global header file.
///
//////////////////////////////////////////////////////////////////////////////

#ifndef __FILE_GLOBALVARS_H_SEEN__
#define __FILE_GLOBALVARS_H_SEEN__

#include <thrust/device_vector.h>

typedef double REAL;

// Class for storing economic and computational parameters of the model
class parameters{
 public:
  REAL eta, beta, alpha, delta, mu, rho, sigma, lambda, tol;
  int nk, nz, howard;
  char maxtype;
  void load(char*);
};

// Function declarations
REAL curr_second (void);
void ar1(const parameters& param, thrust::device_vector<REAL>& Z,
	 thrust::device_vector<REAL>& P);
void kGrid(const parameters& param, const thrust::device_vector<REAL>& Z,
	   thrust::device_vector<REAL>& K);
void vfInit(const parameters& param, const thrust::device_vector<REAL>& Z,
	    thrust::device_vector<REAL>& V);

#endif
