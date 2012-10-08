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

REAL curr_second (void);
void ar1(const REAL& lambda, thrust::device_vector<REAL>& Z,
	 thrust::device_vector<REAL>& P);
void kGrid(const thrust::device_vector<REAL>& Z,
	   thrust::device_vector<REAL>& K);
void vfInit(const thrust::device_vector<REAL>& Z,
	    thrust::device_vector<REAL>& V);

// economic parameters
extern const REAL eta;
extern const REAL beta;
extern const REAL alpha;
extern const REAL delta;
extern const REAL mu;
extern const REAL rho;
extern const REAL sigma;

// computational parameters
extern const int nk;
extern const int nz;
extern const REAL tol;

// maximization parameters
extern const char maxtype;
extern const int howard;

#endif
