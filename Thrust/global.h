//////////////////////////////////////////////////////////////////////////////
///
/// @file global.h
///
/// @brief Global header file.
///
//////////////////////////////////////////////////////////////////////////////

#ifndef __FILE_GLOBALVARS_H_SEEN__
#define __FILE_GLOBALVARS_H_SEEN__

#define GPUSOLVE 1

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

typedef double REAL;

#if GPUSOLVE
  typedef thrust::device_vector<REAL> thrust_vectorXR;
  typedef thrust::device_vector<int> thrust_vectorXi;
  template <typename T>
  T* thrust_ptr(thrust::device_vector<T>& X){return thrust::raw_pointer_cast(&X[0]);}
#else
  typedef thrust::host_vector<REAL> thrust_vectorXR;
  typedef thrust::host_vector<int> thrust_vectorXi;
  template <typename T>
  T* thrust_ptr(thrust::host_vector<T>& X){return &X[0];}
#endif

REAL curr_second (void);
void ar1(const REAL& lambda, thrust_vectorXR& Z, thrust_vectorXR& P);
void kGrid(const thrust_vectorXR& Z, thrust_vectorXR& K);
void vfInit(const thrust_vectorXR& Z, thrust_vectorXR& V);

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
