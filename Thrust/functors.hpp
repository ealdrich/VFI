//////////////////////////////////////////////////////////////////////////////
///
/// @file functors.hpp
///
/// @brief File of Thrust functors and functions.
///
/// @author Eric M. Aldrich \n
///         ealdrich@ucsc.edu
///
/// @version 1.0
///
/// @date 12 July 2012
///
/// @copyright Copyright Eric M. Aldrich 2012 \n
///            Distributed under the Boost Software License, Version 1.0
///            (See accompanying file LICENSE_1_0.txt or copy at \n
///            http://www.boost.org/LICENSE_1_0.txt)
///
//////////////////////////////////////////////////////////////////////////////

#ifndef _FUNCTORS_H_
#define _FUNCTORS_H_

#include <thrust/iterator/zip_iterator.h>
#include <thrust/for_each.h>
#include <thrust/device_vector.h>
#include <cmath>
#include "global.h"
#include <stdio.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////////
///
/// @brief Functor to update the value function.
///
/// @details This functor performs one iteration of the value function
/// iteration algorithm, using V0 as the current value function, maximizing
/// the LHS of the Bellman. Maximization is performed by @link binaryMax
/// @endlink.
///
//////////////////////////////////////////////////////////////////////////////
template <typename T>
struct vfStep
{
  // Attributes
  const parameters params; ///< Object containing parameters.
  const T* K; ///< Pointer to capital grid.
  const T* Z; ///< Pointer to AR1 (TFP) grid.
  const T* P; ///< Pointer to transition matrix.
  const T* V0; ///< Pointer to current iteration of the value function.
  T* V; ///< Pointer to the updated value function.
  T* G; ///< Pointer to current iteration of the capital policy function.

  /// Constructor
  vfStep(parameters _params, T* _K, T* _Z, T* _P, T* _V0, T* _V, T* _G)
    : params(_params), K(_K), Z(_Z), P(_P), V0(_V0), V(_V), G(_G) {}

  /// Kernel to update the value function.
  /// @param hx index of V0 (stored as a flat array).
  /// @return Void.
  __host__ __device__
  void operator()(const int& hx) const 
  {

    // Basic parameters
    const int nk = params.nk;
    const int nz = params.nz;
    const REAL eta = params.eta;
    const REAL beta = params.beta;
    const REAL alpha = params.alpha;
    const REAL delta = params.delta;
    
    // Compute the row and column IDs
    int ix = hx%nk;
    int jx = (hx-ix)/nk;

    // Output and depreciated capital
    const T ydepK = Z[jx]*pow(K[ix],alpha) + (1-delta)*K[ix];

    // impose constraints on grid for future capital
    const int klo = 0;
    int khi = binaryVal(ydepK, nk, K); // nonnegativity of C
    if(K[khi] > ydepK) khi -= 1;
    const int nksub = khi-klo+1;

    // maximization
    binaryMax(klo, nksub, nk, nz, ydepK, eta, beta, K, P+jx,
	      V0+klo, V+ix+jx*nk, G+ix+jx*nk);

  }
};

//////////////////////////////////////////////////////////////////////////////
///
/// @brief Device function to find the location of a value in a monotonic
/// grid.
///
/// @details This function finds the first value X[ix] such that x <= X[ix],
/// where x is a scalar value, X is a monotonic array, and ix is the index
/// of X.
///
/// @param [in] x Value to search for in vector X.
/// @param [in] nx Length of array X.
/// @param [in] X Vector of data to search.
///
/// @return imax Integer ix (<= nx) such that x <= X[ix].
///
//////////////////////////////////////////////////////////////////////////////
template <typename T>
__host__ __device__
int binaryVal(const T x, const int nx, const T* X)
{

  int imax;

  // check if x is out of bounds
  if(x < X[0]){
    imax = 0;
    return imax;
  }
  if(x > X[nx-1]){
    imax = nx-1;
    return imax;
  }

  // otherwise
  int ilo, ihi, imid;
  ilo = 0;
  ihi = nx-1;
  while((ihi-ilo) > 1){
    imid = (ilo + ihi)/2;
    if(X[imid] == x){
      imax = imid;
      return imax;
    } else if(X[imid] > x){
      ihi = imid;
    } else ilo = imid;
  }  
  imax = ihi;
  return imax;
}

//////////////////////////////////////////////////////////////////////////////
///
/// @brief Device function to compute maximum of Bellman objective via grid
/// search.
///
/// @details This function finds the maximum and argmax of the Bellman
/// objective function by using a naive grid search: computing the utility
/// at each value of the grid.
///
/// @param [in] klo Lower index of the capital grid to begin search.
/// @param [in] nksub Number of points in the capital grid to include in
/// search.
/// @param [in] nz Length of TFP grid.
/// @param [in] ydepK value of output plus depreciated capital.
/// @param [in] eta Coefficient of relative risk aversion.
/// @param [in] beta Time discount factor.
/// @param [in] K Grid of capital values.
/// @param [in] P TFP transition matrix.
/// @param [in] V0 Current value function.
/// @param [out] V Updated value function.
/// @param [out] G Updated policy function.
///
/// @returns Void.
///
//////////////////////////////////////////////////////////////////////////////
template <typename T>
__host__ __device__
void gridMax(const int klo, const int nksub, const int nk,
	      const int nz, const T ydepK, const T eta,
	      const T beta, const T* K, const T* P, const T* V0, T* V, T* G)
{
  T Exp = 0.0, w, wmax;
  int l, m, windmax;
  for(m = 0 ; m < nz ; ++m) Exp += (*(P+m*nz))*(*(V0+m*nk));
  w = pow(ydepK-K[klo],1-eta)/(1-eta) + beta*Exp;
  wmax = w;
  windmax = 0;
  for(l = 1 ; l < nksub ; ++l){
    Exp = 0.0;
    for(m = 0 ; m < nz ; ++m) Exp += (*(P+m*nz))*(*(V0+l+m*nk));
    w = pow(ydepK-K[klo+l],1-eta)/(1-eta) + beta*Exp;
    if(w > wmax){
      wmax = w;
      windmax = l;
    }
  }
  *V = wmax;
  *G = klo+windmax;
}

//////////////////////////////////////////////////////////////////////////////
///
/// @brief Device function to compute maximum of Bellman objective via binary
/// search.
///
/// @details This function finds the maximum and argmax of the Bellman
/// objective over a specified subgrid of capital by using a binary search
/// algorithm. The algorithm requires concavity.
///
/// @param [in] klo Lower index of the capital grid to begin search.
/// @param [in] nksub Number of points in the capital grid to include in
/// search.
/// @param [in] nz Length of TFP grid.
/// @param [in] ydepK value of output plus depreciated capital.
/// @param [in] eta Coefficient of relative risk aversion.
/// @param [in] beta Time discount factor.
/// @param [in] K Grid of capital values.
/// @param [in] P TFP transition matrix.
/// @param [in] V0 Current value function.
/// @param [out] V Updated value function.
/// @param [out] G Updated policy function.
///
/// @returns Void.
///
//////////////////////////////////////////////////////////////////////////////
template <typename T>
__host__ __device__
void binaryMax(const int klo, const int nksub, const int nk,
		const int nz, const T ydepK, const T eta,
		const T beta, const T* K, const T* P,const T* V0, T* V, T* G)
{
  // binary search to find the vf max over K'
  // we assume that the value funtion is concave in capital
  int kslo, kshi, ksmid1, ksmid2, l;
  T Exp1, Exp2, w1, w2;
  kslo = 0;
  kshi = nksub-1;
   
  // case 1: capital grid has more than three values
  if(nksub > 3){
    // while the grid has 3 values or more, compute vf at midpoints
    // and revise the bounds of the grid
    while(kshi-kslo > 2){
      ksmid1 = (kslo + kshi)/2;
      ksmid2 = ksmid1+1;
      Exp1 = 0.0;
      Exp2 = 0.0;
      for(l = 0 ; l < nz ; ++l){
  	Exp1 += (*(P+l*nz))*(*(V0+ksmid1+l*nk));
  	Exp2 += (*(P+l*nz))*(*(V0+ksmid2+l*nk));
      }
      w1 = pow(ydepK-K[klo+ksmid1],1-eta)/(1-eta) + beta*Exp1;
      w2 = pow(ydepK-K[klo+ksmid2],1-eta)/(1-eta) + beta*Exp2;
      if(w2 > w1) kslo = ksmid1; else kshi = ksmid2;
    }
    // when the grid is reduced to three values, find the max
    if(w2 > w1){
      Exp1 = 0.0;
      for(l = 0 ; l < nz ; ++l) Exp1 += (*(P+l*nz))*(*(V0+kshi+l*nk));
      w1 = pow(ydepK-K[klo+kshi],1-eta)/(1-eta) + beta*Exp1;
      if(w2 > w1){
  	*V = w2; *G = klo+kslo+1;
      } else {
  	*V = w1; *G = klo+kshi;
      }
    } else {
      Exp2 = 0.0;
      for(l = 0 ; l < nz ; ++l) Exp2 += (*(P+l*nz))*(*(V0+kslo+l*nk));
      w2 = pow(ydepK-K[klo+kslo],1-eta)/(1-eta) + beta*Exp2;
      if(w2 > w1){
  	*V = w2; *G = klo+kslo;
      } else {
  	*V = w1; *G = klo+kslo+1;
      }
    }
  	
  // case 2: capital grid has three values
  } else if(nksub == 3) {
    // evaluate vf at each value and determine max
    Exp1 = 0.0, Exp2 = 0.0;
    T Exp3 = 0.0;
    for(l = 0 ; l < nz ; ++l){
      Exp1 += (*(P+l*nz))*(*(V0+kslo+l*nk));
      Exp2 += (*(P+l*nz))*(*(V0+kslo+1+l*nk));
      Exp3 += (*(P+l*nz))*(*(V0+kshi+l*nk));
    }
    w1 = pow(ydepK-K[klo+kslo],1-eta)/(1-eta) + beta*Exp1;
    w2 = pow(ydepK-K[klo+kslo+1],1-eta)/(1-eta) + beta*Exp2;
    const T w3 = pow(ydepK-K[klo+kshi],1-eta)/(1-eta) + beta*Exp3;
    *V = w1;
    *G = klo+kslo;
    if(w2 > *V){*V = w2; *G = klo+kslo+1;}
    if(w3 > *V){ *V = w3; *G = klo+kshi;}
  	
  // case 3: capital grid has one or two values
  } else {
    Exp1 = 0.0, Exp2 = 0.0;
    for(l = 0 ; l < nz ; ++l){
      Exp1 += (*(P+l*nz))*(*(V0+kslo+l*nk));
      Exp2 += (*(P+l*nz))*(*(V0+kshi+l*nk));
    }
    // evaluate vf at each value and determine max
    w1 = pow(ydepK-K[klo+kslo],1-eta)/(1-eta) + beta*Exp1;
    w2 = pow(ydepK-K[klo+kshi],1-eta)/(1-eta) + beta*Exp2;
    if(w2 > w1){
      *V = w2; *G = klo+kshi;
    } else {
      *V = w1; *G = klo+kslo;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
///
/// @brief Functor to compute the absolute difference between elements of
/// two vectors.
///
//////////////////////////////////////////////////////////////////////////////
template <typename T>
struct absDiff
{

  /// Kernel to compute the absolute difference between elements.
  /// @param x value of first vector element.
  /// @param y value of second vector element.
  /// @return absolute difference between elements.
  __host__ __device__
  T operator()(const T& x, const T& y) const { 
    return fabs(x - y);
  }
};

#endif
