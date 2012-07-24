#include "global.h"
#include "binaryVal.cu"
#include "gridMax.cu"
#include "binaryMax.cu"

//////////////////////////////////////////////////////////////////////////////
///
/// @brief CUDA kernel to update value function.
///
/// @details  This CUDA kernel performs one iteration of the value function
/// iteration algorithm, using V0 as the current value function and either
/// maximizing the LHS of the Bellman if howard = false or using the 
/// concurrent policy function as the argmax if howard = true. Monotonicity
/// of the policy function IS NOT exploited as this creates dependencies
/// among the GPU processors that are not parallelizable. Maximization is
/// performed by either a grid search or binary search algorithm.
///
/// @param nk length of the capital grid.
/// @param nz length of the TFP grid.
/// @param eta risk aversion parameter.
/// @param beta discount factor.	       
/// @param alpha capital share of production.
/// @param delta depreciation rate of capital.
/// @param maxtype maximization method - 'g' corresponds to grid search
/// and 'b' corresponds to binary search.
/// @param howard indicates if the current iteration of the value function
/// will perform a maximization (false) or if it will simply compute the
/// new value function using the concurrent policy function (true).
/// @param K pointer to grid of capital values.
/// @param Z pointer to grid of TFP values.
/// @param P pointer to TFP transition matrix.
/// @param V0 pointer to array storing current value function.
/// @param V pointer to array of storing updated value function.
/// @param G pointer to array of storing policy function (updated if
/// howard = false.
///
/// @returns Void.
///
/// @author Eric M. Aldrich \n
///         ealdrich@ucsc.edu
///
/// @version 1.0
///
/// @date 24 July 2012
///
/// @copyright Copyright Eric M. Aldrich 2012 \n
///            Distributed under the Boost Software License, Version 1.0
///            (See accompanying file LICENSE_1_0.txt or copy at \n
///            http://www.boost.org/LICENSE_1_0.txt)
///
//////////////////////////////////////////////////////////////////////////////
__global__ void vfStep(const int nk, const int nz, const REAL eta,
		       const REAL beta, const REAL alpha, const REAL delta,
		       const char maxtype, const bool howard, const REAL* K,
		       const REAL* Z, const REAL* P, const REAL* V0, REAL* V,
		       REAL* G) 
{
  // thread
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  const int j = blockIdx.y * blockDim.y + threadIdx.y;

  // output and depreciated capital
  const REAL ydepK = Z[j]*pow(K[i],alpha) + (1-delta)*K[i];

  // maximize on non-howard steps
  if(howard == false){

    // impose constraints on grid for future capital
    const int klo = 0;
    int khi = binaryValGPU(ydepK, nk, K); // nonnegativity of C
    if(K[khi] > ydepK) khi -= 1;
    const int nksub = khi-klo+1;

    // maximization either via grid (g), of binary search (b)
    // if binary, turn off policy iteration (to preserve concavity)
    if(maxtype == 'g'){
      gridMaxGPU(klo, nksub, nz, ydepK, eta, beta, K, (P+j*nz), (V0+klo*nz),
		 (V+i*nz+j), (G+i*nz+j));
    } else if(maxtype == 'b'){
      binaryMaxGPU(klo, nksub, nz, ydepK, eta, beta, K, (P+j*nz), (V0+klo*nz),
		   (V+i*nz+j), (G+i*nz+j));
    }

  // iterate on the policy function on non-howard steps
  } else {
    REAL Exp = 0.0;
    for(int m = 0 ; m < nz ; ++m) Exp += P[j*nz+m]*V0[(int)G[i*nz+j]*nz+m];
    V[i*nz+j] = pow(ydepK-K[(int)G[i*nz+j]],1-eta)/(1-eta) + beta*Exp;
  }
}
