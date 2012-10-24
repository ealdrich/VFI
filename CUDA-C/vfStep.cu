//////////////////////////////////////////////////////////////////////////////
///
/// @file vfStep.cpp
///
/// @brief File containing CUDA kernel which performs the main iterative step
/// of the VFI problem.
///
/// @author Eric M. Aldrich \n
///         ealdrich@ucsc.edu
///
/// @version 1.0
///
/// @date 23 Oct 2012
///
/// @copyright Copyright Eric M. Aldrich 2012 \n
///            Distributed under the Boost Software License, Version 1.0
///            (See accompanying file LICENSE_1_0.txt or copy at \n
///            http://www.boost.org/LICENSE_1_0.txt)
///
//////////////////////////////////////////////////////////////////////////////

#include "global.h"
#include "binaryVal.cu"
#include "binaryMax.cu"
#include "gridMax.cu"

//////////////////////////////////////////////////////////////////////////////
///
/// @brief CUDA kernel to update value function.
///
/// @details This function performs one iteration of the value function
/// iteration algorithm, using V0 as the current value function and either
/// maximizing the LHS of the Bellman if @link howard @endlink = false or
/// using the concurrent policy function as the argmax if
/// @link howard @endlink = true. Maximization is performed by either
/// @link gridMax @endlink or @link binaryMax @endlink.
///
/// @param [in] param Object of class parameters.
/// @param [in] howard Indicates if the current iteration of the value
/// function will perform a maximization (false) or if it will simply compute
/// the new value function using the concurrent policy function (true).
/// @param [in] K Grid of capital values.
/// @param [in] Z Grid of TFP values.
/// @param [in] P TFP transition matrix.
/// @param [in] V0 Matrix storing current value function.
/// @param [out] V Matrix storing updated value function.
/// @param [in,out] G Matrix storing policy function (updated if
/// howard = false).
///
/// @returns Void.
///
//////////////////////////////////////////////////////////////////////////////
__global__ void vfStep(const parameters param, const bool howard,
		       const REAL* K, const REAL* Z, const REAL* P,
		       const REAL* V0, REAL* V, REAL* G) 
{
  // thread
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  const int j = blockIdx.y * blockDim.y + threadIdx.y;


  // Basic parameters
  const int nk = param.nk;
  const int nz = param.nz;
  const REAL eta = param.eta;
  const REAL beta = param.beta;
  const REAL alpha = param.alpha;
  const REAL delta = param.delta;
  const char maxtype = param.maxtype;

  // output and depreciated capital
  const REAL ydepK = Z[j]*pow(K[i],alpha) + (1-delta)*K[i];

  // maximize on non-howard steps
  if(howard == false){

    // impose constraints on grid for future capital
    const int klo = 0;
    int khi = binaryVal(ydepK, nk, K); // nonnegativity of C
    if(K[khi] > ydepK) khi -= 1;
    const int nksub = khi-klo+1;

    // maximization either via grid (g), of binary search (b)
    // if binary, turn off policy iteration (to preserve concavity)
    if(maxtype == 'g'){
      gridMax(klo, nksub, nz, ydepK, eta, beta, K, (P+j*nz), (V0+klo*nz),
		 (V+i*nz+j), (G+i*nz+j));
    } else if(maxtype == 'b'){
      binaryMax(klo, nksub, nz, ydepK, eta, beta, K, (P+j*nz), (V0+klo*nz),
		   (V+i*nz+j), (G+i*nz+j));
    }

  // iterate on the policy function on non-howard steps
  } else {
    REAL Exp = 0.0;
    for(int m = 0 ; m < nz ; ++m) Exp += P[j*nz+m]*V0[(int)G[i*nz+j]*nz+m];
    V[i*nz+j] = pow(ydepK-K[(int)G[i*nz+j]],1-eta)/(1-eta) + beta*Exp;
  }
}
