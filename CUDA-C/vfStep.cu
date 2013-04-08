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

//////////////////////////////////////////////////////////////////////////////
///
/// @brief CUDA kernel to update value function.
///
/// @details This function performs one iteration of the value function
/// iteration algorithm, using V0 as the current value function, maximizing
/// the LHS of the Bellman. Maximization is performed by @link binaryMax
/// @endlink.
///
/// @param [in] param Object of class parameters.
/// @param [in] K Grid of capital values.
/// @param [in] Z Grid of TFP values.
/// @param [in] P TFP transition matrix.
/// @param [in] V0 Matrix storing current value function.
/// @param [out] V Matrix storing updated value function.
/// @param [in,out] G Matrix storing policy function.
///
/// @returns Void.
///
//////////////////////////////////////////////////////////////////////////////
__global__ void vfStep(const parameters param, const REAL* K, const REAL* Z,
		       const REAL* P, const REAL* V0, REAL* V, REAL* G) 
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

  // output and depreciated capital
  const REAL ydepK = Z[j]*pow(K[i],alpha) + (1-delta)*K[i];

  // impose constraints on grid for future capital
  const int klo = 0;
  int khi = binaryVal(ydepK, nk, K); // nonnegativity of C
  if(K[khi] > ydepK) khi -= 1;
  const int nksub = khi-klo+1;
  
  // maximization
  binaryMax(klo, nksub, nz, ydepK, eta, beta, K, (P+j*nz), (V0+klo*nz),
	    (V+i*nz+j), (G+i*nz+j));

}
