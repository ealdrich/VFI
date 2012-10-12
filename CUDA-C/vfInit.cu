#include "global.h"

//////////////////////////////////////////////////////////////////////////////
///
/// @brief CUDA kernel to initialize value function.
///
/// @param nz length of the TFP grid.
/// @param eta risk aversion parameter.
/// @param beta discount factor.	       
/// @param alpha capital share of production.
/// @param delta depreciation rate of capital.
/// @param Z pointer to grid of TFP values.
/// @param V pointer to array of value function values.
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
__global__ void vfInit(const parameters param, const REAL* Z, REAL* V)
{ 
  // thread
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  const int j = blockIdx.y * blockDim.y + threadIdx.y;

  // basic parameters
  const int nz = param.nz;
  const REAL alpha = param.alpha;
  const REAL beta = param.beta;
  const REAL delta = param.delta;
  const REAL eta = param.eta;

  // initialize
  const REAL Kj = pow((1/(alpha*Z[j]))*((1/beta)-1+delta),1/(alpha-1));
  V[i*nz+j] = pow(Z[j]*pow(Kj, alpha) - delta*Kj,1-eta)/(1-eta);
}
