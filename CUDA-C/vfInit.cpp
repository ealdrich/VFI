#include "global.h"
#include <math.h>

//////////////////////////////////////////////////////////////////////////////
///
/// @brief CUDA kernel to initialize value function.
///
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
void vfInit(const parameters& param, const REAL* Z, REAL* V)
{ 

  // basic parameters
  const int nk = param.nk;
  const int nz = param.nz;
  const REAL alpha = param.alpha;
  const REAL beta = param.beta;
  const REAL delta = param.delta;
  const REAL eta = param.eta;

  // initialize
  REAL Kj;
  for(int ix = 0 ; ix < nk ; ++ix){
    for(int jx = 0 ; jx < nz ; ++jx){
        Kj = pow((1/(alpha*Z[jx]))*((1/beta)-1+delta),1/(alpha-1));
	V[ix*nz+jx] = pow(Z[jx]*pow(Kj, alpha) - delta*Kj,1-eta)/(1-eta);
    }
  }

}
