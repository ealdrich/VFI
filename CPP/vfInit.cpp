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
void vfInit(const REAL* Z, REAL* V)
{ 

  int i,j;

  // initialize
  REAL Kj;
  for(i = 0 ; i < nk ; ++i){
    for(j = 0 ; j < nz ; ++j){
        Kj = pow((1/(alpha*Z[j]))*((1/beta)-1+delta),1/(alpha-1));
	V[i*nz+j] = pow(Z[j]*pow(Kj, alpha) - delta*Kj,1-eta)/(1-eta);
    }
  }

}
