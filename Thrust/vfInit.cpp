//////////////////////////////////////////////////////////////////////////////
///
/// @file vfInit.cpp
///
/// @brief File containing function to initialize the value function.
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
#include <math.h>
#include <thrust/device_vector.h>

//////////////////////////////////////////////////////////////////////////////
///
/// @brief Function to initialize value function.
///
/// @details This function initializes the value function at the
/// deterministic steady state values for each level of TFP: conditional on
/// a TFP level, the deterministic steady-state value of capital is computed,
/// as well as the associated value function value.
///
/// @param [in] param Object of class parameters.
/// @param [in] Z Grid of TFP values.
/// @param [out] V Matrix of value function values.
///
/// @returns Void.
///
//////////////////////////////////////////////////////////////////////////////
void vfInit(const parameters& param, const thrust::device_vector<REAL>& Z,
	    thrust::device_vector<REAL>& V)
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
	V[ix+jx*nk] = pow(Z[jx]*pow(Kj, alpha) - delta*Kj,1-eta)/(1-eta);
    }
  }

}
