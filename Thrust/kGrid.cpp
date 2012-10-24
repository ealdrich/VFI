//////////////////////////////////////////////////////////////////////////////
///
/// @file kGrid.cpp
///
/// @brief File containing function to create capital grid.
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
/// @brief Function to compute the values of an equally spaced capital grid.
///
/// @details This function computes an equally spaced capital grid. The
/// upper and lower bounds are the deterministic steady-state values of
/// capital at the highest and lowest values of the TFP process
/// (respectively), scaled by 0.95 and 1.05 (respectively).
///
/// @param [in] param Object of class parameters.
/// @param [in] Z Grid of TFP values.
/// @param [out] K Grid of capital values.
///
/// @returns Void.
///
//////////////////////////////////////////////////////////////////////////////
void kGrid(const parameters& param, const thrust::device_vector<REAL>& Z,
	   thrust::device_vector<REAL>& K)
{

  // basic parameters
  const int nk = param.nk;
  const int nz = param.nz;
  const REAL alpha = param.alpha;
  const REAL beta = param.beta;
  const REAL delta = param.delta;

  // initial grid for capital
  REAL kmin = 0.95*pow((1/(alpha*Z[0]))*((1/beta)-1+delta),1/(alpha-1));
  REAL kmax = 1.05*pow((1/(alpha*Z[nz-1]))*((1/beta)-1+delta),1/(alpha-1));
  REAL kstep = (kmax - kmin)/(nk-1);
  for(int ix = 0 ; ix < nk ; ++ix) K[ix] = kmin + ix*kstep;

}
