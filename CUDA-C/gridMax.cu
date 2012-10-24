//////////////////////////////////////////////////////////////////////////////
///
/// @file gridMax.cpp
///
/// @brief File containing grid search maximization CUDA device function.
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

//////////////////////////////////////////////////////////////////////////////
///
/// @brief CUDA device function to compute maximum of Bellman objective via
/// grid search.
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
__device__ void gridMax(const int klo, const int nksub, const int nz,
			const REAL ydepK, const REAL eta,
			const REAL beta, const REAL* K, const REAL* P,
			const REAL* V0, REAL* V, REAL* G)
{
  REAL Exp = 0.0, w, wmax;
  int l, m, windmax;
  for(m = 0 ; m < nz ; ++m) Exp += (*(P+m))*(*(V0+m));
  w = pow(ydepK-K[klo],1-eta)/(1-eta) + beta*Exp;
  wmax = w;
  windmax = 0;
  for(l = 1 ; l < nksub ; ++l){
    Exp = 0.0;
    for(m = 0 ; m < nz ; ++m) Exp += (*(P+m))*(*(V0+l*nz+m));
    w = pow(ydepK-K[klo+l],1-eta)/(1-eta) + beta*Exp;
    if(w > wmax){
      wmax = w;
      windmax = l;
    }
  }
  *V = wmax;
  *G = klo+windmax;
}
