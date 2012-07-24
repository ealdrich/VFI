#include "global.h"

//////////////////////////////////////////////////////////////////////////////
///
/// @brief CUDA device function to compute maximum of Bellman objective via
/// grid search.
///
/// @details This CUDA device function finds the maximum and argmax of
/// utility over a specified subgrid of capital, by using a naive grid
/// search: computing the utility at each value of the grid. The max and
/// argmax are stored in the value and policy functions, respectively.
///
/// @param klo index corresponding to the lowest value of the capital grid
/// over which to maximize.
/// @param nksub length of the subgrid of capital (beginning at klo) over
/// which to maximize.
/// @param nz length of the TFP grid.
/// @param ydepK value of output plus capital net of depreciation.
/// @param eta risk aversion parameter.
/// @param beta discount factor.	       
/// @param K pointer to grid of capital values.
/// @param P pointer to TFP  transition matrix.
/// @param V0 pointer to current value function.
/// @param V pointer to updated value function (output).
/// @param G pointer to updated policy function (output).
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
