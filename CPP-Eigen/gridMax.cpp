#include "global.h"
#include <Eigen/Dense>
#include <math.h>

using namespace Eigen;

//////////////////////////////////////////////////////////////////////////////
///
/// @brief function to compute maximum of Bellman objective via grid search.
///
/// @param klo index corresponding to the lowest value of the capital grid
/// over which to maximize.
/// @param nksub length of the subgrid of capital (beginning at klo) over
/// which to maximize.
/// @param ydepK value of output plus capital net of depreciation.
/// @param K pointer to grid of capital values.
/// @param Exp pointer to expected value function continuation values.
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
void gridMax(const int& klo, const int& nksub, const REAL& ydepK,
		const VectorXR& K, const VectorXR& Exp, REAL& V, int& G)
{
  REAL w = pow(ydepK-K(klo),1-eta)/(1-eta) + beta*Exp(0);
  REAL wmax = w;
  int windmax = 0;
  for(int l = 1 ; l < nksub ; ++l){
    w = pow(ydepK-K(klo+l),1-eta)/(1-eta) + beta*Exp(l);
    if(w > wmax){
      wmax = w;
      windmax = l;         }
  }
  V = wmax;
  G = klo+windmax;
}
