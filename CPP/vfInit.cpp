#include "global.h"
#include <math.h>
#include <Eigen/Dense>
#include <iostream>

using namespace std;
using namespace Eigen;

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
void vfInit(const parameters& param, const VectorXR& Z, MatrixXR& V)
{ 

  // basic parameters
  const int nk = param.nk;
  const REAL alpha = param.alpha;
  const REAL beta = param.beta;
  const REAL delta = param.delta;
  const REAL eta = param.eta;

  // initialize
  ArrayXR Kj = (((alpha*Z).array().pow(-1))*((1/beta)-1+delta)).pow(1/(alpha-1));
  V.row(0) = (((Z.array()*Kj.pow(alpha) - delta*Kj).pow(1-eta))/(1-eta)).matrix();
  for(int ix = 1 ; ix < nk ; ++ix) V.row(ix) = V.row(0);

}
