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
#include <Eigen/Dense>
#include <iostream>

using namespace std;
using namespace Eigen;

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
