#include "global.h"
#include <math.h>
#include <Eigen/Dense>

using namespace Eigen;

//////////////////////////////////////////////////////////////////////////////
///
/// @brief function to compute the values of an equally spaced capital grid.
///
/// @details This function computes an equally spaced capital grid. The
/// upper and lower bounds are the deterministic steady-state values of
/// capital at the highest and lowest values of the TFP process
/// (respectively), scaled by 0.95 and 1.05 (respectively).
///
/// @param Z pointer to grid of TFP values.
/// @param K pointer to grid of capital values.
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
void kGrid(const VectorXR& Z, VectorXR& K)
{

  int i;

  // initial grid for capital
  REAL kmin = 0.95*pow((1/(alpha*Z[0]))*((1/beta)-1+delta),1/(alpha-1));
  REAL kmax = 1.05*pow((1/(alpha*Z[nz-1]))*((1/beta)-1+delta),1/(alpha-1));
  K = VectorXR::LinSpaced(nk, kmin, kmax);

}
