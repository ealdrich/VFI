#include "global.h"
#include "auxFuncs.h"
#include <math.h>
#include <iostream>

using namespace std;

//////////////////////////////////////////////////////////////////////////////
///
/// @brief Function to compute discrete AR1 approximation values and
/// transition matrix.
///
/// @details This function that computes a discrete AR1 approximation and
/// transition matrix using the method of Tauchen (1986).
///
/// @param lambda upper and lower bounds on the AR1 grid in terms of number
/// of standard deviations from the mean.
/// @param Z pointer to array of AR1 grid values.
/// @param P pointer to array of AR1 transition matrix values.
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
void ar1(const REAL& lambda, REAL* Z, REAL* P)
{
  int i,j;

  // grid for TFP
  const REAL sigma_z = sigma/pow(1-pow(rho,2),0.5);
  const REAL mu_z = mu/(1-rho);
  const REAL zmin = mu_z - lambda*sigma_z;
  const REAL zmax = mu_z + lambda*sigma_z;
  const REAL zstep = (zmax - zmin)/(nz-1);
  for(i = 0 ; i < nz ; ++i) Z[i] = exp(zmin + zstep*i);

  // transition matrix
  REAL normarg1, normarg2;
  for(i = 0 ; i < nz ; ++i){
    normarg1 = (zmin - mu - rho*log(Z[i]))/sigma + 0.5*zstep/sigma;
    P[i*nz] = 0.5 + 0.5*erf(normarg1/pow(2,0.5));
    P[i*nz+(nz-1)] = 1 - P[i*nz];
    for(j = 1 ; j < (nz-1) ; ++j){
      normarg1 = (log(Z[j]) - mu - rho*log(Z[i]))/sigma + 0.5*zstep/sigma;
      normarg2 = (log(Z[j]) - mu - rho*log(Z[i]))/sigma - 0.5*zstep/sigma;
      P[i*nz+j] = 0.5*erf(normarg1/pow(2,0.5)) - 0.5*erf(normarg2/pow(2,0.5));
      P[i*nz+(nz-1)] -= P[i*nz+j];
    }
  }
}
