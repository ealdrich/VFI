#include "global.h"
#include <math.h>
#include <Eigen/Dense>

using namespace Eigen;

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
void ar1(const parameters& param, VectorXR& Z, MatrixXR& P)
{

  // basic parameters
  const int nz = param.nz;
  const REAL mu = param.mu;
  const REAL rho = param.rho;
  const REAL sigma = param.sigma;
  const REAL lambda = param.lambda;

  // grid for TFP
  const REAL sigma_z = sigma/pow(1-pow(rho,2),0.5);
  const REAL mu_z = mu/(1-rho);
  const REAL zmin = mu_z - lambda*sigma_z;
  const REAL zmax = mu_z + lambda*sigma_z;
  Z = (VectorXR::LinSpaced(nz, zmin, zmax)).array().exp().matrix();

  // transition matrix
  REAL normarg1, normarg2;
  const REAL zstep = (zmax-zmin)/(nz-1);
  VectorXR ones = VectorXR::Constant(nz,1);
  for(int ix = 0 ; ix < nz ; ++ix){
    normarg1 = (zmin - mu - rho*log(Z[ix]))/sigma + 0.5*zstep/sigma;
    P(ix,0) = 0.5 + 0.5*erf(normarg1/pow(2,0.5));
    for(int jx = 1 ; jx < (nz-1) ; ++jx){
      normarg1 = (log(Z[jx]) - mu - rho*log(Z[ix]))/sigma + 0.5*zstep/sigma;
      normarg2 = (log(Z[jx]) - mu - rho*log(Z[ix]))/sigma - 0.5*zstep/sigma;
      P(ix,jx) = 0.5*erf(normarg1/pow(2,0.5)) - 0.5*erf(normarg2/pow(2,0.5));
    }
    P.col(nz-1) = ones - P.leftCols(nz-1).rowwise().sum();
  }
}
