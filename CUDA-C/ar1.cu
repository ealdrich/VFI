#include "global.h"

//////////////////////////////////////////////////////////////////////////////
///
/// @brief Kernel to compute discrete AR1 approximation values and
/// transition matrix.
///
/// @details This function is a GPU kernel that computes a discrete AR1
/// approximation and transition matrix using the method of Tauchen (1986).
///
/// @param nz number of values in the AR1 approximation grid.
/// @param lambda upper and lower bounds on the AR1 grid in terms of number
/// of standard deviations from the mean.
/// @param mu mean of the AR1 process.
/// @param sigma standard deviation of the AR1 process.
/// @param rho persistence of the AR1 process.
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
__global__ void ar1(const parameters param, REAL* Z, REAL* P)
{ 

  // thread
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j;

  // basic parameters
  const int nz = param.nz;
  const REAL mu = param.mu;
  const REAL rho = param.rho;
  const REAL sigma = param.sigma;
  const REAL lambda = param.lambda;

  // grid for TFP
  const REAL sigma_z = sigma/pow(1-pow(rho,2), 0.5);
  const REAL mu_z = mu*(1/(1-rho));
  const REAL zmin = mu_z - lambda*sigma_z;
  const REAL zmax = mu_z + lambda*sigma_z;
  const REAL zstep = (zmax-zmin)/(nz-1);
  Z[i] = exp(zmin + zstep*i);

  __syncthreads();

  // transition matrix
  REAL normarg1, normarg2;
  normarg1 = (zmin - mu - rho*log(Z[i]))/sigma + 0.5*zstep/sigma;
  P[i*nz] = 0.5 + 0.5*erf(normarg1/sqrt((REAL)2));
  P[i*nz+nz-1] = 1-P[i*nz];
  for(j = 1 ; j < nz-1 ; ++j){
    normarg1 = (log(Z[j]) - mu - rho*log(Z[i]))/sigma + 0.5*zstep/sigma;
    normarg2 = (log(Z[j]) - mu - rho*log(Z[i]))/sigma - 0.5*zstep/sigma;
    P[i*nz+j] = 0.5*erf(normarg1/sqrt((REAL)2)) - 0.5*erf(normarg2/sqrt((REAL)2));
    P[i*nz+nz-1] -= P[i*nz+j];
  }
}

