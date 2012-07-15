/*============================================================================

 Function      ar1

 Usage         ar1(lambda, Z, P)

 Arguments     lambda: reference to constant REAL representing the number of
                       standard deviations from the mean to place the bounds
		       of the AR1 grid.
	           
	       Z:      pointer to array of REALs which stores the AR1 grid
	               values.

	       P:      pointer to array of REALs which stores the AR1
	               transition matrix values.
	              
 Description   This function computes the grid values and transition matrix
               for a finite approximation of an AR1 process according to the
	       method of Tauchen (1986).

 Dependencies  Global variables: mu, rho, sigma, nz (globalvars.h).

 Return value  void.

 =============================================================================

 Author:       Eric M. Aldrich

 Contact:      ealdrich@gmail.com

 Date:         28 July 2011

 ============================================================================*/

#include "globalvars.h"
#include "auxfuncs.h"
#include <math.h>
#include <iostream>

using namespace std;

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
    P[i*nz] = ncdf((zmin - mu - rho*log(Z[i]))/sigma + 0.5*zstep/sigma);
    P[i*nz+(nz-1)] = 1 - P[i*nz];
    for(j = 1 ; j < (nz-1) ; ++j){
      normarg1 = (log(Z[j]) - mu - rho*log(Z[i]))/sigma + 0.5*zstep/sigma;
      normarg2 = (log(Z[j]) - mu - rho*log(Z[i]))/sigma - 0.5*zstep/sigma;
      P[i*nz+j] = ncdf(normarg1) - ncdf(normarg2);
      P[i*nz+(nz-1)] -= P[i*nz+j];
    }
  }
}
