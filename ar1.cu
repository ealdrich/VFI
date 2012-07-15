/*============================================================================

 Function      ar1

 Usage         ar1(nz, lambda, mu, sigma, phi, Z, P)

 Arguments     nz:     constant integer representing number of values in the
                       AR1 approximation grid.

	       lambda: constant REAL representing the number of standard
	               deviations from the mean to place the bounds of the AR1
		       grid.
	           
	       mu:     constant REAL representing the mean of the AR1
	               process.

	       sigma:  constant REAL representing the standard deviation of
	               the AR1 process.

	       rho:    constant REAL representing the persistence of the AR1
	               process.

	       Z:      pointer to array of REALs which stores the AR1 grid
	               grid values.

	       P:      pointer to array of REALs which stores the AR1
	               transition matrix values.
	              
 Description   This function is a CUDA kernel that computes the grid values
               and transition matrix for a finite approximation of an AR1
	       process according to the method of Tauchen (1986).

 Return value  void.

 =============================================================================

 Author:       Eric M. Aldrich

 Contact:      ealdrich@gmail.com

 Date:         28 July 2011

 ============================================================================*/

__global__ void ar1(const int nz, const REAL lambda, const REAL mu,
		    const REAL sigma, const REAL rho, REAL* Z, REAL* P)
{ 

  // thread
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j;

  // grid for TFP
  const REAL sigma_z = sigma/sqrt(1-pow(rho,2));
  const REAL mu_z = mu*(1/(1-rho));
  const REAL zmin = mu_z - lambda*sigma_z;
  const REAL zmax = mu_z + lambda*sigma_z;
  const REAL zstep = (zmax-zmin)/(nz-1);
  Z[i] = exp(zmin + zstep*i);

  __syncthreads();

  // transition matrix
  REAL normarg1, normarg2;
  P[i*nz] = ncdfgpu((zmin - mu - rho*log(Z[i]))/sigma + 0.5*zstep/sigma);
  P[i*nz+nz-1] = 1-P[i*nz];
  for(j = 1 ; j < nz-1 ; ++j){
    normarg1 = (log(Z[j]) - mu - rho*log(Z[i]))/sigma + 0.5*zstep/sigma;
    normarg2 = (log(Z[j]) - mu - rho*log(Z[i]))/sigma - 0.5*zstep/sigma;
    P[i*nz+j] = ncdfgpu(normarg1) - ncdfgpu(normarg2);
    P[i*nz+nz-1] -= P[i*nz+j];
  }
}

