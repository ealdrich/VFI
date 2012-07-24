/*============================================================================

 Function      kGridGPU

 Usage         kGridGPU(nk, nz, beta, alpha, delta, Z, K)

 Arguments     nk:    constant integer representing number of values in
                      capital grid.
                   
	       nz:    constant integer representing number of values in TFP
	              grid.

	       beta:  constant REAL representing the discount factor.

	       alpha: constant REAL representing capital share.

	       delta: constant REAL representing the depreciation rate.

	       Z:     pointer to the array of constant REALs which stores
	              the AR1 grid values.	              

	       K:     pointer to the array of REALs which stores the capital
	              grid values.

 Description   This function is a CUDA kernel that computes the values of an
               equally spaced capital grid. The upper and lower bounds are the
	       deterministic steady-state values of capital at the highest and
	       lowest values of the TFP process (respectively), scaled by 0.95
	       and 1.05 (respectively).

 Return value  void.

 =============================================================================

 Author:       Eric M. Aldrich

 Contact:      ealdrich@gmail.com

 Date:         28 July 2011

 ============================================================================*/

#include "global.h"

__global__ void kGridGPU(const int nk, const int nz, const REAL beta,
			 const REAL alpha, const REAL delta, const REAL* Z,
			 REAL* K) 
{
  // thread
  const int i = blockIdx.x * blockDim.x + threadIdx.x;

  // grid for capital
  const REAL kmin = 0.95*pow((1/(alpha*Z[0]))*((1/beta)-1+delta),1/(alpha-1));
  const REAL kmax = 1.05*pow((1/(alpha*Z[nz-1]))*((1/beta)-1+delta),1/(alpha-1));
  const REAL kstep = (kmax-kmin)/(nk-1);
  K[i] = kmin + kstep*i;
}
