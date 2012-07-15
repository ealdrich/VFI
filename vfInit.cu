/*============================================================================

 Function      vfInit

 Usage         vfInit(nz, eta, beta, alpha, delta, Z, V)

 Arguments     nz:    constant integer representing number of values in TFP
                      grid.

	       eta:   constant REAL representing risk aversion.

	       beta:  constant REAL representing the discount factor.

	       alpha: constant REAL representing capital share.

	       delta: constant REAL representing the depreciation rate.

	       Z:     pointer to the array of constant REALs which stores
	              the AR1 grid values.

	       V:     pointer to the array of REALs which stores the initial
	              value function.

 Description   This function is a CUDA kernel that initializes the value
               function.

 Return value  void.

 =============================================================================

 Author:       Eric M. Aldrich

 Contact:      ealdrich@gmail.com

 Date:         28 July 2011

 ============================================================================*/

__global__ void vfInit(const int nz,  const REAL eta, const REAL beta,
		       const REAL alpha, const REAL delta, const REAL* Z,
		       REAL* V)
{ 
  // thread
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  const int j = blockIdx.y * blockDim.y + threadIdx.y;

  // initialize
  const REAL Kj = pow((1/(alpha*Z[j]))*((1/beta)-1+delta),1/(alpha-1));
  V[i*nz+j] = pow(Z[j]*pow(Kj, alpha) - delta*Kj,1-eta)/(1-eta);
}
