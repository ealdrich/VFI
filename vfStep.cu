/*============================================================================

 Function      vfStep

 Usage         vfStep(nk, nz, eta, beta, alpha, delta, maxtype, howard, K, Z,
                      P, V0, V, G)

 Arguments     nk:      constant integer length of the capital grid.

               nz:      constant integer length of the TFP grid.

	       eta:     constant REAL representing risk aversion.

	       beta:    constant REAL representing the discount factor.	       

	       alpha:   constant REAL representing capital share.

	       delta:   constant REAL representing the depreciation rate.

	       maxtype: constant character representing the maximization
	                method. 'g' corresponds to grid search and 'b'
			corresponds to binary search.

	       howard:  constant boolean that indicates if the current
	                iteration of the value function will perform a
			maximization (false) or if it will simply compute the
			new value function using the concurrent policy
			function (true).

	       K:       pointer to array of constant REALs which stores the
	                grid of capital values.

	       P:       pointer to constant REAL array which stores the TFP
	                transition matrix.

	       Z:       pointer to the array of constant REALs which stores
	                the TFP grid values.

	       V0:      pointer to constant REAL array which stores the
	                current value function.

	       V:       pointer to the array of REALs which represents the
	                value function.

	       G:       pointer to the array of REALs which represents the
	                policy function.
	              
 Description   This function is a CUDA kernel that performs one iteration of
               the value function iteration algorithm, using V0 as the current
	       value function and either maximizing the LHS of the Bellman if
	       howard = false or using the concurrent policy function as the
	       argmax if howard = true. Monotonicity of the policy function
	       IS NOT exploited as this creates dependencies among the GPU
	       processors that are not parallelizable. Maximization is
	       performed by either a grid search or binary search algorithm.

 Dependencies  Kernels:          binary_val (binary_val.cu),
	                         grid_max (grid_max.cu),
				 binary_max (binary_max.cu).
	                         

 Return value  void.

 =============================================================================

 Author:       Eric M. Aldrich

 Contact:      ealdrich@gmail.com

 Date:         28 July 2011

 ============================================================================*/

// kernel for updating the value function
__global__ void vfStep(const int nk, const int nz, const REAL eta,
		       const REAL beta, const REAL alpha,
		       const REAL delta, const char maxtype,
		       const bool howard, const REAL* K, const REAL* Z,
		       const REAL* P, const REAL* V0, REAL* V,
		       REAL* G) 
{
  // thread
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  const int j = blockIdx.y * blockDim.y + threadIdx.y;

  // output and depreciated capital
  const REAL ydepK = Z[j]*pow(K[i],alpha) + (1-delta)*K[i];

  // maximize on non-howard steps
  if(howard == false){

    // impose constraints on grid for future capital
    const int klo = 0;
    int khi = binary_val(ydepK, nk, K); // nonnegativity of C
    if(K[khi] > ydepK) khi -= 1;
    const int nksub = khi-klo+1;

    // maximization either via grid (g), of binary search (b)
    // if binary, turn off policy iteration (to preserve concavity)
    if(maxtype == 'g'){
      grid_max(klo, nksub, nz, ydepK, eta, beta, K, (P+j*nz), (V0+klo*nz),
	       (V+i*nz+j), (G+i*nz+j));
    } else if(maxtype == 'b'){
      binary_max(klo, nksub, nz, ydepK, eta, beta, K, (P+j*nz), (V0+klo*nz),
		 (V+i*nz+j), (G+i*nz+j));
    }

  // iterate on the policy function on non-howard steps
  } else {
    REAL Exp = 0.0;
    for(int m = 0 ; m < nz ; ++m) Exp += P[j*nz+m]*V0[(int)G[i*nz+j]*nz+m];
    V[i*nz+j] = pow(ydepK-K[(int)G[i*nz+j]],1-eta)/(1-eta) + beta*Exp;
  }
}
