/*============================================================================

 Function      gridMaxGPU

 Usage         gridMaxGPU(klo, nksub, nz, ydepK, eta, beta, K, P, V0, V, G)

 Arguments     klo:   constant integer which represents the index
                      corresponding to the lowest value of the capital grid
		      over which to maximize.

	       nksub: constant integer length of the subgrid of capital
	              (beginning at klo) over which to maximize.

	       nz:    constant integer length of the TFP grid.

	       ydepK: constant REAL which stores the value of  output plus
	              capital net of depreciation.

	       eta:   constant REAL representing risk aversion.

	       beta:  constant REAL representing the discount factor.	       

	       K:     pointer to array of constant REALs which stores the
	              grid of capital values.

	       P:     pointer to constant REAL array which stores the TFP
	              transition matrix.

	       V0:    pointer to constant REAL array which stores the
	              current value function.

	       V:     pointer to the array of REALs which represents the
	              value function.

	       G:     pointer to the array of REALs which represents the
	              policy function.
	              
 Description   This is a CUDA device function that finds the maximum and
               argmax of utility over a specified subgrid of capital, by using
	       a naive grid search: computing the utility at each value of the
	       grid. The max and argmax are stored in the value and policy
	       functions, respectively.

 Return value  void.

 =============================================================================

 Author:       Eric M. Aldrich

 Contact:      ealdrich@gmail.com

 Date:         28 July 2011

 ============================================================================*/

#include "global.h"

__device__ void gridMaxGPU(const int klo, const int nksub, const int nz,
			   const REAL ydepK, const REAL eta,
			   const REAL beta, const REAL* K, const REAL* P,
			   const REAL* V0, REAL* V, REAL* G)
{
  REAL Exp = 0.0, w, wmax;
  int l, m, windmax;
  for(m = 0 ; m < nz ; ++m) Exp += (*(P+m))*(*(V0+m));
  w = pow(ydepK-K[klo],1-eta)/(1-eta) + beta*Exp;
  wmax = w;
  windmax = 0;
  for(l = 1 ; l < nksub ; ++l){
    Exp = 0.0;
    for(m = 0 ; m < nz ; ++m) Exp += (*(P+m))*(*(V0+l*nz+m));
    w = pow(ydepK-K[klo+l],1-eta)/(1-eta) + beta*Exp;
    if(w > wmax){
      wmax = w;
      windmax = l;
    }
  }
  *V = wmax;
  *G = klo+windmax;
}
