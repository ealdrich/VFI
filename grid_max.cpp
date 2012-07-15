/*============================================================================

 Function      grid_max

 Usage         grid_max(klo, nksub, l, w, wmax, windmax, ydepK, K, Exp, V, G)

 Arguments     klo:     reference to constant integer which represents the
                        index corresponding to the lowest value of the capital
			grid over which to maximize.

	       nksub:   reference to the constant integer length of the
	                subgrid of capital (beginning at klo) over which to
			maximize.

	       l:       reference to integer value used for incrementing.

	       w:       reference to REAL used to store value of utility
	                at each value of capital subgrid

	       wmax:    reference to REAL that tracks the max of utility
	                over the capital grid.

	       windmax: reference to integer that tracks the index of the
                        capital grid that corresponds to the argmax.

	       ydepK:   reference to constant REAL which stores the value
	                of output plus capital net of depreciation.

	       K:       pointer to array of constant REALs which stores the
	                grid of capital values.

	       Exp:     pointer to array of constant REALs which stores the
	                continuation value of utility at each value of capital
			subgrid.		     

	       V:       pointer to the array of REALs which represents
	                the value function.

	       G:       pointer to the array of REALs which represents
	                the policy function.
	              
 Description   This function finds the maximum and argmax of utility over a
               specified subgrid of capital, by using a naive grid search:
	       computing the utility at each value of the grid. The max and
	       argmax are stored in the value and policy functions,
	       respectively.

 Dependencies  Global variables: eta, beta (globalvars.h).

               Functions:        pow (math.h).

 Return value  void.

 =============================================================================

 Author:       Eric M. Aldrich

 Contact:      ealdrich@gmail.com

 Date:         28 July 2011

 ============================================================================*/

#include "globalvars.h"
#include "auxfuncs.h"
#include <math.h>

// grid search maximization function
void grid_max(const int& klo, const int& nksub, int& l, REAL& w,
	      REAL& wmax, int& windmax, const REAL& ydepK,
	      const REAL* K, const REAL* Exp, REAL* V, REAL* G)
{
  w = pow(ydepK-*(K+klo),1-eta)/(1-eta) + beta*(*Exp);
  wmax = w;
  windmax = 0;
  for(l = 1 ; l < nksub ; ++l){
    w = pow(ydepK-*(K+klo+l),1-eta)/(1-eta) + beta*(*(Exp+l));
    if(w > wmax){
      wmax = w;
      windmax = l;         }
  }
  *V = wmax;
  *G = klo+windmax;
}
