/*============================================================================

 Function      kGridCPU

 Usage         kGridCPU(Z, K)

 Arguments     Z:     pointer to the array of constant REALs which stores
                      the AR1 grid values.	      

	       K:     pointer to the array of REALs which stores the capital
	              grid values.

 Description   This function computes the values of an equally spaced capital
               grid. The upper and lower bounds are the deterministic
	       steady-state values of capital at the highest and lowest values
	       of the TFP process (respectively), scaled by 0.95 and 1.05
	       (respectively).

 Dependencies  Global variables: beta, alpha, delta, nk, nz (global.h).

               Functions:        pow (math.h).

 Return value  void.

 =============================================================================

 Author:       Eric M. Aldrich

 Contact:      ealdrich@gmail.com

 Date:         28 July 2011

 ============================================================================*/

#include "global.h"
#include <math.h>

void kGridCPU(const REAL* Z, REAL* K)
{

  int i;

  // initial grid for capital
  REAL kmin = 0.95*pow((1/(alpha*Z[0]))*((1/beta)-1+delta),1/(alpha-1));
  REAL kmax = 1.05*pow((1/(alpha*Z[nz-1]))*((1/beta)-1+delta),1/(alpha-1));
  REAL kstep = (kmax - kmin)/(nk-1);
  for(i = 0 ; i < nk ; ++i) K[i] = kmin + i*kstep;

}
