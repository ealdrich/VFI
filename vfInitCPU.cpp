/*============================================================================

 Function      vfInitCPU

 Usage         vfInitCPU(Z, V)

 Arguments     Z:     pointer to the array of REALs which stores the AR1
	              grid values.

	       V:     pointer to the array of REALs which stores the initial
	              value function.

 Description   This function initializes the value function.

 Dependencies  Global Variables: eta, beta, alpha, delta, nk, nz
                                 (global.h).

	       Functions:        pow (math.h).

 Return value  void.

 =============================================================================

 Author:       Eric M. Aldrich

 Contact:      ealdrich@gmail.com

 Date:         28 July 2011

 ============================================================================*/

#include "global.h"
#include <math.h>

void vfInitCPU(const REAL* Z, REAL* V)
{ 

  int i,j;

  // initialize
  REAL Kj;
  for(i = 0 ; i < nk ; ++i){
    for(j = 0 ; j < nz ; ++j){
        Kj = pow((1/(alpha*Z[j]))*((1/beta)-1+delta),1/(alpha-1));
	V[i*nz+j] = pow(Z[j]*pow(Kj, alpha) - delta*Kj,1-eta)/(1-eta);
    }
  }

}
