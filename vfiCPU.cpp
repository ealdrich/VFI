/*============================================================================

 Function      vfiCPU

 Usage         vfiCPU(V, G)

 Arguments     V: pointer to array of REALs storing the value function.
                   
               G: pointer to array of REALs storing the policy function.
	              
 Description   This function performs value function iteration on the CPU,
               finding the maximum of the Bellman objective function for each
	       node in the state space and iterating until convergence.

 Dependencies  Global Variables: nk, nz, tol, howard (globalvas.h).

               Functions:        pow, fabs (math.h);
	                         cblas_(s,d)axpy, cblas_i(s,d)amax (cblas.h);
	                         ar1CPU, kGridCPU, vfInitCPU,
				 vfStepCPU (auxfuncs.h).

 Return value  Returns 0 upon successful completion, 1 otherwise.

 =============================================================================

 Author:       Eric M. Aldrich

 Contact:      ealdrich@gmail.com

 Date:         28 July 2011

 ============================================================================*/

#include "global.h"
#include "auxfuncs.h"
#include <math.h>
#include <ctime>
#include <typeinfo>
#include <gsl/gsl_cblas.h>
#include <iostream>

using namespace std;

int vfiCPU(REAL* V, REAL* G)
{

  // admin
  int imax;
  REAL diff = 1.0;
  int i, j, l;
  clock_t start;

  // allocate variables in host memory
  REAL* K = new REAL[nk];
  REAL* Z = new REAL[nz];
  REAL* P = new REAL[(int)pow(nz,2)];
  REAL* V0 = new REAL[nk*nz];
  REAL* Vtemp;

  // compute TFP grid, capital grid and initial VF
  REAL lambda = 3;
  ar1CPU(lambda, Z, P);
  kGridCPU(Z, K);
  vfInitCPU(Z, V0);

  // various arguments to vfStep
  REAL* Exp = NULL;
  int klo, khi, nksub;
  REAL ydepK;
  REAL w, wmax; // if grid
  int windmax; // if grid
  int kslo, kshi, ksmid1, ksmid2; // if binary
  REAL w1, w2, w3; // if binary

  // iterate
  int count = 0;
  bool how = false;
  while(fabs(diff) > tol){
    if(count < 3 | count % howard == 0) how = false; else how = true;
    vfStepCPU(klo, khi, nksub, kslo, kshi, ksmid1, ksmid2, w, wmax, windmax, w1,
	   w2, w3, i, j, l, ydepK, how, K, Z, P, Exp, V0, V, G);
    if(typeid(realtype) == typeid(singletype)){
      cblas_saxpy(nk*nz, -1.0, (float*)V, 1, (float*)V0, 1);
      imax = cblas_isamax(nk*nz, (float*)V0, 1);
    } else if(typeid(realtype) == typeid(doubletype)){
      cblas_daxpy(nk*nz, -1.0, (double*)V, 1, (double*)V0, 1);
      imax = cblas_idamax(nk*nz, (double*)V0, 1);
    }	
    diff = *(V0+imax);
    Vtemp = V0;
    V0 = V;
    V = Vtemp;
    ++count;
  }
  //V = V0; // this assignment doesn't work sometimes (it works in GPU code)

  // i resort to a full copy in lieu of pointer reassignemnt
  for(i = 0 ; i < nk ; ++i){
    for(j = 0 ; j < nz ; ++j){
      V[i*nz+j] = V0[i*nz+j];
    }
  }

  return 0;

}
