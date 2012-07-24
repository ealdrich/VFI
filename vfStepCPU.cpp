/*============================================================================

 Function      vfStepCPU

 Usage         vfStepCPU(klo, khi, nksub, kslo, kshi, ksmid1, ksmid2, w, wmax,
                         windmax, w1, w2, w3, i, j, l, ydepK, howard, K, Z, P,
		         Exp, V0, V, G)

 Arguments     klo:     reference to integer which represents the index
                        corresponding to the lowest value of the capital grid
			over which to maximize.

               khi:     reference to integer which represents the index
                        corresponding to the highest value of the capital
			grid over which to maximize.

	       nksub:   reference to the integer length of the subgrid of
	                capital (beginning at klo) over which to maximize.

	       kslo:    reference to integer representing the index
	                corresponding to the lower bound of the region within
		        the subgrid of capital which the binary search method
		        isolates for maximization.

	       kshi:    reference to integer representing the index
	                corresponding to the upper bound of the region within
		        the subgrid of capital which the binary search method
		        isolates for maximization.

	       ksmid1:  reference to integer representing the index
	                corresponding to one of two midpoints of the region
		        within the subgrid of capital which the binary search
		        method isolates for maximization.

	       ksmid1:  reference to integer representing the index
	                corresponding to one of two midpoints of the region
		        within the subgrid of capital which the binary search
		        method isolates for maximization. ksmid2 = ksmid1 + 1.

	       w:       reference to REAL used to store value of utility
	                at each value of capital subgrid (grid search).

	       wmax:    reference to REAL that tracks the max of utility
	                over the capital grid (grid search).

	       windmax: reference to integer that tracks the index of the
                        capital grid that corresponds to the argmax (grid
			search).

	       w1:      reference to REAL used to store the value of utility
	                at a particular value of capital during binary search.

	       w2:      reference to REAL used to store the value of utility
	                at a particular value of capital during binary search.

	       w3:      reference to REAL used to store the value of utility
	                at a particular value of capital during binary search.
		        This variable is only used when the subgrid of capital
		        has three values.

	       i:       reference to integer value used for incrementing.

	       j:       reference to integer value used for incrementing.

	       l:       reference to integer value used for incrementing.

	       ydepK:   reference to REAL which stores the value of output
	                plus capital net of depreciation.

	       howard:  reference to constant boolean that indicates if the
	                current iteration of the value function will perform
			a maximization (false) or if it will simply compute
			the new value function using the concurrent policy
			function (true).

	       K:       pointer to array of constant REALs which stores the
	                grid of capital values.

	       Z:       pointer to array of constant REALs which stores the
	                TFP grid values.
	               
	       P:       pointer to array of constant REALs which stores the
	                TFP transition matrix values.
	              
	       Exp:     pointer to array of REALs which stores the
	                continuation value of utility at each value of capital
			subgrid (grid search).		     

	       V0:      pointer to the array of constant REALs which
	                represents the initial value function.

	       V:       pointer to the array of REALs which represents the 
	                value function.

	       G:       pointer to the array of REALs which represents the 
	                policy function.

 Description   This function performs one iteration of the value function
               iteration algorithm, using V0 as the current value function and
	       either maximizing the LHS of the Bellman if howard = false or
	       using the concurrent policy function as the argmax if
	       howard = true. Monotonicity of the policy function is exploited
	       to restrict the region of maximization. Maximization is
	       performed by either a grid search or binary search algorithm.

 Dependencies  Global variables: eta, beta, alpha, delta, nk, nz,
                                 maxtype, howard (global.h).

               Functions:        pow (math.h);
				 cblas_(s,d)gemm, cblas_(s,d)dot (cblas.h);
	                         binaryValCPU, gridMaxCPU,
				 binaryMaxCPU (auxfuncs.h).

 Return value  void.

 =============================================================================

 Author:       Eric M. Aldrich

 Contact:      ealdrich@gmail.com

 Date:         28 July 2011

 ============================================================================*/

#include "global.h"
#include "auxfuncs.h"
#include <math.h>
#include <iostream>
#include <typeinfo>
#include <gsl/gsl_cblas.h>

using namespace std;

void vfStepCPU(int& klo, int& khi, int& nksub, int& kslo, int& kshi, int& ksmid1,
	       int& ksmid2, REAL& w, REAL& wmax, int& windmax, REAL& w1,
	       REAL& w2, REAL& w3, int& i, int& j, int& l, REAL& ydepK,
	       const bool& howard, const REAL* K,const REAL* Z,
	       const REAL* P, REAL* Exp, const REAL* V0, REAL* V,
	       REAL* G)
{
  for(i = 0 ; i < nk ; ++i){
    for(j = 0 ; j < nz ; ++j){

      // output and depreciated capital
      ydepK = Z[j]*pow(K[i],alpha) + (1-delta)*K[i];

      // maximize on non-howard steps
      if(howard == false){

	// impose constraints on grid for future capital
	klo = 0;
	khi = binaryValCPU(ydepK, nk, K); // consumption nonnegativity
	if(K[khi] > ydepK) khi -= 1;

	// further restrict capital grid via monotonicity (CPU only)
	if(i > 0){
	  if(G[(i-1)*nz+j] > klo & G[(i-1)*nz+j] < khi) klo = (int)G[(i-1)*nz+j];
	}
	nksub = khi-klo+1;

	// continuation value for subgrid
	Exp = (REAL*)realloc(Exp, nksub*sizeof(REAL));
	if(typeid(realtype) == typeid(singletype)){
	  cblas_sgemv(CblasRowMajor, CblasNoTrans, nksub, nz, 1.0, ((float*)V0+klo*nz),
		      nz, ((float*)P+j*nz), 1, 0.0, (float*)Exp, 1);
	} else if(typeid(realtype) == typeid(doubletype)){
	  cblas_dgemv(CblasRowMajor, CblasNoTrans, nksub, nz, 1.0, ((double*)V0+klo*nz),
		      nz, ((double*)P+j*nz), 1, 0.0, (double*)Exp, 1);
	}

	// maximization either via grid (g), of binary search (b)
	// if binary, turn off policy iteration (to preserve concavity)
	if(maxtype == 'g'){
	  gridMaxCPU(klo, nksub, l, w, wmax, windmax, ydepK, K, Exp, V+i*nz+j,
		     G+i*nz+j);
	} else if (maxtype == 'b'){
	  binaryMaxCPU(klo, nksub, kslo, kshi, ksmid1, ksmid2, w1, w2, w3,
		       ydepK, K, Exp, V+i*nz+j, G+i*nz+j);
	}

      // iterate on the policy function on non-howard steps
      } else {
	Exp = (REAL*)realloc(Exp, sizeof(REAL));
	if(typeid(realtype) == typeid(singletype)){
	  Exp[0] = cblas_sdot(nz, ((float*)V0+(int)G[i*nz+j]*nz), 1, ((float*)P+j*nz), 1);
	} else if(typeid(realtype) == typeid(doubletype)){
	  Exp[0] = cblas_ddot(nz, ((double*)V0+(int)G[i*nz+j]*nz), 1, ((double*)P+j*nz), 1);
	}	
	V[i*nz+j] = pow(ydepK-K[(int)G[i*nz+j]],1-eta)/(1-eta) + beta*Exp[0];
      }
    }
  }
}
