/*============================================================================

 Function      binary_max

 Usage         binary_max(klo, kslo, kshi, nksub, ksmid1, ksmid2, w1, w2,
                          w3, ydepK, K, Exp, V, G)

 Arguments     klo:    reference to constant integer which represents the
                       lowest value of the capital grid over which to
		       maximize.

	       nksub:  reference to the constant integer length of the
	               subgrid of capital (beginning at klo) over which to
		       maximize.

	       kslo:   reference to integer representing the index
	               corresponding to the lower bound of the region within
		       the subgrid of capital which the binary search method
		       isolates for maximization.

	       kshi:   reference to integer representing the index
	               corresponding to the upper bound of the region within
		       the subgrid of capital which the binary search method
		       isolates for maximization.

	       ksmid1: reference to integer representing the index
	               corresponding to one of two midpoints of the region
		       within the subgrid of capital which the binary search
		       method isolates for maximization.

	       ksmid1: reference to integer representing the index
	               corresponding to one of two midpoints of the region
		       within the subgrid of capital which the binary search
		       method isolates for maximization. ksmid2 = ksmid1 + 1.

	       w1:     reference to REAL used to store the value of utility
	               at a particular value of capital during binary search.

	       w2:     reference to REAL used to store the value of utility
	               at a particular value of capital during binary search.

	       w3:     reference to REAL used to store the value of utility
	               at a particular value of capital during binary search.
		       This variable is only used when the subgrid of capital
		       has three values.

	       ydepK:  reference to constant REAL which stores the value of
	               output plus capital net of depreciation.

	       K:      pointer to array of constant REALs which stores the
	               grid of capital values.

	       Exp:    pointer to array of constant REALs which stores the
	               continuation value of utility at each value of capital
		       subgrid.		     

	       V:      pointer to the array of REALs which represents the 
	               value function.

	       G:      pointer to the array of REALs which represents the 
	               policy function.
	              
 Description   This function finds the maximum and argmax of utility over a
               specified subgrid of capital by using a binary search
	       algorithm.  The algorithm requires concavity and cannot be
	       used with the howard improvement method.The max and argmax are
	       stored in the value and policy functions, respectively.

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

// binary search maximization function
void binary_max(const int& klo, const int& nksub, int& kslo, int& kshi,
		int& ksmid1, int& ksmid2, REAL& w1, REAL& w2, REAL& w3,
		const REAL& ydepK, const REAL* K, const REAL* Exp,
		REAL* V, REAL* G)
{

  // binary search to find the vf max over K'
  // we assume that the value funtion is concave in capital
  kslo = 0;
  kshi = nksub-1;
  
  // case 1: capital grid has more than three values
  if(nksub > 3){
    // while the grid has 3 values or more, compute vf at midpoints
    // and revise the bounds of the grid
    while(kshi-kslo > 2){
      ksmid1 = (kslo + kshi)/2;
      ksmid2 = ksmid1+1;
      w1 = pow(ydepK-*(K+klo+ksmid1),1-eta)/(1-eta) + beta*(*(Exp+ksmid1));
      w2 = pow(ydepK-*(K+klo+ksmid2),1-eta)/(1-eta) + beta*(*(Exp+ksmid2));
      if(w2 > w1) kslo = ksmid1; else kshi = ksmid2;
    }
    // when the grid is reduced to three values, find the max
    if(w2 > w1){
      w1 = pow(ydepK-*(K+klo+kshi),1-eta)/(1-eta) + beta*(*(Exp+kshi));
      if(w2 > w1){
	*V = w2; *G = klo+kslo+1;
      } else {
	*V = w1; *G = klo+kshi;
      }
    } else {
      w2 = pow(ydepK-*(K+klo+kslo),1-eta)/(1-eta) + beta*(*(Exp+kslo));
      if(w2 > w1){
	*V = w2; *G = klo+kslo;
      } else {
	*V = w1; *G = klo+kslo+1;
      }
    }
    
    // case 2: capital grid has three values
  } else if(nksub == 3) {
    // evaluate vf at each value and determine max
    w1 = pow(ydepK-*(K+klo+kslo),1-eta)/(1-eta) + beta*(*(Exp+kslo));
    w2 = pow(ydepK-*(K+klo+kslo+1),1-eta)/(1-eta) + beta*(*(Exp+kslo+1));
    w3 = pow(ydepK-*(K+klo+kshi),1-eta)/(1-eta) + beta*(*(Exp+kshi));
    *V = w1;
    *G = klo+kslo;
    if(w2 > *V){*V = w2; *G = klo+kslo+1;}
    if(w3 > *V){*V = w3; *G = klo+kshi;}
    
    // case 3: capital grid has one or two values
  } else {
    // evaluate vf at each value and determine max
    w1 = pow(ydepK-*(K+klo+kslo),1-eta)/(1-eta) + beta*(*(Exp+kslo));
    w2 = pow(ydepK-*(K+klo+kshi),1-eta)/(1-eta) + beta*(*(Exp+kshi));
    if(w2 > w1){
      *V = w2; *G = klo+kshi;
    } else {
      *V = w1; *G = klo+kslo;
    }
  }
}
