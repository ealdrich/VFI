/*============================================================================

 Function      binary_val

 Usage         binary_val(x, nx, X)

 Arguments     x:  reference to constant REAL representing the value to
	           search for in X.
                   
	       nx: reference to constant integer representing length of X.	           
	              
	       X:  pointer to array of constant REALs.

 Description   This function finds the index, ind, of X such that x <= X[ind].
               We assume that X is increasing.

 Return value  integer, representing an index of the array X.

 =============================================================================

 Author:       Eric M. Aldrich

 Contact:      ealdrich@gmail.com

 Date:         28 July 2011

 ============================================================================*/

#include "globalvars.h"
#include "auxfuncs.h"

int binary_val(const REAL& x, const int& nx, const REAL* X)
{
  // storage for integer return value
  int imax;

  // check if x is out of bounds
  if(x < *X){
    imax = 0;
    return imax;
  }
  if(x > *(X+nx-1)){
    imax = nx-1;
    return imax;
  }

  // otherwise
  int ilo, ihi, imid;
  ilo = 0;
  ihi = nx-1;
  while((ihi-ilo) > 1){
    imid = (ilo + ihi)/2;
    if(*(X+imid) == x){
      imax = imid;
      return imax;
    } else if(*(X+imid) > x){
      ihi = imid;
    } else ilo = imid;
  }  
  imax = ihi;
  return imax;
}
