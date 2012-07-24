/*============================================================================

 Function      binaryValGPU

 Usage         binaryValGPU(x, nx, X)

 Arguments     x:  constant REAL representing the value to search for in X.
                   
	       nx: constant integer representing length of X.	           
	              
	       X:  pointer to array of constant REALs.

 Description   This is a CUDA device function that finds the index, ind, of X
               such that x <= X[ind]. We assume that X is increasing.

 Return value  integer, representing an index of the array X.

 =============================================================================

 Author:       Eric M. Aldrich

 Contact:      ealdrich@gmail.com

 Date:         28 July 2011

 ============================================================================*/

#include "global.h"

// function to find the index, ind, of X such that x <= X[ind]
// assume that X is increasing
__device__ int binaryValGPU(const REAL x, const int n, const REAL* X)
{

  int imax;

  // check if x is out of bounds
  if(x < X[0]){
    imax = 0;
    return imax;
  }
  if(x > X[n-1]){
    imax = n-1;
    return imax;
  }

  // otherwise
  int ilo, ihi, imid;
  ilo = 0;
  ihi = n-1;
  while((ihi-ilo) > 1){
    imid = (ilo + ihi)/2;
    if(X[imid] == x){
      imax = imid;
      return imax;
    } else if(X[imid] > x){
      ihi = imid;
    } else ilo = imid;
  }  
  imax = ihi;
  return imax;
}
