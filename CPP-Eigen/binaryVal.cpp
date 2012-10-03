#include "global.h"

//////////////////////////////////////////////////////////////////////////////
///
/// @brief function to find the index, ind, of X such that x <= X[ind].
/// We assume that X is increasing.
///
/// @param x value to search for in array X.
/// @param nx length of array X.
/// @param X pointer to array of data to search.
///
/// @return Void.
///
/// @author Eric M. Aldrich \n
///         ealdrich@ucsc.edu
///
/// @version 1.0
///
/// @date 24 July 2012
///
/// @copyright Copyright Eric M. Aldrich 2012 \n
///            Distributed under the Boost Software License, Version 1.0
///            (See accompanying file LICENSE_1_0.txt or copy at \n
///            http://www.boost.org/LICENSE_1_0.txt)
///
//////////////////////////////////////////////////////////////////////////////
int binaryVal(const REAL& x, const int& nx, const VectorXR& X)
{
  // storage for integer return value
  int imax;

  // check if x is out of bounds
  if(x < X(0)){
    imax = 0;
    return imax;
  }
  if(x > X(nx-1)){
    imax = nx-1;
    return imax;
  }

  // otherwise
  int ilo, ihi, imid;
  ilo = 0;
  ihi = nx-1;
  while((ihi-ilo) > 1){
    imid = (ilo + ihi)/2;
    if(X(imid) == x){
      imax = imid;
      return imax;
    } else if(X(imid) > x){
      ihi = imid;
    } else ilo = imid;
  }  
  imax = ihi;
  return imax;
}
