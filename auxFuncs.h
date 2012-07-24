//////////////////////////////////////////////////////////////////////////////
///
/// @file auxFuncs.h
///
/// @brief Simple auxiliary functions.
///
/// @author Eric M. Aldrich \n
///         ealdrich@ucsc.edu
///
/// @version 1.0
///
/// @date 18 July 2012
///
/// @copyright Copyright Eric M. Aldrich 2012 \n
///            Distributed under the Boost Software License, Version 1.0
///            (See accompanying file LICENSE_1_0.txt or copy at \n
///            http://www.boost.org/LICENSE_1_0.txt)
///
//////////////////////////////////////////////////////////////////////////////

#ifndef __FILE_AUX_FUNCS_H_SEEN__
#define __FILE_AUX_FUNCS_H_SEEN__

//#include "global.h"
//#include <math.h>
#include <iostream>
#include <iomanip>

//////////////////////////////////////////////////////////////////////////////
///
/// @brief Function to print the elements of a matrix.
///
/// @details This functions prints a subset of the elements of a matrix to
/// the screen.
///
/// @param colMaj boolean indicating if the matrix is stored in
/// column-major format.
/// @param M number of rows in the data matrix.
/// @param N number of columns in the data matrix.
/// @param X pointer to data matrix.
/// @param printRows number of rows to print.
/// @param printCols number of columns to print.
/// @param precision number of significant digits to print.
///
/// @return Void.
///
//////////////////////////////////////////////////////////////////////////////
template<class T>
void printMatrix(const bool colMaj, const int M, const int N, T* X,
		 const int printRows, const int printCols, const int digits)
{
  std::cout.precision(digits);
  for(int ix = 0 ; ix < printRows ; ++ix){
    for(int jx = 0 ; jx < printCols ; ++jx){
      if(colMaj){
        std::cout << std::setw(2*digits) << *(X+ix+jx*M) << " ";
      } else {
        std::cout << std::setw(2*digits) << *(X+ix*N+jx) << " ";
      }
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////
///
/// @brief Function to print the elements of a vector.
///
/// @details This functions prints a subset of the elements of a vector to
/// the screen.
///
/// @param N number of elements in the data matrix.
/// @param X pointer to data matrix.
/// @param precision number of significant digits to print.
///
/// @return Void.
///
//////////////////////////////////////////////////////////////////////////////
template<class T>
void printVector(const int N, T* X, const int digits)
{
  std::cout.precision(digits);
  for(int ix = 0 ; ix < N ; ++ix){
    std::cout << *(X+ix) << std::endl;
  }
}

#endif
