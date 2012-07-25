#include "global.h"
#include "auxFuncs.h"
#include <math.h>
#include <iostream>
#include <typeinfo>
#include <gsl/gsl_cblas.h>
#include <stdlib.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////////
///
/// @brief function to update value function.
///
/// @details  This function performs one iteration of the value function
/// iteration algorithm, using V0 as the current value function and either
/// maximizing the LHS of the Bellman if howard = false or using the 
/// concurrent policy function as the argmax if howard = true. Monotonicity
/// of the policy function IS NOT exploited as this creates dependencies
/// among the GPU processors that are not parallelizable. Maximization is
/// performed by either a grid search or binary search algorithm.
///
/// @param howard indicates if the current iteration of the value function
/// will perform a maximization (false) or if it will simply compute the
/// new value function using the concurrent policy function (true).
/// @param K pointer to grid of capital values.
/// @param Z pointer to grid of TFP values.
/// @param P pointer to TFP transition matrix.
/// @param V0 pointer to array storing current value function.
/// @param V pointer to array of storing updated value function.
/// @param G pointer to array of storing policy function (updated if
/// howard = false.
///
/// @returns Void.
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
void vfStep(const bool& howard, const REAL* K, const REAL* Z,
	    const REAL* P, const REAL* V0, REAL* V, REAL* G)
{
  int klo, khi, nksub;
  REAL ydepK;
  REAL* Exp;
  for(int i = 0 ; i < nk ; ++i){
    for(int j = 0 ; j < nz ; ++j){

      // output and depreciated capital
      ydepK = Z[j]*pow(K[i],alpha) + (1-delta)*K[i];

      // maximize on non-howard steps
      if(howard == false){

	// impose constraints on grid for future capital
	klo = 0;
	khi = binaryVal(ydepK, nk, K); // consumption nonnegativity
	if(K[khi] > ydepK) khi -= 1;

	// further restrict capital grid via monotonicity (CPU only)
	if(i > 0){
	  if(G[(i-1)*nz+j] > klo & G[(i-1)*nz+j] < khi) klo = (int)G[(i-1)*nz+j];
	}
	nksub = khi-klo+1;

	// continuation value for subgrid
	Exp = NULL;
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
	  gridMax(klo, nksub, ydepK, K, Exp, V+i*nz+j, G+i*nz+j);
	} else if (maxtype == 'b'){
	  binaryMax(klo, nksub, ydepK, K, Exp, V+i*nz+j, G+i*nz+j);
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
