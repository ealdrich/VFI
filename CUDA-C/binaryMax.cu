//////////////////////////////////////////////////////////////////////////////
///
/// @file binaryMax.cpp
///
/// @brief File containing binary search maximization CUDA device function.
///
/// @author Eric M. Aldrich \n
///         ealdrich@ucsc.edu
///
/// @version 1.0
///
/// @date 23 Oct 2012
///
/// @copyright Copyright Eric M. Aldrich 2012 \n
///            Distributed under the Boost Software License, Version 1.0
///            (See accompanying file LICENSE_1_0.txt or copy at \n
///            http://www.boost.org/LICENSE_1_0.txt)
///
//////////////////////////////////////////////////////////////////////////////

#include "global.h"

//////////////////////////////////////////////////////////////////////////////
///
/// @brief CUDA device function to compute maximum of Bellman objective via
/// binary search.
///
/// @details This function finds the maximum and argmax of the Bellman
/// objective over a specified subgrid of capital by using a binary search
/// algorithm. The algorithm requires concavity.
///
/// @param [in] klo Lower index of the capital grid to begin search.
/// @param [in] nksub Number of points in the capital grid to include in
/// search.
/// @param [in] nz Length of TFP grid.
/// @param [in] ydepK value of output plus depreciated capital.
/// @param [in] eta Coefficient of relative risk aversion.
/// @param [in] beta Time discount factor.
/// @param [in] K Grid of capital values.
/// @param [in] P TFP transition matrix.
/// @param [in] V0 Current value function.
/// @param [out] V Updated value function.
/// @param [out] G Updated policy function.
///
/// @returns Void.
///
//////////////////////////////////////////////////////////////////////////////
__device__ void binaryMax(const int klo, const int nksub, const int nz,
			  const REAL ydepK, const REAL eta,
			  const REAL beta, const REAL* K, const REAL* P,
			  const REAL* V0, REAL* V, REAL* G)
{
  // binary search to find the vf max over K'
  // we assume that the value funtion is concave in capital
  int kslo, kshi, ksmid1, ksmid2, l;
  REAL Exp1, Exp2, w1, w2;
  kslo = 0;
  kshi = nksub-1;
   
  // case 1: capital grid has more than three values
  if(nksub > 3){
    // while the grid has 3 values or more, compute vf at midpoints
    // and revise the bounds of the grid
    while(kshi-kslo > 2){
      ksmid1 = (kslo + kshi)/2;
      ksmid2 = ksmid1+1;
      Exp1 = 0.0;
      Exp2 = 0.0;
      for(l = 0 ; l < nz ; ++l){
  	Exp1 += (*(P+l))*(*(V0+ksmid1*nz+l));
  	Exp2 += (*(P+l))*(*(V0+ksmid2*nz+l));
      }
      w1 = pow(ydepK-K[klo+ksmid1],1-eta)/(1-eta) + beta*Exp1;
      w2 = pow(ydepK-K[klo+ksmid2],1-eta)/(1-eta) + beta*Exp2;
      if(w2 > w1) kslo = ksmid1; else kshi = ksmid2;
    }
    // when the grid is reduced to three values, find the max
    if(w2 > w1){
      Exp1 = 0.0;
      for(l = 0 ; l < nz ; ++l) Exp1 += (*(P+l))*(*(V0+kshi*nz+l));
      w1 = pow(ydepK-K[klo+kshi],1-eta)/(1-eta) + beta*Exp1;
      if(w2 > w1){
  	*V = w2; *G = klo+kslo+1;
      } else {
  	*V = w1; *G = klo+kshi;
      }
    } else {
      Exp2 = 0.0;
      for(l = 0 ; l < nz ; ++l) Exp2 += (*(P+l))*(*(V0+kslo*nz+l));
      w2 = pow(ydepK-K[klo+kslo],1-eta)/(1-eta) + beta*Exp2;
      if(w2 > w1){
  	*V = w2; *G = klo+kslo;
      } else {
  	*V = w1; *G = klo+kslo+1;
      }
    }
  	
  // case 2: capital grid has three values
  } else if(nksub == 3) {
    // evaluate vf at each value and determine max
    Exp1 = 0.0, Exp2 = 0.0;
    REAL Exp3 = 0.0;
    for(l = 0 ; l < nz ; ++l){
      Exp1 += (*(P+l))*(*(V0+kslo*nz+l));
      Exp2 += (*(P+l))*(*(V0+kslo+1*nz+l));
      Exp3 += (*(P+l))*(*(V0+kshi*nz+l));
    }
    w1 = pow(ydepK-K[klo+kslo],1-eta)/(1-eta) + beta*Exp1;
    w2 = pow(ydepK-K[klo+kslo+1],1-eta)/(1-eta) + beta*Exp2;
    const REAL w3 = pow(ydepK-K[klo+kshi],1-eta)/(1-eta) + beta*Exp3;
    *V = w1;
    *G = klo+kslo;
    if(w2 > *V){*V = w2; *G = klo+kslo+1;}
    if(w3 > *V){ *V = w3; *G = klo+kshi;}
  	
  // case 3: capital grid has one or two values
  } else {
    Exp1 = 0.0, Exp2 = 0.0;
    for(l = 0 ; l < nz ; ++l){
      Exp1 += (*(P+l))*(*(V0+kslo*nz+l));
      Exp2 += (*(P+l))*(*(V0+kshi*nz+l));
    }
    // evaluate vf at each value and determine max
    w1 = pow(ydepK-K[klo+kslo],1-eta)/(1-eta) + beta*Exp1;
    w2 = pow(ydepK-K[klo+kshi],1-eta)/(1-eta) + beta*Exp2;
    if(w2 > w1){
      *V = w2; *G = klo+kshi;
    } else {
      *V = w1; *G = klo+kslo;
    }
  }
}
