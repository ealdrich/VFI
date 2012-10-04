#include "global.h"
#include <Eigen/Dense>
#include <math.h>

using namespace Eigen;

//////////////////////////////////////////////////////////////////////////////
///
/// @brief function to compute maximum of Bellman objective via binary
/// search.
///
/// @details This function finds the maximum and argmax of utility over a
/// specified subgrid of capital by using a binary search algorithm. The
/// algorithm requires concavity and cannot be used with the howard
/// improvement method. The max and argmax are stored in the value and policy
/// functions, respectively.
///
/// @param klo index corresponding to the lowest value of the capital grid
/// over which to maximize.
/// @param nksub length of the subgrid of capital (beginning at klo) over
/// which to maximize.
/// @param ydepK value of output plus capital net of depreciation.
/// @param K pointer to grid of capital values.
/// @param Exp pointer to expected value function continuation values.
/// @param V pointer to updated value function (output).
/// @param G pointer to updated policy function (output).
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
void binaryMax(const int& klo, const int& nksub, const REAL& ydepK,
	       const VectorXR& K, const VectorXR& Exp, REAL& V, int& G)
{

  // binary search to find the vf max over K'
  // we assume that the value funtion is concave in capital
  int kslo = 0;
  int kshi = nksub-1;
  int ksmid1, ksmid2;
  REAL w1, w2, w3;
  
  // case 1: capital grid has more than three values
  if(nksub > 3){
    // while the grid has 3 values or more, compute vf at midpoints
    // and revise the bounds of the grid
    while(kshi-kslo > 2){
      ksmid1 = (kslo + kshi)/2;
      ksmid2 = ksmid1+1;
      w1 = pow(ydepK-K(klo+ksmid1),1-eta)/(1-eta) + beta*Exp(ksmid1);
      w2 = pow(ydepK-K(klo+ksmid2),1-eta)/(1-eta) + beta*Exp(ksmid2);
      if(w2 > w1) kslo = ksmid1; else kshi = ksmid2;
    }
    // when the grid is reduced to three values, find the max
    if(w2 > w1){
      w1 = pow(ydepK-K(klo+kshi),1-eta)/(1-eta) + beta*Exp(kshi);
      if(w2 > w1){
	V = w2; G = klo+kslo+1;
      } else {
	V = w1; G = klo+kshi;
      }
    } else {
      w2 = pow(ydepK-K(klo+kslo),1-eta)/(1-eta) + beta*Exp(kslo);
      if(w2 > w1){
	V = w2; G = klo+kslo;
      } else {
	V = w1; G = klo+kslo+1;
      }
    }
    
    // case 2: capital grid has three values
  } else if(nksub == 3) {
    // evaluate vf at each value and determine max
    w1 = pow(ydepK-K(klo+kslo),1-eta)/(1-eta) + beta*Exp(kslo);
    w2 = pow(ydepK-K(klo+kslo+1),1-eta)/(1-eta) + beta*Exp(kslo+1);
    w3 = pow(ydepK-K(klo+kshi),1-eta)/(1-eta) + beta*Exp(kshi);
    V = w1;
    G = klo+kslo;
    if(w2 > V){V = w2; G = klo+kslo+1;}
    if(w3 > V){V = w3; G = klo+kshi;}
    
    // case 3: capital grid has one or two values
  } else {
    // evaluate vf at each value and determine max
    w1 = pow(ydepK-K(klo+kslo),1-eta)/(1-eta) + beta*Exp(kslo);
    w2 = pow(ydepK-K(klo+kshi),1-eta)/(1-eta) + beta*Exp(kshi);
    if(w2 > w1){
      V = w2; G = klo+kshi;
    } else {
      V = w1; G = klo+kslo;
    }
  }
}
