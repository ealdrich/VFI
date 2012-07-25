//////////////////////////////////////////////////////////////////////////////
///
/// @file main.cu
///
/// @brief File containing main main function for the VFI problem.
///
//////////////////////////////////////////////////////////////////////////////

#include "global.h"
#include "auxFuncs.h"
#include <iostream>
#include <ctime>
#include "functors.hpp"
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/for_each.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////////
///
/// @fn main()
///
/// @brief Main function for the VFI problem.
///
/// @details This function solves a standard neoclassical growth model with
/// value function iteration, using Thrust. Parallelization occurs at the
/// grid of values for the state space, with each thread finding the
/// maximum of the Bellman objective function for a pair of state values.
///
/// @details See Aldrich, Eric M., Jesus Fernandez-Villaverde,
/// A. Ronald Gallant and Juan F. Rubio-Ramirez (2011), "Tapping the
/// supercomputer under your desk: Solving dynamic equilibrium models with
/// graphics processors", Journal of Economic Dynamics & Control, 35, 386-393.
///
/// @see functors.hpp
///
/// @returns 0 upon successful completion, 1 otherwise.
///
/// @author Eric M. Aldrich \n
///         ealdrich@ucsc.edu
///
/// @version 1.0
///
/// @date 12 July 2012
///
/// @copyright Copyright Eric M. Aldrich 2012 \n
///            Distributed under the Boost Software License, Version 1.0
///            (See accompanying file LICENSE_1_0.txt or copy at \n
///            http://www.boost.org/LICENSE_1_0.txt)
///
//////////////////////////////////////////////////////////////////////////////
int main()
{ 

  // Admin
  //int imax;
  REAL diff = 1.0;

  // Allocate variables in device memory
  thrust::device_vector<REAL> K(nk);
  thrust::device_vector<REAL> Z(nz);
  thrust::device_vector<REAL> P(nz*nz);
  thrust::device_vector<REAL> V(nk*nz);
  thrust::device_vector<REAL> G(nk*nz);
  thrust::device_vector<REAL> V0(nk*nz);
  thrust::device_vector<REAL>::iterator maxIter;

  // Compute TFP grid (Z)
  double lambda = 3;
  thrust::counting_iterator<int> counter(0);
  thrust::transform(counter, counter+nz,
		    Z.begin(), // output destination
		    ar1Vals<REAL>(nz, lambda, mu, sigma, rho));

  // Compute transition matrix (P)
  thrust::for_each(counter, counter+nz,
		   transMat<REAL>(nz, mu, sigma, rho,
				  raw_pointer_cast(&Z[0]), raw_pointer_cast(&P[0])));

  // Compute capital grid (K)
  thrust::transform(counter, counter+nk,
		    K.begin(), // output destination
		    kGrid<REAL>(nk, nz, beta, alpha, delta, raw_pointer_cast(&Z[0])));

  // Initialize value function
  thrust::for_each(counter, counter+nz,
  		   vfInit<REAL>(nk, eta, beta, alpha, delta,
				raw_pointer_cast(&Z[0]), raw_pointer_cast(&V0[0])));

  // iterate on the value function
  int count = 0;
  bool how = false;
  REAL tic = curr_second(); // Start counting time needed to compute the solution
  while(fabs(diff) > tol){
    if(count < 3 | count % howard == 0) how = false; else how = true;
    thrust::for_each(counter, counter+nk*nz,
		     vfStep<REAL>(nk, nz, eta, beta, alpha, delta, maxtype, how,
				  raw_pointer_cast(&K[0]), raw_pointer_cast(&Z[0]),  
				  raw_pointer_cast(&P[0]), raw_pointer_cast(&V0[0]), 
				  raw_pointer_cast(&V[0]), raw_pointer_cast(&G[0])));
    thrust::transform(V.begin(), V.end(), V0.begin(), V0.begin(), abs_diff<REAL>());
    maxIter = thrust::max_element(V0.begin(), V0.end());
    diff = *maxIter;
    V0 = V;
    ++count;
    //cout << "Iteration: " << count << ", Diff: " << diff << endl;
  }
  REAL toc = curr_second();
  cout << "Solution Time: " << toc - tic << endl;
  V = V0;

  thrust::host_vector<REAL> hV = V;
  printMatrix<REAL>(1, nk, nz, &hV[0], 4, nz, 10);
  
  return 0;

}
