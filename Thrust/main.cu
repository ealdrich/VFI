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
#include <thrust/host_vector.h>
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
  REAL tic = curr_second(); // Start timer
  thrust_vectorXR K(nk);
  thrust_vectorXR Z(nz);
  thrust_vectorXR P(nz*nz);
  thrust_vectorXR V(nk*nz);
  thrust_vectorXR G(nk*nz);
  thrust_vectorXR V0(nk*nz);
  thrust_vectorXi seq_vec(nk*nz);
  thrust::sequence(seq_vec.begin(), seq_vec.end());
  thrust_vectorXR::iterator maxIter;

  // Compute TFP grid (Z)
  REAL lambda = 3;
  ar1(lambda, Z, P);
  kGrid(Z, K);
  vfInit(Z, V0);

  // iterate on the value function
  int count = 0;
  bool how = false;
  while(fabs(diff) > tol){
    if(count < 3 | count % howard == 0) how = false; else how = true;
    thrust::for_each(seq_vec.begin(), seq_vec.end(),
		     vfStep<REAL>(nk, nz, eta, beta, alpha, delta, maxtype, how,
				  thrust_ptr(K), thrust_ptr(Z), thrust_ptr(P),
				  thrust_ptr(V0), thrust_ptr(V), thrust_ptr(G)));
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
