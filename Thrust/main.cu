//////////////////////////////////////////////////////////////////////////////
///
/// @file main.cu
///
/// @brief File containing main function for the VFI problem.
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
#include "auxFuncs.h"
#include <iostream>
#include <fstream>
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
/// Thrust can utilize a CUDA (GPU) or OpenMP (multi-core CPU) backend.
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
//////////////////////////////////////////////////////////////////////////////
int main()
{ 

  // Admin
  int ix, jx;
  REAL diff = 1.0;

  // Load parameters
  parameters params;
  params.load("../parameters.txt");
  int nk = params.nk;
  int nz = params.nz;

  // Allocate variables in device memory
  REAL tic = curr_second(); // Start timer
  thrust::device_vector<REAL> K(nk);
  thrust::device_vector<REAL> Z(nz);
  thrust::device_vector<REAL> P(nz*nz);
  thrust::device_vector<REAL> V(nk*nz);
  thrust::device_vector<REAL> G(nk*nz);
  thrust::device_vector<REAL> V0(nk*nz);
  thrust::device_vector<int> seq_vec(nk*nz);
  thrust::sequence(seq_vec.begin(), seq_vec.end());
  thrust::device_vector<REAL>::iterator maxIter;

  // compute TFP grid, capital grid and initial VF
  ar1(params, Z, P);
  kGrid(params, Z, K);
  vfInit(params, Z, V0);

  // iterate on the value function
  int count = 0;
  bool how = false;
  while(fabs(diff) > params.tol){
    if(count < 3 | count % params.howard == 0) how = false; else how = true;
    thrust::for_each(seq_vec.begin(), seq_vec.end(),
		     vfStep<REAL>(params, how,
				  raw_pointer_cast(&K[0]), raw_pointer_cast(&Z[0]),
				  raw_pointer_cast(&P[0]), raw_pointer_cast(&V0[0]),
				  raw_pointer_cast(&V[0]), raw_pointer_cast(&G[0])));
    thrust::transform(V.begin(), V.end(), V0.begin(), V0.begin(), absDiff<REAL>());
    maxIter = thrust::max_element(V0.begin(), V0.end());
    diff = *maxIter;
    V0 = V;
    ++count;
    //cout << "Iteration: " << count << ", Diff: " << diff << endl;
  }

  // Compute solution time
  REAL toc = curr_second();
  REAL solTime  = toc - tic;

  // write to file (column major)
  ofstream fileSolTime, fileValue, filePolicy;
  fileSolTime.open("solTimeThrust.dat");
  fileValue.open("valFunThrust.dat");
  filePolicy.open("polFunThrust.dat");
  fileSolTime << solTime << endl;
  fileValue << nk << endl;
  fileValue << nz << endl;
  filePolicy << nk << endl;
  filePolicy << nz << endl;
  for(jx = 0 ; jx < nz ; ++jx){
    for(ix = 0 ; ix < nk ; ++ix){
      fileValue << V[ix+jx*nk] << endl;
      filePolicy << G[ix+jx*nk] << endl;
    }
  }  
  fileSolTime.close();
  fileValue.close();
  filePolicy.close();

  return 0;

}
