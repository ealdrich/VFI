//////////////////////////////////////////////////////////////////////////////
///
/// @file main.cpp
///
/// @brief File containing main main function for the VFI problem.
///
//////////////////////////////////////////////////////////////////////////////

#include "global.h"
#include <math.h>
#include <ctime>
#include <typeinfo>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>

using namespace std;
using namespace Eigen;

//////////////////////////////////////////////////////////////////////////////
///
/// @fn main()
///
/// @brief Main function for the VFI problem.
///
/// @details Performs value function iteration on the CPU, finding the
/// maximum of the Bellman objective function for each node in the state
/// space and iterating until convergence.
///
/// @returns 0 upon successful complete, 1 otherwise.
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
int main()
{

  // admin
  int imax;
  REAL diff = 1.0;
  int i, j, l;
  double tic = curr_second(); // Start time

  // allocate variables in host memory
  VectorXR K(nk);
  VectorXR Z(nz);
  MatrixXR P(nz, nz);
  MatrixXR V0(nk, nz);
  MatrixXR V(nk, nz);
  MatrixXi G(nk, nz);

  // compute TFP grid, capital grid and initial VF
  REAL lambda = 3;
  ar1(lambda, Z, P);
  kGrid(Z, K);
  vfInit(Z, V0);

  // iterate
  int count = 0;
  bool how = false;
  while(fabs(diff) > tol){
    if(count < 3 | count % howard == 0) how = false; else how = true;
    vfStep(how, K, Z, P, V0, V, G);
    diff = (V-V0).array().abs().maxCoeff();
    V0 = V;
    ++count;
    //cout << "Iteration: " << count << ", Max Value Function Diff: " << diff << endl;
  }

  // Compute solution time
  REAL toc = curr_second();
  REAL solTime  = toc - tic;

  // write to file (column major)
  ofstream fileSolTime, fileValue, filePolicy;
  fileSolTime.open("solutionTime.dat");
  fileValue.open("valueFunc.dat");
  filePolicy.open("policyFunc.dat");
  fileSolTime << solTime << endl;
  fileValue << nk << endl;
  fileValue << nz << endl;
  filePolicy << nk << endl;
  filePolicy << nz << endl;
  for(j = 0 ; j < nz ; ++j){
    for(i = 0 ; i < nk ; ++i){
      fileValue << V(i,j) << endl;
      filePolicy << G(i,j) << endl;
    }
  }  
  fileSolTime.close();
  fileValue.close();
  filePolicy.close();

  return 0;

}
