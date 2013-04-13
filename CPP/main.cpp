//////////////////////////////////////////////////////////////////////////////
///
/// @file main.cpp
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
/// @details This function solves a standard neoclassical growth model with
/// value function iteration on a CPU.
///
/// @details See Aldrich, Eric M., Jesus Fernandez-Villaverde,
/// A. Ronald Gallant and Juan F. Rubio-Ramirez (2011), "Tapping the
/// supercomputer under your desk: Solving dynamic equilibrium models with
/// graphics processors", Journal of Economic Dynamics & Control, 35, 386-393.
///
/// @returns 0 upon successful completion, 1 otherwise.
///
//////////////////////////////////////////////////////////////////////////////
int main()
{

  // admin
  int imax;
  REAL diff = 1.0;
  int i, j, l;
  double tic = curr_second(); // Start time

  // Load parameters
  parameters params;
  params.load("../parameters.txt");
  int nk = params.nk;
  int nz = params.nz;

  // allocate variables in host memory
  VectorXR K(nk);
  VectorXR Z(nz);
  MatrixXR P(nz, nz);
  MatrixXR V0(nk, nz);
  MatrixXR V(nk, nz);
  MatrixXi G(nk, nz);

  // compute TFP grid, capital grid and initial VF
  ar1(params, Z, P);
  kGrid(params, Z, K);
  vfInit(params, Z, V0);

  // iterate
  int count = 0;
  while(fabs(diff) > params.tol){
    vfStep(params, K, Z, P, V0, V, G);
    diff = (V-V0).array().abs().maxCoeff();
    V0 = V;
    ++count;
    //cout << "Iteration: " << count << ", Max Value Function Diff: " << diff << endl;
  }

  // Compute solution time
  double toc = curr_second();
  double solTime  = toc - tic;

  // write to file (column major)
  ofstream fileSolTime, fileValue, filePolicy;
  fileValue.precision(10);
  filePolicy.precision(10);
  fileSolTime.open("solTimeCPP.dat");
  fileValue.open("valFunCPP.dat");
  filePolicy.open("polFunCPP.dat");
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
