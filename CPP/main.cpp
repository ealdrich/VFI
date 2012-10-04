//////////////////////////////////////////////////////////////////////////////
///
/// @file main.cpp
///
/// @brief File containing main main function for the VFI problem.
///
//////////////////////////////////////////////////////////////////////////////

#include "global.h"
#include "auxFuncs.h"
#include <math.h>
#include <ctime>
#include <typeinfo>
#include <gsl/gsl_cblas.h>
#include <iostream>
#include <fstream>

using namespace std;

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
  REAL* K = new REAL[nk];
  REAL* Z = new REAL[nz];
  REAL* P = new REAL[(int)pow(nz,2)];
  REAL* V0 = new REAL[nk*nz];
  REAL* Vtemp;
  REAL* V = new REAL[nk*nz];
  REAL* G = new REAL[nk*nz];

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
    if(typeid(realtype) == typeid(singletype)){
      cblas_saxpy(nk*nz, -1.0, (float*)V, 1, (float*)V0, 1);
      imax = cblas_isamax(nk*nz, (float*)V0, 1);
    } else if(typeid(realtype) == typeid(doubletype)){
      cblas_daxpy(nk*nz, -1.0, (double*)V, 1, (double*)V0, 1);
      imax = cblas_idamax(nk*nz, (double*)V0, 1);
    }	
    diff = *(V0+imax);
    Vtemp = V0;
    V0 = V;
    V = Vtemp;
    ++count;
    //cout << "Iteration: " << count << ", Max Value Function Diff: " << diff << endl;
  }

  REAL toc = curr_second();
  REAL solTime  = toc - tic;
  cout << endl;
  cout << "Solution Time: " << solTime << endl;
  cout << endl;

  //V = V0; // this assignment doesn't work sometimes (it works in GPU code)
  // i resort to a full copy in lieu of pointer reassignemnt
  for(i = 0 ; i < nk ; ++i){
    for(j = 0 ; j < nz ; ++j){
      V[i*nz+j] = V0[i*nz+j];
    }
  }

  // write to file (column major)
  ofstream fileValue, filePolicy;
  fileValue.open("valueFunc.dat");
  filePolicy.open("policyFunc.dat");
  fileValue << nk << endl;
  fileValue << nz << endl;
  filePolicy << nk << endl;
  filePolicy << nz << endl;
  for(j = 0 ; j < nz ; ++j){
    for(i = 0 ; i < nk ; ++i){
      fileValue << V[i*nz+j] << endl;
      filePolicy << G[i*nz+j] << endl;
    }
  }  
  fileValue.close();
  filePolicy.close();

  return 0;

}
