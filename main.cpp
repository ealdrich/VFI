/*============================================================================

 Function      main

 Description   Solves a standard stochastic ramsey model via value function
               iteration on a CPU and/or GPU. The flags 'gpusolve' and
	       'cpusolve' allow the user to solve the the problem on only the
	       CPU, the GPU or both. The maximization step at each  iteration
	       of the algorithm is carried out by a grid or binary search
	       method. The variable 'max' specifies the type of maximation
	       method to use, and is defined in 'global.cpp'. Economic
	       variables and grid sizes are also defined in 'global.cpp'.

 Dependencies  Global variables: nk, nz (global.h).

               Functions:        clock() (ctime.h);
	                         cblas_(s,d)axpy, cblas_i(s,d)amax (cblas.h);
				 vfiCPU, vfiGPU (auxfuncs.h).

 Return value  Returns 0 upon successful completion, 1 otherwise.

 =============================================================================

 Author:       Eric M. Aldrich

 Contact:      ealdrich@gmail.com

 Date:         28 July 2011

 ============================================================================*/

#include "global.h"
#include "auxfuncs.h"
#include <iostream>
#include <ctime>
#include <typeinfo>
#include <gsl/gsl_cblas.h>

using namespace std;

int main()
{

  // admin
  int i,j;
  clock_t start;
  const bool gpusolve = true;
  const bool cpusolve = false;

  // declare value (V) and policy (G) function matrices
  REAL* Vgpu = new REAL[nk*nz];
  REAL* Ggpu = new REAL[nk*nz];
  REAL* Vcpu = new REAL[nk*nz];
  REAL* Gcpu = new REAL[nk*nz];

  // vfi on gpu
  if(gpusolve){
    start = clock();
    vfiGPU(Vgpu, Ggpu);
    cout << "GPU Total Time: " << (clock() - start)/(REAL)CLOCKS_PER_SEC << endl;
  }

  // vfi on cpu
  if(cpusolve){
    start = clock();
    vfiCPU(Vcpu, Gcpu);
    cout << "CPU Solve Time: " << (clock() - start)/(REAL)CLOCKS_PER_SEC << endl;
  }

  // report max abs diff between GPU and CPU value functions and policy functions
  if(gpusolve & cpusolve){
    int maxVdiff, maxGdiff;
    if(typeid(realtype) == typeid(singletype)){
      cblas_saxpy(nk*nz, -1.0, (float*)Vcpu, 1, (float*)Vgpu, 1);
      cblas_saxpy(nk*nz, -1.0, (float*)Gcpu, 1, (float*)Ggpu, 1);
      maxVdiff = cblas_isamax(nk*nz, (float*)Vgpu, 1);
      maxGdiff = cblas_isamax(nk*nz, (float*)Ggpu, 1);
    } else if(typeid(realtype) == typeid(doubletype)){
      cblas_daxpy(nk*nz, -1.0, (double*)Vcpu, 1, (double*)Vgpu, 1);
      cblas_daxpy(nk*nz, -1.0, (double*)Gcpu, 1, (double*)Ggpu, 1);
      maxVdiff = cblas_idamax(nk*nz, (double*)Vgpu, 1);
      maxGdiff = cblas_idamax(nk*nz, (double*)Ggpu, 1);
    }
    cout << "Max Absolute Difference for CPU and GPU Value Functions : " << Vgpu[maxVdiff] << endl;
    cout << "Max Absolute Difference for CPU and GPU Policy Functions : " << Ggpu[maxGdiff] << endl;
  }  

  return 0;
}
