/*============================================================================

 Function      vfiGPU

 Usage         vfiGPU(hV, hG)

 Arguments     hV: pointer to array of REALs storing the value function.
                   
               hG: pointer to array of REALs storing the policy function.
	              
 Description   This function performs value function iteration on the GPU,
               finding the maximum of the Bellman objective function for each
	       node in the state space and iterating until convergence.

 Dependencies  Global Variables: eta, beta, alpha, delta, mu, rho, sigma,
                                 block_size, nk, nz, tol, maxtype, howard
				 (globalvas.h).

               Functions:        pow (math.h);
	                         cblas(S,D)axpy, cblasI(s,d)amax (cblas.h).

	       Kernels:          ar1GPU (ar1.cu), kGridGPU (kGrid.cu),
	                         vfInitGPU(vfInit.cu), vfStepGPU (vfStep.cu).

 Return value  Returns 0 upon successful completion, 1 otherwise.

 =============================================================================

 Author:       Eric M. Aldrich

 Contact:      ealdrich@gmail.com

 Date:         28 July 2011

 ============================================================================*/

#include "global.h"
#include "auxFuncs.h"
#include "cublas_v2.h"
#include <iostream>
#include <ctime>
#include <typeinfo>

using namespace std;

#include "ar1GPU.cu"
#include "kGridGPU.cu"
#include "vfInitGPU.cu"
#include "vfStepGPU.cu"

int vfiGPU(REAL* hV, REAL* hG) 
{ 

  int imax;
  REAL diff = 1.0;
  cublasHandle_t handle;
  cublasStatus_t status;
  status = cublasCreate(&handle);
  REAL negOne = -1.0;

  // pointers to variables in device memory
  REAL* K;
  REAL* Z;
  REAL* P;
  REAL* V0;
  REAL* V;
  REAL* G;
  REAL* Vtemp;

  // allocate variables in device memory
  size_t sizeK = nk*sizeof(REAL);
  size_t sizeZ = nz*sizeof(REAL);
  size_t sizeP = nz*nz*sizeof(REAL);
  size_t sizeV = nk*nz*sizeof(REAL);
  size_t sizeG = nk*nz*sizeof(REAL);
  clock_t start = clock();
  cudaMalloc((void**)&K, sizeK);
  cudaMalloc((void**)&Z, sizeZ);
  cudaMalloc((void**)&P, sizeP);
  cudaMalloc((void**)&V0, sizeV);
  cudaMalloc((void**)&Vtemp, sizeV);
  cudaMalloc((void**)&V, sizeV);
  cudaMalloc((void**)&G, sizeG);
  cout << "GPU Memory Allocation: " << (clock() - start)/(REAL)CLOCKS_PER_SEC << endl;

  // blocking
  dim3 dimBlockZ(nz, 1);
  dim3 dimBlockK(block_size,1);
  dim3 dimBlockV(block_size, nz);
  dim3 dimGridZ(1,1);
  dim3 dimGridK(nk/block_size,1);
  dim3 dimGridV(nk/block_size,1);

  // compute TFP grid, capital grid and initial VF
  REAL lambda = 3;
  ar1GPU<<<dimGridZ,dimBlockZ>>>(nz,lambda,mu,sigma,rho,Z,P);
  kGridGPU<<<dimGridK,dimBlockK>>>(nk,nz,beta,alpha,delta,Z,K);
  vfInitGPU<<<dimGridV,dimBlockV>>>(nz,eta,beta,alpha,delta,Z,V0);

  // iterate on the value function
  int count = 0;
  bool how = false;
  start = clock();
  while(fabs(diff) > tol){
    if(count < 3 | count % howard == 0) how = false; else how = true;
    vfStepGPU<<<dimGridV,dimBlockV>>>(nk,nz,eta,beta,alpha,delta,maxtype,how,K,Z,P,V0,V,G);
    if(typeid(realtype) == typeid(singletype)){
      status = cublasSaxpy(handle, nk*nz, (float*)&negOne, (float*)V, 1, (float*)V0, 1);
      status = cublasIsamax(handle, nk*nz, (float*)V0, 1, &imax);
    } else if(typeid(realtype) == typeid(doubletype)){
      status = cublasDaxpy(handle, nk*nz, (double*)&negOne, (double*)V, 1, (double*)V0, 1);
      status = cublasIdamax(handle, nk*nz, (double*)V0, 1, &imax);
    }
    cudaMemcpy(&diff, V0+imax, sizeof(REAL), cudaMemcpyDeviceToHost);
    Vtemp = V0;
    V0 = V;
    V = Vtemp;
    ++count;
  }
  cout << "GPU Solve Time: " << (clock() - start)/(REAL)CLOCKS_PER_SEC << endl;
  V = V0;
  
  // copy value and policy functions to host memory
  cudaMemcpy(hV, V, sizeV, cudaMemcpyDeviceToHost);
  cudaMemcpy(hG, G, sizeG, cudaMemcpyDeviceToHost);

  // copy state variable grids and transition matrix to host memory
  REAL* hK = new REAL[nk];
  REAL* hZ = new REAL[nz];
  REAL* hP = new REAL[nz*nz];
  cudaMemcpy(hK, K, sizeK, cudaMemcpyDeviceToHost);
  cudaMemcpy(hZ, Z, sizeZ, cudaMemcpyDeviceToHost);
  cudaMemcpy(hP, P, sizeP, cudaMemcpyDeviceToHost);

  // free variables in device memory
  cudaFree(K);
  cudaFree(Z);
  cudaFree(P);
  cudaFree(V0);
  cudaFree(V);
  cudaFree(Vtemp);
  cudaFree(G);
  cublasDestroy(handle);

  printMatrix<REAL>(0, nk, nz, &hV[0], 4, nz, 10);

  return 0;

}
