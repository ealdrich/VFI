//////////////////////////////////////////////////////////////////////////////
///
/// @file global.h
///
/// @brief Global header file.
///
//////////////////////////////////////////////////////////////////////////////

#ifndef __FILE_GLOBALVARS_H_SEEN__
#define __FILE_GLOBALVARS_H_SEEN__

#include <Eigen/Dense>

typedef double REAL;
typedef Eigen::Matrix<REAL, Eigen::Dynamic, 1> VectorXR;
typedef Eigen::Matrix<REAL, 1, Eigen::Dynamic> RowVectorXR;
typedef Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> MatrixXR;
typedef Eigen::Array<REAL, Eigen::Dynamic, 1> ArrayXR;
typedef Eigen::Array<REAL, Eigen::Dynamic, Eigen::Dynamic> ArrayXXR;

// Class for storing economic and computational parameters of the model
class parameters{
 public:
  REAL eta, beta, alpha, delta, mu, rho, sigma, lambda, tol;
  int nk, nz, howard;
  char maxtype;
  void load(const char*);
};

// Function declarations
REAL curr_second (void);
void ar1(const parameters& param, VectorXR& Z, MatrixXR& P);
void kGrid(const parameters& param, const VectorXR& Z, VectorXR& K);
void vfInit(const parameters& param, const VectorXR& Z, MatrixXR& V);
void vfStep(const parameters& param, const bool& howard, const VectorXR& K,
	    const VectorXR& Z, const MatrixXR& P, const MatrixXR& V0,
	    MatrixXR& V, Eigen::MatrixXi& G);
int binaryVal(const REAL& x, const int& nx, const VectorXR& X);
void gridMax(const int& klo, const int& nksub, const REAL& ydepK,
	     const REAL eta, const REAL beta, const VectorXR& K,
	     const VectorXR& Exp, REAL& V, int& G);
void binaryMax(const int& klo, const int& nksub, const REAL& ydepK,
	       const REAL eta, const REAL beta, const VectorXR& K,
	       const VectorXR& Exp, REAL& V, int& G);

#endif
