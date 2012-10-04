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

REAL curr_second (void);
void ar1(const REAL& lambda, VectorXR& Z, MatrixXR& P);
void kGrid(const VectorXR& Z, VectorXR& K);
void vfInit(const VectorXR& Z, MatrixXR& V);
void vfStep(const bool& howard, const VectorXR& K, const VectorXR& Z,
	    const MatrixXR& P, const MatrixXR& V0, MatrixXR& V,
	    Eigen::MatrixXi& G);
int binaryVal(const REAL& x, const int& nx, const VectorXR& X);
void gridMax(const int& klo, const int& nksub, const REAL& ydepK,
	     const VectorXR& K, const VectorXR& Exp, REAL& V, int& G);
void binaryMax(const int& klo, const int& nksub, const REAL& ydepK,
	       const VectorXR& K, const VectorXR& Exp, REAL& V, int& G);

// economic parameters
extern const REAL eta;
extern const REAL beta;
extern const REAL alpha;
extern const REAL delta;
extern const REAL mu;
extern const REAL rho;
extern const REAL sigma;

// computational parameters
extern const int nk;
extern const int nz;
extern const REAL tol;

// maximization parameters
extern const bool eigenMax;
extern const char maxtype;
extern const int howard;

// to determine whether single or double precision is being used
extern const float singletype;
extern const double doubletype;
extern const REAL realtype;

#endif
