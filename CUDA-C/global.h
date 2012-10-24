//////////////////////////////////////////////////////////////////////////////
///
/// @file global.h
///
/// @brief Global header file.
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

#ifndef __FILE_GLOBALVARS_H_SEEN__
#define __FILE_GLOBALVARS_H_SEEN__

typedef double REAL;

//////////////////////////////////////////////////////////////////////////////
///
/// @class parameters
///
/// @brief Object to store parameter values for VFI problem.
///
//////////////////////////////////////////////////////////////////////////////
class parameters{
 public:
  REAL eta; ///< Coefficient of relative risk aversion.
  REAL beta; ///< Time discount factor.
  REAL alpha; ///< Share of capital in the production function.
  REAL delta; ///< Rate of capital depreciation.
  REAL mu; ///< TFP mean.
  REAL rho; ///< TFP persistence.
  REAL sigma; ///< TFP volatility.
  REAL lambda; ///< Number of standard deviations for AR1 approximation.
  int nk; ///< Number of values in capital grid.
  int nz; ///< Number of values in TFP grid.
  REAL tol; ///< Tolerance for convergence.
  char maxtype; ///< @brief Maximization method - choices are `g' (grid) and `b' (binary search).
  int howard; ///< @brief Number of howard steps to perform between maximizations - set howard = 1 if max = `b'.
  void load(const char*);
};

// Function declarations
double curr_second (void);
void ar1(const parameters& param, REAL* Z, REAL* P);
void kGrid(const parameters& param, const REAL* Z, REAL* K);
void vfInit(const parameters& param, const REAL* Z, REAL* V);

// to determine whether single or double precision is being used
extern const float singletype;
extern const double doubletype;
extern const REAL realtype;

#endif
