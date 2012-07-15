/*============================================================================

 Description   Global variables for the value function iteration problem.

 =============================================================================

 Author:       Eric M. Aldrich

 Contact:      ealdrich@gmail.com

 Date:         28 July 2011

 ============================================================================*/
#include "globalvars.h"

// economic parameters
const int np = 6; // number of parameters
const REAL eta = 2; // risk aversion
const REAL beta = 0.984; // discount factor
const REAL alpha = 0.35; // capital share
const REAL delta = 0.01; // depreciation
const REAL mu = 0.0; // TFP mean
const REAL rho = 0.95; // TFP persistence
const REAL sigma = 0.005; // TFP volatility

// computational parameters
const int block_size = 4;
const int nk = 256; // this should be integer multiple of block_size
const int nz = 4;
const REAL tol = 0.00000001*(1-beta);

// maximization parameters
const char maxtype = 'g'; // maximization method - choices are 'g' (grid) and 'b' (binary search)
const int howard = 20; // maximize every how many steps? set howard = 1 if max = 'b'.

// to determine whether single or double precision is being used
//const float singletype;
//const double doubletype;
//const REAL realtype;
