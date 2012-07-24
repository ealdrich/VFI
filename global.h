#ifndef __FILE_GLOBALVARS_H_SEEN__
#define __FILE_GLOBALVARS_H_SEEN__

typedef double REAL;

// economic parameters
extern const int np;
extern const REAL eta;
extern const REAL beta;
extern const REAL alpha;
extern const REAL delta;
extern const REAL mu;
extern const REAL rho;
extern const REAL sigma;

// computational parameters
extern const int block_size;
extern const int nk;
extern const int nz;
extern const REAL tol;

// maximization parameters
extern const char maxtype;
extern const int howard;

// to determine whether single or double precision is being used
extern const float singletype;
extern const double doubletype;
extern const REAL realtype;

#endif
