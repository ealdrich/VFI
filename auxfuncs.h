#ifndef __FILE_AUXFUNCS_H_SEEN__
#define __FILE_AUXSFUNCS_H_SEEN__

REAL ncdf(const REAL q);

int binary_val(const REAL& x, const int& nx, const REAL* X);

void ar1(const REAL& lambda, REAL* Z, REAL* P);

void kGrid(const REAL* Z, REAL* K);

void vfInit(const REAL* Z, REAL* V);

void grid_max(const int& klo, const int& nksub, int& l, REAL& w,
	      REAL& wmax, int& windmax, const REAL& ydepK,
	      const REAL* K, const REAL* Exp, REAL* V, REAL* G);

void binary_max(const int& klo, const int& nksub, int& kslo, int& kshi,
		int& ksmid1, int& ksmid2, REAL& w1, REAL& w2, REAL& w3,
		const REAL& ydepK, const REAL* K, const REAL* Exp,
		REAL* V, REAL* G);

void vfStep(int& klo, int& khi, int& nksub, int& kslo, int& kshi, int& ksmid1,
	    int& ksmid2, REAL& w, REAL& wmax, int& windmax, REAL& w1,
	    REAL& w2, REAL& w3, int& i, int& j, int& l, REAL& ydepK,
	    const bool& howard, const REAL* K,const REAL* Z,
	    const REAL* P, REAL* Exp, const REAL* V0, REAL* V,
	    REAL* G);

int vfiCPU(REAL* V, REAL* G);

int vfiGPU(REAL* hV, REAL* hG);

#endif
