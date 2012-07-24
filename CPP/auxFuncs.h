//////////////////////////////////////////////////////////////////////////////
///
/// @file auxFuncs.h
///
/// @brief Function declarations.
///
/// @author Eric M. Aldrich \n
///         ealdrich@ucsc.edu
///
/// @version 1.0
///
/// @date 18 July 2012
///
/// @copyright Copyright Eric M. Aldrich 2012 \n
///            Distributed under the Boost Software License, Version 1.0
///            (See accompanying file LICENSE_1_0.txt or copy at \n
///            http://www.boost.org/LICENSE_1_0.txt)
///
//////////////////////////////////////////////////////////////////////////////

#ifndef __FILE_AUXFUNCS_H_SEEN__
#define __FILE_AUXSFUNCS_H_SEEN__

int binaryVal(const REAL& x, const int& nx, const REAL* X);

void ar1(const REAL& lambda, REAL* Z, REAL* P);

void kGrid(const REAL* Z, REAL* K);

void vfInit(const REAL* Z, REAL* V);

void gridMax(const int& klo, const int& nksub, const REAL& ydepK,
		const REAL* K, const REAL* Exp, REAL* V, REAL* G);

void binaryMax(const int& klo, const int& nksub, const REAL& ydepK,
		  const REAL* K, const REAL* Exp, REAL* V, REAL* G);

void vfStep(const bool& howard, const REAL* K, const REAL* Z,
	       const REAL* P, const REAL* V0, REAL* V, REAL* G);

#endif
