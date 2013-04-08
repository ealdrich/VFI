//////////////////////////////////////////////////////////////////////////////
///
/// @file parameters.cpp
///
/// @brief File containing parameters class method for loading VFI
/// parameter values.
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

#include "global.h"
#include <stdlib.h>
#include <vector>
#include <fstream>

using namespace std;

//////////////////////////////////////////////////////////////////////////////
///
/// @brief Function to load VFI parameter values to parameters object.
///
/// @details This function is a parameters class method which loads
/// parameter values from a text file for storage in the object. The input
/// file must have 13 lines, each line beginning with a parameter value,
/// followed by a comma and a character string describing the parameter. The
/// order of the parameters must correspond to the order in the parameters
/// class description.
///
/// @param [in] fileName Name of file storing parameter values.
///
/// @returns Void.
///
//////////////////////////////////////////////////////////////////////////////
void parameters::load(const char* fileName)
{
  int nParam = 11;
  ifstream fileIn;
  fileIn.open(fileName);
  std::vector<std::string> params(nParam);
  std::string trash;
  for(int ix = 0 ; ix < nParam ; ++ix){
    getline(fileIn, params[ix], ',');
    getline(fileIn, trash);
  }
  eta = atof(params[0].c_str());
  beta = atof(params[1].c_str());
  alpha = atof(params[2].c_str());
  delta = atof(params[3].c_str());
  mu = atof(params[4].c_str());
  rho = atof(params[5].c_str());
  sigma = atof(params[6].c_str());
  lambda = atof(params[7].c_str());
  nk = atoi(params[8].c_str());
  nz = atoi(params[9].c_str());
  tol = atof(params[10].c_str());
}
