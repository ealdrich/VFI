#include "global.h"
#include <fstream>
#include <vector>

using namespace std;

void parameters::load(const char* fileName)
{
  int nParam = 13;
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
  maxtype = params[11][0];
  howard = atoi(params[12].c_str());
}
