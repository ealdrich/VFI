#include "global.h"
#include <fstream>

using namespace std;

void parameters::load(char* fileName)
{
  int nParam = 13;
  ifstream fileIn;
  fileIn.open(fileName);
  std::vector<std::string> parameters(nParam);
  std::string trash;
  for(int ix = 0 ; ix < nParam ; ++ix){
    getline(fileIn, parameters[ix], ',');
    getline(fileIn, trash);
  }
  eta = atof(parameters[0].c_str());
  beta = atof(parameters[1].c_str());
  alpha = atof(parameters[2].c_str());
  delta = atof(parameters[3].c_str());
  mu = atof(parameters[4].c_str());
  rho = atof(parameters[5].c_str());
  sigma = atof(parameters[6].c_str());
  lambda = atoi(parameters[7].c_str());
  nk = atoi(parameters[8].c_str());
  nz = atoi(parameters[9].c_str());
  tol = atof(parameters[10].c_str());
  maxtype = parameters[11][0];
  howard = atoi(parameters[12].c_str());
}
