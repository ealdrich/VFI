##############################################################################
#
# @file solutionDiff.R
#
# @brief File containing R code to compute differences between value and
# and policy functions for different implementations of the VFI problem.
#
# @details This code takes two command line arguments that correspond to
# differing software implementations used to solve a standard neoclassical
# growth model with value function iteration. Acceptable arguments are
# described in README.dox. The output consists of a file which lists two
# floating point values: the maximum absolute difference of their value and
# policy functions, respectively.
#
# @details See Aldrich, Eric M., Jesus Fernandez-Villaverde,
# A. Ronald Gallant and Juan F. Rubio-Ramirez (2011), "Tapping the
# supercomputer under your desk: Solving dynamic equilibrium models with
# graphics processors", Journal of Economic Dynamics & Control, 35, 386-393.
#
# @author Eric M. Aldrich \n
#         ealdrich@ucsc.edu
#
# @version 1.0
#
# @date 23 Oct 2012
#
# @copyright Copyright Eric M. Aldrich 2012 \n
#            Distributed under the Boost Software License, Version 1.0
#            (See accompanying file LICENSE_1_0.txt or copy at \n
#            http://www.boost.org/LICENSE_1_0.txt)
#
##############################################################################
  
# Determine which methods to compare from command line args
Methods = commandArgs(TRUE);
Method1 = Methods[1];
Method2 = Methods[2];
Dir1 = Method1;
Dir2 = Method2;

# Take care of ThrustGPU and ThrustOMP directory names
if(Method1 == 'ThrustGPU' | Method1 == 'ThrustOMP') Dir1 = 'Thrust';
if(Method2 == 'ThrustGPU' | Method2 == 'ThrustOMP') Dir2 = 'Thrust';

# Names of data files
vDat1 = paste('/valFun', Method1, '.dat', sep='');
vDat2 = paste('/valFun', Method2, '.dat', sep='');
pDat1 = paste('/polFun', Method1, '.dat', sep='');
pDat2 = paste('/polFun', Method2, '.dat', sep='');

# Import value functions and compute difference
V1 = scan(paste(Dir1, vDat1, sep=''), skip=2);
V2 = scan(paste(Dir2, vDat2, sep=''), skip=2);
VDiff = max(abs(V1 - V2));
  
# Import value functions and compute difference
G1 = scan(paste(Dir1, pDat1, sep=''), skip=2);
G2 = scan(paste(Dir2, pDat2, sep=''), skip=2);
GDiff = max(abs(G1 - G2));

# Write out
fileName = paste('Errors_', Method1, '_', Method2, '.dat', sep='');
write(c(VDiff, GDiff), fileName, ncolumns=1);

