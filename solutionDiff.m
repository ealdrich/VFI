%=============================================================================
%>
%> @file solutionDiff.R
%>
%> @brief File containing R code to compute differences between value and
%> and policy functions for different implementations of the VFI problem.
%>
%> @details This code takes two command line arguments that correspond to
%> differing software implementations used to solve a standard neoclassical
%> growth model with value function iteration. Acceptable arguments are
%> described in README.dox. The output consists of a file which lists two
%> floating point values: the maximum absolute difference of their value and
%> policy functions, respectively.
%>
%> @details See Aldrich, Eric M., Jesus Fernandez-Villaverde,
%> A. Ronald Gallant and Juan F. Rubio-Ramirez (2011), "Tapping the
%> supercomputer under your desk: Solving dynamic equilibrium models with
%> graphics processors", Journal of Economic Dynamics & Control, 35, 386-393.
%>
%> @author Eric M. Aldrich \n
%>         ealdrich@ucsc.edu
%>
%> @version 1.0
%>
%> @date 23 Oct 2012
%>
%> @copyright Copyright Eric M. Aldrich 2012 \n
%>            Distributed under the Boost Software License, Version 1.0
%>            (See accompanying file LICENSE_1_0.txt or copy at \n
%>            http://www.boost.org/LICENSE_1_0.txt)
%>
%=============================================================================

function solutionDiff(Method1, Method2)
    
    % Names of value and policy function files
    vDat1 = ['/valFun' Method1 '.dat'];
    vDat2 = ['/valFun' Method2 '.dat'];
    pDat1 = ['/polFun' Method1 '.dat'];
    pDat2 = ['/polFun' Method2 '.dat'];

    % Take care of ThrustGPU and ThrustOMP directory names
    Dir1 = Method1;
    Dir2 = Method2;
    if strcmp(Method1, 'ThrustGPU') | strcmp(Method1, 'ThrustOMP')
      Dir1 = 'Thrust';
    end
    if strcmp(Method2, 'ThrustGPU') | strcmp(Method2, 'ThrustOMP')
      Dir2 = 'Thrust';
    end

    % Import first value function
    fileID = fopen([Dir1 vDat1]);
    V1 = textscan(fileID, '%f');
    fclose(fileID);

    % Import second value function
    fileID = fopen([Dir2 vDat2]);
    V2 = textscan(fileID, '%f');
    fclose(fileID);

    % Compute value function differences
    VDiff = max(abs(V1{1}(3:end) - V2{1}(3:end)));

    % Import first policy function
    fileID = fopen([Dir1 pDat1]);
    G1 = textscan(fileID, '%f');
    fclose(fileID);

    % Import second policy function
    fileID = fopen([Dir2 pDat2]);
    G2 = textscan(fileID, '%f');
    fclose(fileID);

    % Compute policy function differences
    GDiff = max(abs(G1{1}(3:end) - G2{1}(3:end)));
    
    % Write out
    fileOut = ['Errors_' Method1 '_' Method2 '.dat'];
    dlmwrite(fileOut, [VDiff; GDiff], '');

end
