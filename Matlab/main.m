%=============================================================================
%>
%> @file main.m
%>
%> @brief File containing Matlab code to solve the VFI problem.
%>
%> @details This code solves a standard neoclassical growth model with value
%> function iteration on a CPU.
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

format long;

% load parameters
params = parameters;
load(params, '../parameters.txt')
matlabMax = 0;

% admin
diff = 1.0;
tic;

% compute TFP grid, capital grid and initial VF
lambda = 3.0;
[Z, P] = ar1(params);
K = kGrid(params, Z);
V0 = vfInit(params, Z);

% iterate
count = 0;
how = false;
G0 = ones(params.nk,params.nz);
while(abs(diff) > params.tol)
    if(count <= 3 | mod(count, params.howard) == 0)
        how = false;
    else
        how = true;
    end
    [V, G] = vfStep(params, matlabMax, how, K, Z, P, V0,G0);
    diff = max(max(abs(V-V0)));
    V0 = V;
    G0 = G;
    count = count+1;
end

% compute solution time
solTime = toc;

% write to file
dlmwrite('solTimeMatlab.dat', solTime, '');
dlmwrite('valFunMatlab.dat', [params.nk; params.nz; V(:)], '');
dlmwrite('polFunMatlab.dat', [params.nk; params.nz; G(:)], '');
