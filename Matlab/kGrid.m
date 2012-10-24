%=============================================================================
%>
%> @file kGrid.m
%>
%> @brief File containing function to create capital grid.
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

%=============================================================================
%>
%> @brief Function to compute the values of an equally spaced capital grid.
%>
%> @details This function computes an equally spaced capital grid. The
%> upper and lower bounds are the deterministic steady-state values of
%> capital at the highest and lowest values of the TFP process
%> (respectively), scaled by 0.95 and 1.05 (respectively).
%>
%> @param [in] param Object of class parameters.
%> @param [in] Z Grid of TFP values.
%>
%> @retval K Grid of capital values.
%>
%=============================================================================
function K = kGrid(param, Z)

    % basic parameters
    nk = param.nk;
    nz = param.nz;
    alpha = param.alpha;
    beta = param.beta;
    delta = param.delta;
    eta = param.eta;

    % initial grid for capital
    kmin = 0.95*(((1/(alpha*Z(1)))*((1/beta)-1+delta))^(1/(alpha-1)));
    kmax = 1.05*(((1/(alpha*Z(end)))*((1/beta)-1+delta))^(1/(alpha-1)));
    kstep = (kmax - kmin)/(nk-1);
    K = kmin:kstep:kmax;

end
