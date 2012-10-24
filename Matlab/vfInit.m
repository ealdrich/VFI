%=============================================================================
%>
%> @file vfInit.m
%>
%> @brief File containing function to initialize the value function.
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
%> @brief Function to initialize value function.
%>
%> @details This function initializes the value function at the
%> deterministic steady state values for each level of TFP: conditional on
%> a TFP level, the deterministic steady-state value of capital is computed,
%> as well as the associated value function value.
%>
%> @param [in] param Object of class parameters.
%> @param [in] Z Grid of TFP values.
%>
%> @retval V Matrix of value function values.
%>
%=============================================================================
function V = vfInit(param, Z)

    % basic parameters
    nk = param.nk;
    nz = param.nz;
    alpha = param.alpha;
    beta = param.beta;
    delta = param.delta;
    eta = param.eta;

    % initialize
    Kj = ((1./(alpha*Z))*((1/beta)-1+delta)).^(1/(alpha-1));
	V(1,:) = ((Z.*(Kj.^alpha) - delta*Kj).^(1-eta))/(1-eta);
    V(2:nk,:) = V(ones(1,nk-1), :);
    
end
