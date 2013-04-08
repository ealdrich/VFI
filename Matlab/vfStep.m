%=============================================================================
%>
%> @file vfStep.m
%>
%> @brief File containing main iterative step of the VFI problem.
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
%> @brief Function to update value function.
%>
%> @details This function performs one iteration of the value function
%> iteration algorithm, using V0 as the current value function, maximizing
%> the LHS of the Bellman. Maximization is performed by @link binaryMax
%> @endlink.
%>
%> @param [in] param Object of class parameters.
%> @param [in] matlabMax Boolean which determines the Matlab max function is
%> used for maximization instead of @link gridMax @endlink or
%> @link binaryMax @endlink.
%> @param [in] K Grid of capital values.
%> @param [in] Z Grid of TFP values.
%> @param [in] P TFP transition matrix.
%> @param [in] V0 Matrix storing current value function.
%> @param [in] G0 Matrix storing current policy function.
%>
%> @retval V Matrix storing updated value function.
%> @retval G Matrix storing policy function.
%>
%=============================================================================
function [V,G] = vfStep(param, matlabMax, K, Z, P, V0, G0)

    % Basic parameters
    nk = param.nk;
    nz = param.nz;
    eta = param.eta;
    beta = param.beta;
    alpha = param.alpha;
    delta = param.delta;
    
    % output and depreciated capital
    ydepK = (K.^alpha)'*Z + (1-delta)*K(ones(1,nz),:)';

    for(i = 1:nk)
        for(j = 1:nz)
        
            % impose constraints on grid of future capital
            klo = 1;
            khi = binaryVal(ydepK(i,j), K); % consumption nonnegativity
	    if(K(khi) > ydepK(i,j))
                khi = khi - 1;
            end
            nksub = khi-klo+1;
        
            % continuation value for subgrid
	    % note that this computes more values than necessary for
	    % the maximization methods, but the Matlab matrix multiply
	    % is so efficient that it is faster to compute all possible
	    % continuation values outside of the max routine rather than
	    % only the necessary values inside the routine.
            Exp = V0(klo:khi, :)*P(j,:)';
           
            % maximization
            if(matlabMax)
                w = ((ydepK(i,j)*ones(nksub,1) - K(klo:khi)').^(1-eta))/(1-eta) + beta*Exp;
                [V(i,j), G(i,j)] = max(w);
                G(i,j) = G(i,j)+klo-1;
            else
	        [V(i,j), G(i,j)] = binaryMax(klo, nksub, ydepK(i,j), eta, beta, K, Exp);
            end
        end
    end
end
