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
%> iteration algorithm, using V0 as the current value function and either
%> maximizing the LHS of the Bellman if @link howard @endlink = false or
%> using the concurrent policy function as the argmax if
%> @link howard @endlink = true. Maximization is performed by either
%> @link gridMax @endlink or @link binaryMax @endlink.
%>
%> @param [in] param Object of class parameters.
%> @param [in] matlabMax Boolean which determines the Matlab max function is
%> used for maximization instead of @link gridMax @endlink or
%> @link binaryMax @endlink.
%> @param [in] howard Indicates if the current iteration of the value
%> function will perform a maximization (false) or if it will simply compute
%> the new value function using the concurrent policy function (true).
%> @param [in] K Grid of capital values.
%> @param [in] Z Grid of TFP values.
%> @param [in] P TFP transition matrix.
%> @param [in] V0 Matrix storing current value function.
%> @param [in] G0 Matrix storing current policy function.
%>
%> @retval V Matrix storing updated value function.
%> @retval G Matrix storing policy function (updated if howard = false).
%>
%=============================================================================
function [V,G] = vfStep(param, matlabMax, howard, K, Z, P, V0, G0)

    % Basic parameters
    nk = param.nk;
    nz = param.nz;
    eta = param.eta;
    beta = param.beta;
    alpha = param.alpha;
    delta = param.delta;
    maxtype = param.maxtype;
    
    % output and depreciated capital
    ydepK = (K.^alpha)'*Z + (1-delta)*K(ones(1,nz),:)';

    for(i = 1:nk)
        for(j = 1:nz)
        
            % maximize on non-howard steps
            if(howard == false)
          
                % impose constraints on grid of future capital
                klo = 1;
                khi = binaryVal(ydepK(i,j), K); % consumption nonnegativity
                if(K(khi) > ydepK(i,j))
                    khi = khi - 1;
                end
                
                % further restrict capital grid via monotonicity (CPU only)
                if(i > 1)
                    if(G(i-1,j) > klo & G(i-1,j) < khi)
                        klo = G(i-1,j);
                    end
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
                    if(maxtype == 'g')
		      [V(i,j), G(i,j)] = gridMax(klo, nksub, ydepK(i,j), eta, beta, K, Exp);
                    elseif(maxtype == 'b')
		      [V(i,j), G(i,j)] = binaryMax(klo, nksub, ydepK(i,j), eta, beta, K, Exp);
                    end
                end
            
            else
                Exp = V0(G0(i,j),:)*P(j,:)';
                V(i,j) = ((ydepK(i,j)-K(G0(i,j)))^(1-eta))/(1-eta) + beta*Exp(1);
            end
        
        end
    end
    if(howard == true)
        G = G0;
    end
end
