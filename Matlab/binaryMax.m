%=============================================================================
%>
%> @file binaryMax.m
%>
%> @brief File containing binary search maximization function.
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
%> @brief Function to compute maximum of Bellman objective via binary search.
%>
%> @details This function finds the maximum and argmax of the Bellman
%> objective over a specified subgrid of capital by using a binary search
%> algorithm. The algorithm requires concavity and cannot be used with the
%> howard improvement method.
%>
%> @param [in] klo Lower index of the capital grid to begin search.
%> @param [in] nksub Number of points in the capital grid to include in
%> search.
%> @param [in] ydepK value of output plus depreciated capital.
%> @param [in] eta Coefficient of relative risk aversion.
%> @param [in] beta Time discount factor.
%> @param [in] K Grid of capital values.
%> @param [in] Exp Expected value function continuation values.
%>
%> @retval V Updated value function.
%> @retval G Updated policy function.
%>
%=============================================================================
function [V,G] = binaryMax(klo, nksub, ydepK, eta, beta, K, Exp)

    % binary search to find the vf max over K'
    % we assume that the value funtion is concave in capital
    kslo = 1;
    kshi = nksub;
  
    % case 1: capital grid has more than three values
    if(nksub > 3)
        % while the grid has 3 values or more, compute vf at midpoints
        % and revise the bounds of the grid
        while(kshi-kslo > 2)
            ksmid1 = floor((kslo + kshi)/2);
            ksmid2 = ksmid1+1;
            w1 = ((ydepK-K(klo+ksmid1-1))^(1-eta))/(1-eta) + beta*Exp(ksmid1);
            w2 = ((ydepK-K(klo+ksmid2-1))^(1-eta))/(1-eta) + beta*Exp(ksmid2);
            if(w2 > w1) kslo = ksmid1; else kshi = ksmid2; end
        end
        % when the grid is reduced to three values, find the max
        if(w2 > w1)
            w2 = ((ydepK-K(klo+kshi-1))^(1-eta))/(1-eta) + beta*Exp(kshi);
            if(w2 > w1)
                V = w2; G = klo+kslo;
            else
                V = w1; G = klo+kshi-1;
            end
        else
            w2 = ((ydepK-K(klo+kslo-1))^(1-eta))/(1-eta) + beta*Exp(kslo);
            if(w2 > w1)
                V = w2; G = klo+kslo-1;
            else
                V = w1; G = klo+kslo;
            end
        end
    
    % case 2: capital grid has three values
    elseif(nksub == 3)
        % evaluate vf at each value and determine max
        w1 = ((ydepK-K(klo+kslo-1))^(1-eta))/(1-eta) + beta*Exp(kslo);
        w2 = ((ydepK-K(klo+kslo))^(1-eta))/(1-eta) + beta*Exp(kslo+1);
        w3 = ((ydepK-K(klo+kshi-1))^(1-eta))/(1-eta) + beta*Exp(kshi);
        V = w1;
        G = klo+kslo-1;
        if(w2 > V)
            V = w2; G = klo+kslo;
        end
        if(w3 > V)
            V = w3; G = klo+kshi-1;
        end
    
    % case 3: capital grid has one or two values
    else
        % evaluate vf at each value and determine max
        w1 = ((ydepK-K(klo+kslo-1))^(1-eta))/(1-eta) + beta*Exp(kslo);
        w2 = ((ydepK-K(klo+kshi-1))^(1-eta))/(1-eta) + beta*Exp(kshi);
        if(w2 > w1)
            V = w2; G = klo+kshi-1;
        else
            V = w1; G = klo+kslo-1;
        end
    end
end
