%=============================================================================
%>
%> @file gridMax.m
%>
%> @brief File containing grid search maximization function.
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
%> @brief Function to compute maximum of Bellman objective via grid search.
%>
%> @details This function finds the maximum and argmax of the Bellman
%> objective function by using a naive grid search: computing the utility
%> at each value of the grid.
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
%>
%=============================================================================
function [V, G] = gridMax(klo, nksub, ydepK, eta, beta, K, Exp)

  w = ((ydepK-K(klo))^(1-eta))/(1-eta) + beta*Exp(1);
  wmax = w;
  windmax = 1;
  for(l = 2:nksub)
    w = ((ydepK-K(klo+l-1))^(1-eta))/(1-eta) + beta*Exp(l);
    if(w > wmax)
      wmax = w;
      windmax = l;
    end
  end
  V = wmax;
  G = klo+windmax-1;

end
