%=============================================================================
%>
%> @file binaryVal.m
%>
%> @brief File containing a function which finds the approximate location of
%> a value in a vector with monotonically increasing values.
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
%> @brief Function to find the location of a value in a monotonic grid.
%>
%> @details This function finds the first value X[ix] such that x <= X[ix],
%> where x is a scalar value, X is a monotonic array, and ix is the index
%> of X.
%>
%> @param [in] x Value to search for in vector X.
%> @param [in] X Vector of data to search.
%>
%> @retval imax Integer ix (<= nx) such that x <= X[ix].
%>
%=============================================================================
function imax = binaryVal(x, X)
    if(x > X(end))
        imax = length(X);
    else
        imax = find(X > x, 1);
    end
end
