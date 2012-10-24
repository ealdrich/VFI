%=============================================================================
%>
%> @file ar1.m
%>
%> @brief File containing AR1 function for the VFI problem.
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
%> @brief Function to compute discrete AR1 approximation values and
%> transition matrix.
%>
%> @details This function that computes a discrete AR1 approximation and
%> transition matrix using the method of Tauchen (1986).
%>
%> @param [in] param Object of class parameters.
%>
%> @retval Z Grid of AR1 values.
%> @retval P AR1 transition matrix values.
%>
%=============================================================================
function [Z, P] = ar1(param)

    % basic parameters
    nz = param.nz;
    mu = param.mu;
    rho = param.rho;
    sigma = param.sigma;
    lambda = param.lambda;

    % grid for TFP
    sigma_z = sigma/sqrt(1-rho^2);
    mu_z = mu/(1-rho);
    zmin = mu_z - lambda*sigma_z;
    zmax = mu_z + lambda*sigma_z;
    zstep = (zmax - zmin)/(nz-1);
    Z = exp(zmin:zstep:zmax);

    % transition matrix
    normarg1 = (zmin - mu - rho*log(Z))/sigma + 0.5*zstep/sigma;
    P(:,1) = 0.5 + 0.5*erf(normarg1/sqrt(2));
    for(j = 2:(nz-1))
        normarg1 = (log(Z(j)) - mu - rho*log(Z))/sigma + 0.5*zstep/sigma;
        normarg2 = (log(Z(j)) - mu - rho*log(Z))/sigma - 0.5*zstep/sigma;
        P(:,j) = 0.5*erf(normarg1/sqrt(2)) - 0.5*erf(normarg2/sqrt(2));
    end
    P(:,nz) = 1 - sum(P(:,1:(nz-1))')';
end
