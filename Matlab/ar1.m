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