% economic parameters
eta = 2; % Coefficient of relative risk aversion.
beta = 0.984; % Time discount factor.
alpha = 0.35; % Share of capital in the production function.
delta = 0.01; % Rate of capital depreciation.
mu = 0.0; % TFP mean.
rho = 0.95; % TFP persistence.
sigma = 0.005; % TFP volatility.

% computational parameters
nk = 64; % Number of values in capital grid.
nz = 4; % Number of values TFP grid.
tol = 1e-8*(1-beta); % Tolerance for convergence.

% maximization parameters
matlabMax = 0;
maxtype = 'b'; % Maximization method - choices are 'g'
               %(grid) and 'b' (binary search).
howard = 1; % Number of howard steps to perform between
            % maximizations - set howard = 1 if max = 'b'.

% admin
diff = 1.0;
tic;

% compute TFP grid, capital grid and initial VF
lambda = 3.0;
[Z, P] = ar1(rho, nz, mu, sigma, lambda);
K = kGrid(alpha, beta, delta, nk, Z);
V0 = vfInit(alpha, beta, delta, eta, nk, Z);

% iterate
count = 0;
how = false;
G0 = ones(nk,nz);
while(abs(diff) > tol)
    if(count <= 3 | mod(count, howard) == 0)
        how = false;
    else
        how = true;
    end
    [V, G] = vfStep(alpha, beta, delta, eta, matlabMax, maxtype, how, K, Z, P, V0,G0);
    diff = max(max(abs(V-V0)));
    V0 = V;
    G0 = G;
    count = count+1;
end

toc