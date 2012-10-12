% load parameters
params = parameters;
load(params, '../parameters.txt')
matlabMax = 0;

% admin
diff = 1.0;
tic;

% compute TFP grid, capital grid and initial VF
lambda = 3.0;
[Z, P] = ar1(params);
K = kGrid(params, Z);
V0 = vfInit(params, Z);

% iterate
count = 0;
how = false;
G0 = ones(params.nk,params.nz);
while(abs(diff) > params.tol)
    if(count <= 3 | mod(count, params.howard) == 0)
        how = false;
    else
        how = true;
    end
    [V, G] = vfStep(params, matlabMax, how, K, Z, P, V0,G0);
    diff = max(max(abs(V-V0)));
    V0 = V;
    G0 = G;
    count = count+1;
end

% compute solution time
solTime = toc;

% write to file
dlmwrite('solutionTime.dat', solTime, '');
dlmwrite('valueFunc.dat', [params.nk; params.nz; V(:)], '');
dlmwrite('policyFunc.dat', [params.nk; params.nz; G(:)], '');
