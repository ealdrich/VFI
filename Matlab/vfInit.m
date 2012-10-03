function V = vfInit(alpha, beta, delta, eta, nk, Z)

    % initialize
    Kj = ((1./(alpha*Z))*((1/beta)-1+delta)).^(1/(alpha-1));
	V(1,:) = ((Z.*(Kj.^alpha) - delta*Kj).^(1-eta))/(1-eta);
    V(2:nk,:) = V(ones(1,nk-1), :);
    
end