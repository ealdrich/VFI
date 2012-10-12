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
                Exp = V0(klo:khi, :)*P(j,:)';
           
                % maximization
                if(matlabMax)
                    w = ((ydepK(i,j)*ones(nksub,1) - K(klo:khi)').^(1-eta))/(1-eta) + beta*Exp;
                    [V(i,j), G(i,j)] = max(w);
                    G(i,j) = G(i,j)+klo-1;
                else
                    if(maxtype == 'g')
                        [V(i,j), G(i,j)] = gridMax(beta, eta, klo, nksub, ydepK(i,j), K, Exp);
                    elseif(maxtype == 'b')
                        [V(i,j), G(i,j)] = binaryMax(beta, eta, klo, nksub, ydepK(i,j), K, Exp);
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
