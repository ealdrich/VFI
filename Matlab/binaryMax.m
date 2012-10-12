function [V,G] = binaryMax(beta, eta, klo, nksub, ydepK, K, Exp)

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
