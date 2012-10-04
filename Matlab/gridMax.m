function [V, G] = gridMax(beta, eta, klo, nksub, ydepK, K, Exp)

  w = ((ydepK-K(klo))^(1-eta))/(1-eta) + beta*Exp(1);
  windmax = 1;
  for(l = 2:nksub)
    w = ((ydepK-K(klo+l))^(1-eta))/(1-eta) + beta*Exp(l);
    if(w > wmax)
      wmax = w;
      windmax = l;
    end
  end
  V = wmax;
  G = klo+windmax;

end
