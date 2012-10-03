function K = kGrid(alpha, beta, delta, nk, Z)

    % initial grid for capital
    kmin = 0.95*(((1/(alpha*Z(1)))*((1/beta)-1+delta))^(1/(alpha-1)));
    kmax = 1.05*(((1/(alpha*Z(end)))*((1/beta)-1+delta))^(1/(alpha-1)));
    kstep = (kmax - kmin)/(nk-1);
    K = kmin:kstep:kmax;

end
