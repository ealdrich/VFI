function imax = binaryVal(x, X)
    if(x > X(end))
        imax = length(X);
    else
        imax = find(X > x, 1);
    end
end