classdef parameters<handle
    properties
        eta
        beta
        alpha
        delta
        mu
        rho
        sigma
        lambda
        nk
        nz
        tol
        maxtype
        howard
    end
    methods
        function load(obj, fileName)
            nParam = 13;
            params = textread(fileName, '%s', 'delimiter', ',');
            obj.eta = str2num(char(params(1)));
            obj.beta = str2num(char(params(3)));
            obj.alpha = str2num(char(params(5)));
            obj.delta = str2num(char(params(7)));
            obj.mu = str2num(char(params(9)));
            obj.rho = str2num(char(params(11)));
            obj.sigma = str2num(char(params(13)));
            obj.lambda = str2num(char(params(15)));
            obj.nk = str2num(char(params(17)));
            obj.nz = str2num(char(params(19)));
            obj.tol = str2num(char(params(21)));
            obj.maxtype = char(params(23));
            obj.howard = str2num(char(params(25)));
        end
    end
end
        