function [VError, GError] = solutionDiff(Method1, Method2)
    
    % Names of value and policy function files
    vDat = '/valueFunc.dat';
    pDat = '/policyFunc.dat';

    % Import first value function
    fileID = fopen([Method1 vDat]);
    V1 = textscan(fileID, '%f');
    fclose(fileID);

    % Import second value function
    fileID = fopen([Method2 vDat]);
    V2 = textscan(fileID, '%f');
    fclose(fileID);

    % Compute value function error
    VError = max(abs(V1{1}(3:end) - V2{1}(3:end)));

    % Import first policy function
    fileID = fopen([Method1 pDat]);
    G1 = textscan(fileID, '%f');
    fclose(fileID);

    % Import second policy function
    fileID = fopen([Method2 pDat]);
    G2 = textscan(fileID, '%f');
    fclose(fileID);

    % Compute policy function error
    GError = max(abs(G1{1}(3:end) - G2{1}(3:end)));
    
    % Write out
    fileOut = ['Errors_' Method1 '_' Method2 '.dat'];
    dlmwrite(fileOut, [VError; GError], '');

end
