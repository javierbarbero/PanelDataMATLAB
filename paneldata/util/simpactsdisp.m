function simpactsdisp( estsimpacts )
%SIMPACTSDISP Display computed spatial impacts
%   Display computed spatial impacts
%
%   SIMPACTSDISP( est ) display computed spatial impacts.
% 
%   Example
%     
%      est  = spanel(id, time, y, X, W, 'fe', 'slagy', 1);
%      simps = simpacts(est);
%      simpactsdisp(simps)
%
%   See also S2SLS, SPANEL, SIMPACTS 
%
%   Copyright 2013-2020 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 24, November, 2020

    xnames = estsimpacts.xnames;
    k = estsimpacts.k;
    
    ADI = estsimpacts.ADI;
    AII = estsimpacts.AII;
    ATI = estsimpacts.ATI;
    
    seADI = estsimpacts.seADI;
    seAII = estsimpacts.seAII;
    seATI = estsimpacts.seATI;   
    
    pADI = estsimpacts.pADI;
    pAII = estsimpacts.pAII;
    pATI = estsimpacts.pATI;

    % Code adapted from: Arbia, Bera, Dogan, and Taspinar (2020)

    % Print table
    fprintf('\n');
    fprintf('-------------------------------------------\n');
    fprintf('%12s%10s%10s%10s \n','Variable','Direct','Indirect','Total');
    fprintf('-------------------------------------------\n');
    for jj = 1:k
        fprintf('%12s%10.6f%10.6f%10.6f \n',string(xnames(jj)),ADI(jj),AII(jj),ATI(jj));
        fprintf('%12s%10.6f%10.6f%10.6f \n','Std. Error',seADI(jj),seAII(jj),seATI(jj));
        fprintf('%12s%10.6f%10.6f%10.6f \n','p-value',pADI(jj),pAII(jj),pATI(jj));
        fprintf('-------------------------------------------\n');
    end

end

