function [ simps ] = simpacts( est )
%SIMPACTS Computes spatial impacts of estimated model
%   Computes spatial impacts of estimated model
%
%   est = SIMPACTS( est ) computes the spatial impacts of the estimated
%   spatial mode.
%
%   Standard errors are computed using the Delta Method following Arbia,
%   Bera, Dogan, and Taspinar (2020), "Testing Impact Measures in Spatial
%   Autoregressive Models". International Refional Science Review. Vol
%   (43). doi: 10.1177/0160017619826264
% 
%   Example
%     
%      est  = spanel(id, time, y, X, W, 'fe', 'slagy', 1);
%      simps = simpacts(est);
%      simpactsdisp(simps)
%
%   See also S2SLS, SPANEL SIMPACTSDISP
%
%   Copyright 2013-2020 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 24, November, 2020
%

    % Check that model has spatial lag of the dependent variable
    if est.slagy ~= 1
        error("Model does not have spatial lag of dependent variable.")
    end
    
    % Error if additional spatial lags of X variables
    if ~isempty(est.slagX)
        error("Not yet implemented when additional spatial lags of X are included in the model.")
    end

    % Get W matrix
    if est.isPanel
        W = kron(sparse(est.W), eye(est.T));
    else
        W = est.W;
    end

    % Code adapted from: Arbia, Bera, Dogan, and Taspinar (2020)
    
    % Get date from estimation    
    k = est.k - 1 - est.hasConstant;
    coefs_beta = est.coef(1:k);
    coefs_stderr = est.stderr(1:k);
    lambda = est.coef(k+1); % Spatial lag estimated coefficient
    lambda_stderr = est.stderr(k+1); % Spatial lag estimated standard error
    N = est.N; % Total number of observations     
    xnames = est.xnames;
    
    % Compute spatial impacts impacts
    Si = (speye(N) - lambda * W)\speye(N);
    SSi = trace(Si);
    Ki = sum(sum(Si));

    ADI = full(SSi * coefs_beta/N);
    ATI = full(Ki * coefs_beta/N);
    AII = full(ATI - ADI);

    SSSi = trace(Si * W * Si);
    KKi = sum(sum(Si * W * Si));
    
    % Delta Method
    hessi = est.varcoef;
    hessi = circshift(hessi,[1,1]); % Put spatial lag in first place

    seADI = zeros(k,1);
    seAII = zeros(k,1);
    seATI = zeros(k,1);

    ind = 1;
    for jj = 1:k
        seADI(jj) = full(sqrt([SSSi*coefs_beta(jj)/N, SSi/N]*(N*hessi([1,jj+ind],[1,jj+ind]))*[SSSi*coefs_beta(jj)/N, SSi/N]'/N));
        seATI(jj) = full(sqrt([KKi*coefs_beta(jj)/N, Ki/N]*(N*hessi([1,jj+ind],[1,jj+ind]))*[KKi*coefs_beta(jj)/N, Ki/N]'/N));
        seAII(jj) = full(sqrt(([KKi*coefs_beta(jj)/N, Ki/N]-[SSSi*coefs_beta(jj)/N, SSi/N])*(N*hessi([1,jj+ind],[1,jj+ind]))*([KKi*coefs_beta(jj)/N, Ki/N]-[SSSi*coefs_beta(jj)/N, SSi/N])'/N));
    end
    
    % p-value
    pADI = (1 - normaldist(abs(ADI ./ seADI)))' * 2;
    pAII = (1 - normaldist(abs(AII ./ seAII)))' * 2;
    pATI = (1 - normaldist(abs(ATI ./ seATI)))' * 2;
    
    % Return structure with results
    simps.xnames = xnames;
    simps.k = k;
    
    simps.ADI = ADI;
    simps.AII = AII;
    simps.ATI = ATI;
    
    simps.seADI = seADI;
    simps.seAII = seAII;
    simps.seATI = seATI;
    
    simps.pADI = pADI';
    simps.pAII = pAII';
    simps.pATI = pATI';

end

