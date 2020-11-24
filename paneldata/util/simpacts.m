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

    D_d = zeros(k,1);
    D_i = zeros(k,1);
    D_t = zeros(k,1);

    ind = 1;
    for jj = 1:k
        D_d(jj) = full(sqrt([SSSi*coefs_beta(jj)/N, SSi/N]*(N*hessi([1,jj+ind],[1,jj+ind]))*[SSSi*coefs_beta(jj)/N, SSi/N]'/N));
        D_t(jj) = full(sqrt([KKi*coefs_beta(jj)/N, Ki/N]*(N*hessi([1,jj+ind],[1,jj+ind]))*[KKi*coefs_beta(jj)/N, Ki/N]'/N));
        D_i(jj) = full(sqrt(([KKi*coefs_beta(jj)/N, Ki/N]-[SSSi*coefs_beta(jj)/N, SSi/N])*(N*hessi([1,jj+ind],[1,jj+ind]))*([KKi*coefs_beta(jj)/N, Ki/N]-[SSSi*coefs_beta(jj)/N, SSi/N])'/N));
    end
    
    % p-value
    pADI = (1 - normaldist(abs(ADI ./ D_d)))' * 2;
    pAII = (1 - normaldist(abs(AII ./ D_i)))' * 2;
    pATI = (1 - normaldist(abs(ATI ./ D_t)))' * 2;
    
    % Return structure with results
    simps.xnames = xnames;
    simps.k = k;
    
    simps.ADI = ADI;
    simps.AII = AII;
    simps.ATI = ATI;
    
    simps.seADI = D_d;
    simps.seAII = D_i;
    simps.seATI = D_t;
    
    simps.pADI = pADI';
    simps.pAII = pAII';
    simps.pATI = pATI';

end

