function [ test ] = bsjksatest( est )
%BSJKSATEST BSJK test of individual effects
%   Computes the Baltagi, Song, Jung and Koh's test of spatial
%   autocorrelation, serial correlation and random effects.
%
%   testo = BSJKSATEST( est ) Computes the BSJK test for the specified 
%   estimation output. Returns a test output structure, testout.
%
%   Example
%     
%      test = bsjksatest(est);
%
%   See also TESTOUT
%
%   Copyright 2013-2015 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 17, October, 2015
%

    if est.isMultiEq
        error('Baltagi, Song, Jung amd Koh''s test not available for multi-equation models')
    end
    
    if ~est.isPanel
        error('Baltagi, Song, Jung amd Koh''s test requieres a panel data estimation')
    end

    if est.isRobust
        error('Baltagi, Song, Jung amd Koh''s test not valid for robust estimations');
    end   
    
    % Create otuput structure
    test = testout();
    
    % Get used data
    n = est.n;
    T = est.T;
    W = est.W;
    
    % Check size of W matrix
    if size(W,1) ~= n
        error('BSJK test can only be computed if the model is estimated with a W matrix of size n.')
    end
    
    % Get OLS residuals
    pols = ols(est.y,est.X);
    res = pols.res;
    
    % Build A
    A = (res' * kron(eye(n),ones(T)) * res)/(res'*res) - 1;
    % A = (res' * kron(ones(T), eye(n)) * res)/(res'*res) - 1;
    
    % Create bi-diagonal matrix
    G = zeros(T);
    for i=1:T-1
        G(i,i+1) = 1;
        G(i+1,i) = 1;
    end
    
    % Build F
    % The change in the order of the kron is becuase of the article order
    % the data
    F = 1/2 * (((res' * kron(eye(n),G)) * res) / (res'*res));

    % Build H
    % The change in the order of the kron is becuase of the article order
    % the data
    H = 1/2 * ((res' * (kron(W'+W,eye(T))) * res) / (res'*res));
    
    % Build b
    b = trace(W^2 + W'*W);
    
    % Compute LM
    LM = (n*T.^2)/(2*(T-1)*(T-2)) * (A^2 - 4*A*F + 2*T*F^2) + (n^2*T)/b * H^2;
    df = 3;
    p = 1 - chi2dist(LM,df);
    
    % Store results
    test.test = 'BSJKSA';
    test.value = LM;
    test.p = p;
    test.df = df;
    test.isAsymptotic = 0;
    test.isRobust = 0;
    
end


