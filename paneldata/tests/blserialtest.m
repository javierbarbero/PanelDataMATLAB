function [ test ] = blserialtest( est )
%BLSERIALTEST Baltagi and Li's serial correlation test.
%   Computes the Baltagi and Li's test for serial correlation and random
%   effects.
%
%   testo = BLSERIALTEST( est ) Computes the Baltagi and Li's test for the 
%   specified estimation output. Returns a test output structure, testout.
%
%   Example
%     
%      test = blserialtest(est);
%
%   See also TESTOUT
%
%   Copyright 2013-2015 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 17, June, 2015
%

    if est.isMultiEq
        error('Baltagi and Li''s test not available for multi-equation models')
    end

    if est.isRobust
        error('Baltagi and Li''s test not valid for robust estimations');
    end
    
    if ~strcmpi(est.options.method,'re')
        error('Baltagi and Li''s test must be performed after a random effects estimation');
    end
    
    if est.isInstrumental
        error('Baltagi and Li''s test not valid for instrumental estimation')
    end
    
    
    % Create otuput structure
    test = testout();
    
    % Get used data
    n = est.n;
    T = est.T;
    id = est.id;
    uid = est.uid;
    
    % Get OLS residuals
    pols = ols(est.y,est.X);
    res = pols.res;
    
    % Build A
    A = (res' * kron(eye(n),ones(T)) * res)/(res'*res) - 1;
    
    % Build B
    rescrossLag = nan(length(uid),1);
    rescross = nan(length(uid),1);
    for i=1:1:length(uid)
        resi = res(id == uid(i));
        rescross(i) = resi(2:end)'*resi(2:end);
        rescrossLag(i) = resi(1:end-1)'*resi(2:end);
    end
    rescross = sum(rescross);
    rescrossLag = sum(rescrossLag);
    
    B = rescrossLag/rescross;

    % Compute LM
    LM = (n*T.^2)/(2*(T-1)*(T-2)) * (A^2 - 4*A*B + 2*T*B^2);
    df = 2;
    p = 1 - chi2dist(LM,df);
    
    % Store results
    test.test = 'BLSERIAL';
    test.value = LM;
    test.p = p;
    test.df = df;
    test.isAsymptotic = 0;
    test.isRobust = 0;
    
end

