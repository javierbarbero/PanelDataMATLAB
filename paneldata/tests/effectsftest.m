function [ test ] = effectsftest( est )
%EFFECTSFTEST F test of individual effects
%   Computes the F test of individual effects
%
%   testo = EFFECTSFTEST( est ) Computes the F test of individual effects
%   for the specified estimation output. Returns a test output structure, 
%   testout.
%
%   Example
%     
%      test = effectsftest(est);
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
        error('F tests of individual effects not available for multi-equation models')
    end
    %{
    if est.isInstrumental
        error('F tests of individual effects is not valid for instrumental estimations')
    end
    %}
    
    if est.isRobust
        error('F tests of individual effects is not available after a robust estimatin')
    end
    
    if ~strcmpi(est.options.method, 'fe')
        error('F tests of individual effects can only be applied after a FE estimation')
    end

    % Create otuput structure
    test = testout();

    % Get used data
    y = est.y;
    X = est.X;
    res = est.res;
    resdf = est.resdf;
    n = est.n;
    N = est.N;
    k = est.k;
    
    % Add constant term to X
    X = [X ones(N,1)];
    
    % Compute residuals
    res_pool = (y - X*((X'*X)\X'*y));
    RRSS = res_pool'*res_pool;
    URSS = res'*res;

    % Value of the tests
    value = ((RRSS - URSS)/(n-1)) / (URSS/(N-n-k)); % N-n instead of n*(T-1) (for Unbalanced compatibility)
    df = [n-1, resdf];
    p = 1 - fdist(value, df(1), df(2));          
    
    % Save test results
    test.test = 'EFFECTSF';
    test.value = value;
    test.p = p;
    test.df = df;
    test.isAsymptotic = 0;
    test.isRobust = 0;


end

