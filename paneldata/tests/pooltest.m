function [ test ] = pooltest( est )
%POOLTEST Test for poolability of the data
%   Computes the test for poolability of the data
%
%   testo = POOLTEST( est ) Computes the test for poolability of the data 
%   for the specified estimation output. Returns a test output structure, 
%   testout.
%
%   Example
%     
%      test = pooltest(est);
%
%   See also TESTOUT
%
%   Copyright Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 17, June, 2015
%

    if est.isMultiEq
        error('Pool''s test not available for multi-equation models')
    end

    if est.isRobust
        error('Pool''s test not valid for robust estimations');
    end
    
    if strcmpi(est.options.method,'be')
        error('Pool''s test can not be performed after between estimation');
    end
    
    if est.isInstrumental
        error('Pool''s test not valid for instrumetnal estimation')
    end
    
    % Check if one T == 1 in unbalanced
    if est.isBalanced == 0
        if min(est.Tid) == 1
            error('Pool''s test can not be applied if some group have only 1 observation');
        end
    end
    
    % Create otuput structure
    test = testout();
    
    % Get used data
    n = est.n;
    T = est.T;
    y = est.y;
    X = est.X;
    k = est.k;
    id = est.id;
    uid = est.uid;
    
    % Estimate pool ols
    poolols = ols(y, X);
        
    % Regress each individual
    res2i = nan(n,1);
    for i=1:length(uid)
        tempols = ols(y(id == uid(i)),X(id == uid(i),:));
        res2i(i) = tempols.res'*tempols.res;
    end
    
    % Compute F statistic
    df = [(n-1)*(k+1),  n*(T-(k+1))];
    Fnum = (poolols.res'*poolols.res - sum(res2i)) / df(1);
    Fdenom = sum(res2i) / df(2);
    F = Fnum/Fdenom;
    p = 1 - fdist(F, df(1), df(2));
    
    
    % 
    test.test = 'POOL';
    test.value = F;
    test.p = p;
    test.df = df;
    test.isAsymptotic = 0;
    test.isRobust = 0;
    
end

