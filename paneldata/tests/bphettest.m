function [ test ] = bphettest( est )
%BPHETTEST Breusch and Pagan heteroskedasticity test
%   Computes the Breusch and Pagan heteroskedasticity test
%
%   testo = BPHETTEST( est ) Computes the Breusch and Pagan heteroskedasticity
%   test for the specified estimation output. Returns a test output
%   structure, testout.
%
%   Example
%     
%      test = bphettest(est);
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
        error('Breusch-Pagan''s test not available for multi-equation models')
    end
    
    if est.isPanel
        error('Breusch-Pagan''s test not available for panel data estiamtions')
    end

    %{
    if est.isInstrumental
        error('Breusch-Pagan''s test is not valid for instrumental estimations')
    end
    %}
    
    % Create otuput structure
    test = testout();
    
    % Compute squared residualds
    res2 = est.res.^2;
    
    % Regress squared residuals on X
    if ~est.isInstrumental
        bpreg = ols(res2, est.X);
    else
        bpreg = ols(res2, [est.X(:,est.exog), est.Z]);
    end
    
    % Compute Score, df and p-value
    Score = bpreg.N * bpreg.r2;
    df = bpreg.k - 1;
    p = 1 - chi2dist(Score, df);
    
    % Save test results
    test.test = 'BPHET';
    test.value = Score;
    test.p = p;
    test.df = df;
    test.isAsymptotic = 1;
    test.isRobust = 0;
    test.isInstrumental = est.isInstrumental;

end

