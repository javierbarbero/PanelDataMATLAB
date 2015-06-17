function [ test ] = resettest( est, uppowers )
%RESETTEST Ramsey RESET test
%   Computes the Ramsey RESET test
%
%   testo = RESETTEST( est ) Computes the Ramsey RESET test for the specified
%   estimation output. Returns a test output structure, testout.
%   testo = RESETTEST( est, uppowerd ) Specify the powers to compute the
%   test. Default 4.
%
%   Example
%     
%      test = resettest(est);
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
        error('Ramsey test not available for multi-equation models')
    end
    
    if est.isPanel
        error('Ramsey test not available for panel dataestimations')
    end
    
    if est.isInstrumental
        error('Ramsey test not implemented for IV estimations');
    end
    
    
    % Create otuput structure
    test = testout();

    % Check input
    if nargin < 2
        uppowers = 4;
    else
        if uppowers < 2
            error('At least yhat.^2 requiered');
        end
    end
    
    % Compute matrix of yPowers
    ypow = [];
    for i=2:uppowers      
        ypow(1:est.N,end+1) = est.yhat.^i;
    end
    
    % Regresseion    
    ramseyols = ols(est.y, [est.X, ypow]);
    
    % F test
    R = [zeros(uppowers-1,est.k-1), eye(uppowers-1), zeros(uppowers-1,1)];
    r = zeros(uppowers-1,1);
    wsigtest = waldsigtest(ramseyols,R,r);
    
    value = wsigtest.value;
    p = wsigtest.p;
    df = wsigtest.df;
    
    % Reset test results
    test.test = 'RESET';
    test.value = value;
    test.p = p;
    test.df = df;
    test.isAsymptotic = 0;
    test.isRobust = est.isRobust;
    
    % resettest specific
    test.powers = uppowers;

end

