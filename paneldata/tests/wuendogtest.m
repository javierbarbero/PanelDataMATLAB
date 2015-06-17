function [ test ] = wuendogtest( est )
%WUENDOGTEST Wu's variable addition test of endogeneity
%   Computes the Wu's variables addition test of endogeneity.
%
%   testo = WUENDOGTEST( est ) Computes the Wu's variable addition test for
%   the specified estimation output. Returns a test output structure,
%   testout.
%
%   Example
%     
%      test = wuendogtest(est);
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
        error('Wu''s test not available for multi-equation models')
    end
    
    if est.isPanel
        error('Wu''s test not implemented for panel data models')
    end
    
    % Create otuput structure
    test = testout();
    
    % Get used data from the regression
    nendog = est.nendog;
    nexog = est.nexog;
    N = est.N;
    isRobust = est.isRobust;
    
    % Get endogenous and exogenous variables
    Xendog = est.X(:,est.endog);
    Xexog = est.X(:,est.exog);
    
    % Build Z: Xexog + Z + const
    Z = [est.Z Xexog ones(N,1)];
    
    % Get fitted values from the endogenous variables
    Xendoghat = (Z*((Z'*Z)\Z'*Xendog));
    
    % Test regression
    if ~isRobust
        treg = ols(est.y,[Xendog Xexog Xendoghat]);
    else
        treg = ols(est.y,[Xendog Xexog Xendoghat],'vartype','robust');
    end
    
    % Compute Wald joint sig test for the coefficients of res1
    R = [zeros(nendog,nendog+nexog) eye(nendog) zeros(nendog,1)];
    r = zeros(nendog,1);
    
    wtest = waldsigtest(treg,R,r);
    
    % Save test results
    test.test = 'WUENDOG';
    test.value = wtest.value;
    test.p = wtest.p;
    test.df = wtest.df;
    test.isAsymptotic = 0;
    test.isRobust = isRobust;

end

