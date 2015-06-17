function [ test ] = hausmantest( estA, estB )
%HAUSMANTEST Hausman test
%   Computes the Hausman test
%
%   testo = HAUSMANTEST( estA, estB ) Computes the Hausman test comparing
%   the two specified estimation output. Returns a test output structure, 
%   testout.
%
%   Example
%     
%      test = hausmantest(estA, estB);
%
%   See also TESTOUT
%
%   Copyright 2013-2015 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 17, June, 2015
%

    if estA.isMultiEq || estB.isMultiEq
        error('Hausman''s test not available for multi-equation models')
    end

    if estA.isRobust || estB.isRobust
        error('Hausman''s test not valid for robust estimations');
    end
        
    % Create otuput structure
    test = testout();
    
    % Get coefficients and covariance matrices    
    kA = estA.k - estA.hasConstant;
    kB = estB.k - estB.hasConstant;
    if kA ~= kB
        error('Number of estimated coefficients must be the same in both estimations');
    end
    ktest = kA;
    %{
    if estA.hasConstant & estB.hasConstant;
        ktest = estA.k-1; % Test without the constant
    else
        ktest = estA.k; 
    end
    %}
    
    % Set default var names
    estA = defaultVarNames(estA);
    estB = defaultVarNames(estB);
    
    % Check if variable names match bwetten estimations
    if any(~strcmp(estA.ynames,estB.ynames))
        warning('yvarnames differs bwetten estimations');
    end
    
    if any(~strcmp(estA.xnames(1:ktest),estB.xnames(1:ktest)))
        warning('xvarnames differ bwetten estimations');
    end
    
    % Extract coefficients and variance matrix
    coefA = estA.coef(1:ktest);
    coefB = estB.coef(1:ktest);
    varcoefA = estA.varcoef(1:ktest,1:ktest);
    varcoefB = estB.varcoef(1:ktest,1:ktest);
    
    
    % Compute Hausman's statistic
    H = (coefA - coefB)' * ((varcoefA - varcoefB)\eye(ktest)) * (coefA - coefB);
    
    if H <= 0
        warning('Invalid estimates. H <= 0');
    end
    
    % Distributed as a Chi-Square with ktest degrees of freedom
    df = ktest;
    p = 1 - chi2dist(H, df);
    
    % Save test results
    test.test = 'HAUSMAN';
    test.value = H;
    test.p = p;
    test.df = df;
    test.isAsymptotic = 0;
    test.isRobust = 0;
    
    % hasumantest specific
    test.estA.k = estA.k;
    test.estB.k = estB.k;    
    test.estA.method = estA.method;
    test.estB.method = estB.method;    
    test.estA.coef = estA.coef;
    test.estB.coef = estB.coef;    
    test.estA.varcoef = estA.varcoef;
    test.estB.varcoef = estB.varcoef;    
    test.estA.xnames = estA.xnames;

end

