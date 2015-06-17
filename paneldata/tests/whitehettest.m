function [ test ] = whitehettest( est )
%WHITEHETTEST White's heteroskedasticity test.
%   Computes the White's heteroskedasticity test.
%
%   testo = WHITEHETTEST( est ) Computes the White's heteroskedasticity test
%   for the specified estimation output. Returns a test output structure,
%   testout.
%
%   Example
%     
%      test = whitehettest(est);
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
        error('White''s test not available for multi-equation models')
    end

    if est.isPanel
        error('White''s test not available for panel data estimation')
    end
    
    %{
    if est.isInstrumental
        error('White''s test is not valid for instrumental estimations')
    end
    %}
    
    % Create otuput structure
    test = testout();
    
    if est.hasConstant
        k = est.k-1;
    else
        k = est.k;
    end
    
    if est.isInstrumental
        % Add the instruments to the number of variables
        % and substract the endogenous variables
        k = k + size(est.Z,2) - est.nendog;
    end
        
    % Compute squared residualds
    res2 = est.res.^2;
    
    % Build X matrix for the test
    if ~est.isInstrumental
        X = est.X;
    else
        Xexog = est.X(:,est.exog);
        X = [Xexog, est.Z];
    end
    
    Xtest = X;
    for i=1:k
        for j=i:k
            Xtest = [Xtest X(:,i).*X(:,j)];
        end
    end
    
    % Remove duplicated variables (squares of dummies, etc.)
    kall = size(Xtest,2);
    Xtest = unique(Xtest','rows');
    Xtest = Xtest';
    
    % Number of duplicates removed.
    duprm = kall - size(X,2);
    
    % Regress squared residuals on X
    whitereg = ols(res2, Xtest);
    
    % Compute Score, df and p-value
    Score = whitereg.N * whitereg.r2;
    df = size(Xtest,2);
    p = 1 - chi2dist(Score, df);
    
    % Save test results
    test.test = 'WHITEHET';
    test.value = Score;
    test.p = p;
    test.df = df;
    test.isAsymptotic = 1;
    test.isRobust = 0;
    test.isInstrumental = 1;
    
    % whitehet specific
    test.duprm = duprm;

end

