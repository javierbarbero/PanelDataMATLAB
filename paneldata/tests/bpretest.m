function [ test ] = bpretest( est )
%BPRETEST Breusch and Pagan LM test for random effects
%   Computes the Baltagi and Li (1990) version of the Breusch and Pagan (1980)
%   test for random effects.
%
%   testo = BPRETEST( est ) Computes the Breusch and Pagan LM test for 
%   random effects for the specified estimation output. Returns a test output
%   structure, testout.
%
%   Example
%     
%      test = bpretest(est);
%
%   See also TESTOUT
%
%   Copyright 2013-2015 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 17, June, 2015
%

    if ~est.isPanel
        error('Breusch-Pagan''s test for random effects is for panel data estiamtions')
    end
    
    if ~strcmpi(est.options.method,'re');
        error('Breusch-Pagan''s test for random effects is for panel data RE estiamtions')
    end

    if est.isInstrumental
        error('Breusch-Pagan''s test for random effects is not valid for instrumental estimations')
    end
    
    % Create otuput structure
    test = testout();
    
    % Run pooled regression
    poest = ols(est.y, est.X);
    res = poest.res;
    n = est.n;
    
    % If balanced
    if est.isBalanced
        T = est.T;
        
        % Compute numerator
        res2num = 0;
        for i=1:n
            s = (i-1)*T+1;
            e = i*T;
            res2num = res2num + sum(res(s:e)).^2;
        end
        % Compute denominator
        res2den = sum(res.^2);
        
        % Compute LM
        if est.sigma2_mu >= 0
            lambdaLM = (n*T)/(2*(T-1)) * (res2num/res2den - 1)^2;
        else
            lambdaLM = 0;
        end
    else
        % Unbalanced Panel
        
        Tid = est.Tid;
        Tmean = est.Tmean;
        
        % Compute numerator
        res2num = 0;
        lastt = 1;
        for i=1:n
            s = lastt;
            e = lastt + Tid(i) - 1;
            res2num = res2num + sum(res(s:e)).^2;
            lastt = lastt + Tid(i);
        end
        % Compute denominator
        res2den = sum(res.^2);
        
        % Compute LM
        if est.sigma2_mu >= 0
            A = 1 - (res2num/res2den);
            lambdaLM = (n*Tmean).^2/2 * ((A.^2)/(sum(Tid.^2) - n*Tmean));
        else
            lambdaLM = 0;
        end
    end
    
    df = 1;
    p = 1 - chi2dist(lambdaLM, df);
    
    % Save test results
    test.test = 'BPRE';
    test.value = lambdaLM;
    test.p = p;
    test.df = df;
    test.isAsymptotic = 1;
    test.isRobust = 0;

end

