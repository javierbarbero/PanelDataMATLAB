function [ test ] = sarganoitest( est )
%SARGANOITEST Sargan overidentification test.
%   Computes the Sargan overidentification test.
%
%   testo = SARGANOITEST( est ) Computes the Sargan overidentification test
%   for the specified estimation output. Returns a test output structure, 
%   testout.
%
%   Example
%     
%      test = sarganoitest(est);
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
        error('Sargan''s test not available for multi-equation models')
    end
    
    %{
    if est.isPanel
        error('Sargan''s test not available for panel data estimations')
    end
    %}
    if est.isPanel && strcmpi(est.options.method,'be')
        error('Sargan''s test not available for between estimation');
    end
    
    if ~est.isInstrumental
        error('Sargan''s test is for instrumental estimations')
    end   
        
    if est.lnew <= est.nendog
        error('Sargan''s test of overidentification requieres the model to be over-identified');
    end
    
    % Create otuput structure
    test = testout();
    
    % Get used data
    Xexog = est.X(:,est.exog);
    Z = est.Z;
    res = est.res;
    
    if est.isRobust
        % Wooldridge Score
        
        % Perform individual regressions for each new instrument
        % Get residuals from regressen each new instrument on Xexog, Xhat and
        % 1
        Xsar = est.Xhat;
        resr = Z - (Xsar*((Xsar'*Xsar)\Xsar'*Z));
        % Compute ur
        ur = repmat(est.res,1,est.lnew) .* resr;
        
        % Score regression
        sarreg = ols(ones(est.N,1), ur,'constant',0);
        
        % Compute Score, df and p-value
        Score = sarreg.N - sarreg.RSS;
        df = est.lnew - est.nendog;
        p = 1 - chi2dist(Score,df);       
                
    else
        % Regress residuals on exogneous and new instruments
        if ~est.isPanel
            sarreg = ols(res, [Xexog, Z]);
        else
            sarreg = panel(est.id, est.time, res, [Xexog, Z],lower(est.options.method));
        end

        % Compute Score, df and p-value
        if ~est.isPanel || strcmpi(est.method,'re')
            Score = sarreg.N * sarreg.r2;
        else
            Score = (sarreg.N - sarreg.n) * sarreg.r2;
        end
        df = est.lnew - est.nendog;
        p = 1 - chi2dist(Score,df);
    end

    % Save test results
    test.test = 'SARGANOI';
    test.value = Score;
    test.p = p;
    test.df = df;
    test.isAsymptotic = 1;
    test.isRobust = est.isRobust;
    test.isPanel = est.isPanel;

end

