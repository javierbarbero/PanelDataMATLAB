function [ test ] = mundlakvatest( est )
%MUNDLAKVATEST Mundlak variable addition test
%   Computes the Mundlak variable addition test
%
%   testo = MUNDLAKVATEST( est ) Computes the Mundlak variable addition test for
%   the specified estimation output. Returns a test output structure, testout.
%
%   Example
%     
%      test = mundlakvatest(est);
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
        error('Mundlak''s variable addition test not available for multi-equation models')
    end
    
    if ~est.isPanel
        error('Mundlak''s variable addition test requieres a panel data estimation')
    end
    
    if est.isSpatial
        error('Mundlak''s variable addition test not implemented for spatial panels')
    end
    
    % This is to ensure that only time-invariant variables are included
    if ~strcmpi(est.options.method,'fe')
        error('Mundlak''s variable addition test must be applied after an FE estimation')
    end
    
    if est.isInstrumental
        error('Mundlak''s variable addition test not for instrumental variables estimation')
    end
    
    % Create otuput structure
    test = testout();
    
    % Get variables of interest
    id = est.id;
    time = est.time;
    k = est.k;
    isRobust = est.isRobust;
        
    % Compute group means variables and replicate for all observations
    Xbar = groupmeans(id,est.X,'replicate',1);
    
    % Test regression ignoring thr Xbar in the 'first step' of the RE
    % regression
    if ~isRobust
        treg = panel(id, time, est.y, [est.X Xbar],'re');
    else
        treg = panel(id, time, est.y, [est.X Xbar],'re','vartype','robust');
    end
    
    % Compute Wald joint sig test for the coefficients of res1
    R = [zeros(k,k) eye(k) zeros(k,1)];
    r = zeros(k,1);
    
    wtest = waldsigtest(treg,R,r);
    
    % Save test results
    test.test = 'MUNDLAKVA';
    test.value = wtest.value;
    test.p = wtest.p;
    test.df = wtest.df;
    test.isAsymptotic = 0;
    test.isRobust = isRobust;

end

