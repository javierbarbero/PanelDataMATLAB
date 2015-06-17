function [ test ] = testout(  )
%TESTOUT Generates a default testout structure to store results of a test
%   Generates a default testout structure to store results of a test
%
%   Common fields:
%   - test: name of the test performed.
%   - value: value of the test.
%   - df: degrees of freedom of the test.
%   - p: p value of the test.
%   - isAsymptotic: 1 if the test statistic is asymptotically distributed.
%   - isRobust: 1 if the test have been performed robust to
%   heteroskedasticity.
%
%   Example
%     
%      test = testout();
%
%   See also TESTDISP
%
%   Copyright 2013-2015 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 17, June, 2015
%

    test.test = NaN;            % Test name
    test.value = NaN;           % Value/score of the test
    test.df    = NaN;           % Degrees of freedom
    test.p     = NaN;           % Associated p-value    
    test.isAsymptotic = NaN;    % Asymptotic
    test.isRobust = NaN;        % Test robust to heteroskedasticity

end

