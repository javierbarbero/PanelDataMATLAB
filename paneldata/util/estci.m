function [ lower, upper ] = estci( est, sig )
%ESTCI Computes confidence intervals
%   Computes confidence intervals for the specified significance level
%
%   [lower, upper] = ESTCI( est, sig ) computes confidence intervals for
%   the specified significance level, sig. Default 0.05. Returns the lower
%   and the upper.
%
%   Example
%     
%      [l, u] = estci(est, 0.05);
%
%   See also ESTOUT, ESTCIDISP
%
%   Copyright 2013-2015 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 17, June, 2015
%

    % Default sig = 0.05 (95%)
    if nargin < 2
        sig = 0.05;
    end
    
    % Report warning if sig < 0.01
    if sig < 0.01
        warning('Very low significance < 0.01');
    end

    % Compute statistic
    if est.isAsymptotic
        statistic = abs(normalinvdist(sig/2));
    else
        statistic = tinvdist(sig/2, est.resdf);
    end
    
    % Compute Lower and Upper CI
    lower = est.coef - statistic .* est.stderr;
    upper = est.coef + statistic .* est.stderr;

end

