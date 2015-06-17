function [  ] = estcidisp( est, sig )
%ESTCIDISP Display confidence intervals
%   Display confidence intervals for the specified significance level
%
%   [lower, upper] = ESTCI( est, sig ) display confidence intervals for
%   the specified significance level, sig. Default 0.05. 
%
%   Example
%     
%      estcidisp(est, 0.05);
%
%   See also ESTOUT, ESTCI
%
%   Copyright 2013-2015 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 17, June, 2015
%

    % Fill missing variables names. In case the user has change some names
    % and not specify all of them
    est = defaultVarNames(est);
    
    % Default sig = 0.05 (95%)
    if nargin < 2
        sig = 0.05;
    end
    
    % Report warning if sig < 0.01
    if sig < 0.01
        warning('Very low significance < 0.01');
    end

    % Compue CI
    [lower, upper] = estci(est,sig);
    
    % Std err name
    if est.isRobust
        stderrname = 'Rob.Std.Err';        
    else
        stderrname = ' Std. Error';
    end
    
    % CI Table
    fprintf('<strong>Confidence Intervals at sig=%3.2f (%2.0f%%)</strong> \n',sig,(100-sig*100));
    fprintf('---------------------------------------------------------------------------\n');
    fprintf('%15.5s |   Coefficient   %s          Lower          Upper\n',char(est.ynames),stderrname);
    fprintf('---------------------------------------------------------------------------\n');
    for i=1:1:est.k   
        fprintf('%15.15s | %13.6f  %12.6f  %13.6f  %13.6f \n',char(est.xnames(i)),est.coef(i),est.stderr(i),lower(i), upper(i));
    end
    fprintf('---------------------------------------------------------------------------\n');
    
    disp(' ');

end

