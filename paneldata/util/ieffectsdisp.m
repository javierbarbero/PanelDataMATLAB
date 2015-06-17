function [ ] = ieffectsdisp( est, options )
%IEFFECTSDISP Displau individual effects after a panel data estimation
%   Displau individual effects after a panel data estimation
%
%   [ieff, se, t, p] = IEFFECTSDISP( est ) display individual effects after a
%   panel data estimation. Returns the individual effets, ieff, standard
%   errors, se, t-statistics, t, and p-values, p.
%   [ieff, se, t, p] = IEFFECTSDISP( est, 'overall' ) display the 'overall' 
%   individual effect.
%
%   Example
%     
%      [ieff, se, t, p] = ieffectsdisp(est);
%      [ieff, se, t, p] = ieffectsdisp(est, 'overall');
%
%   See also ESTOUT, ESTCIDISP
%
%   Copyright 2013-2015 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 17, June, 2015
%

    if nargin == 2
        if strcmp(options,'overall')
            overall = 1;
        else
            error('Invalid option');
        end
    else
        overall = 0;
        options = [];
    end

    % Extract fixxed effects
    [ ieff, se, t, p ] = ieffects( est, options );
    
    % Get needed data from estimation
    if overall
        n = 1;
        sid_uniq(1) = cellstr('OVERALL');
    else
        n = est.n;
        sid_uniq = cellstr(num2str(unique(est.id)));
    end
    
    % Header
    fprintf('<strong>Individual Effects</strong> \n');
    
    % Coeffiient significance
    test = waldsigtest(est);
    if est.isAsymptotic
        statname = 'z-stat';
    else
        statname = 't-stat';
    end
    
    % Std err name
    switch(est.vartype)
        case 'homo'
            stderrname = ' Std. Error';
        case 'robust'
            fprintf('Standard errors robust to heteroskedasticity\n');
            stderrname = 'Rob.Std.Err';  
        case 'cluster'
            fprintf('Standard errors robust to heteroskedasticity adjusted for %d clusters\n',length(unique(est.options.clusterid)));
            stderrname = 'Rob.Std.Err';  
        otherwise
            error('Unknown ''vartype''');
    end
    
    % Display table    
    fprintf('---------------------------------------------------------------\n');
    fprintf('     id  |       ieffect   %s      %s    p-value\n',stderrname,statname);
    fprintf('---------------------------------------------------------------\n');
    for i=1:1:n
        fprintf(' %7s | %13.6f  %12.6f  %10.4f  %7.3f ',char(sid_uniq(i)),ieff(i),se(i),t(i), p(i));
        
        % p-value stars
        if p(i) < 0.01
            fprintf('***');
        elseif p(i) < 0.05
            fprintf('** ');
        elseif p(i) < 0.10
            fprintf('*  ');
        end

            % New line
            fprintf('\n');
    end
    fprintf('---------------------------------------------------------------\n');

end

