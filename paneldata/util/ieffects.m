function [ ieff, se, t, p ] = ieffects( est, options )
%IEFFECTS Computes individual effects after a panel data estimation
%   Computes individual effects after a panel data estimation
%
%   [ieff, se, t, p] = IEFFECTS( est ) computes individual effects after a
%   panel data estimation. Returns the individual effets, ieff, standard
%   errors, se, t-statistics, t, and p-values, p.
%   [ieff, se, t, p] = IEFFECTS( est, 'overall' ) computes the 'overall' 
%   individual effect.
%
%   Example
%     
%      [ieff, se, t, p] = ieffects(est);
%      [ieff, se, t, p] = ieffects(est, 'overall');
%
%   See also ESTOUT, ESTCIDISP
%
%   Copyright 2013-2015 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 17, June, 2015
%

    if ~strcmpi(est.options.method,'fe')
        error('Fixed effects can only be extracted after FE estimation');
    end
    
    if nargin == 2
        if strcmp(options,'overall')
            overall = 1;
        elseif isempty(options)
            overall = 0;
        else
            error('Invalid option');
        end
    else
        overall = 0;
    end
    
    id = est.id;
    n = est.n;
    Tid = est.Tid;
    y = est.y;
    X = est.X;
    coef = est.coef;
    resvar = est.resvar;
    varcoef = est.varcoef;
    resdf = est.resdf;
    
    % Compute group means of y and X
    ymeans = groupmeans(id, y);
    Xmeans = groupmeans(id, X);
    
    % Compute individual effects
    ieff = ymeans - Xmeans * coef;
    
    if overall 
        % Compute overall fixed effects
        ieff = mean(ieff);
        
        % Standard error
        se = mean(Xmeans) * varcoef * mean(Xmeans)';
        se = sqrt(se);
        
    else
        % Individual fixed effects      

        % Standard errors of the individual effect (Greene, p.402)
        se = nan(n,1);
        for i=1:n
            se(i) = resvar/Tid(i) + Xmeans(i,:) * varcoef * Xmeans(i,:)';
        end
        se = sqrt(se);

    end
    
    % T statistic
    t = ieff ./ se;

    % p-value
    p = (1 - tdist(abs(t), resdf))' * 2;


end

