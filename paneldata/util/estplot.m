function [ fig ] = estplot( est, type, options )
%ESTPLOT Function in progress
%   Function in progress
%
%   Copyright 2013-2015 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 17, June, 2015
%
    
    warning('Function in progress. Results may be worng and syntax can change.')

    % If no options
    if nargin < 3
        options = NaN;
    end

    % Switch plot type
    switch type
        case 'fit'
            fig = estplot_Fit(est);
        case 'res'
            fig = estplot_Res(est);
        case 'resfit'
            fig = estplot_ResFit(est,options);
        otherwise
            error('Unknown estplot type');
            
    end

end

function [ fig ] = estplot_Fit(est)

    fig = figure;
    hold on
    plot(est.y,'b');
    plot(est.yhat,'r--');
    xlim([0 est.N]);
    xlabel('Observation');
    ylabel('Observed and fitted');
    title('Observed and fitted values');
    legend('Observed','Fitted','Location','Best');
    hold off
    
end

function [ fig ] = estplot_Res(est)

    fig = figure;
    hold on
    plot(est.res,'b.');
    plot([0,est.N],[0,0],'k--');
    xlim([0 est.N]);
    xlabel('Observation');
    ylabel('Residuals');
    title('Residuals');
    hold off

end

function [ fig ] = estplot_ResFit(est, options)

    if isnan(options)
        options = 'none';
    end  
    
    fig = figure;
    hold on
    plot(est.yhat, est.res,'b.');
    xlim([min(est.yhat) max(est.yhat)]);    
    % Add a dotted line at 0
    plot([min(est.yhat) max(est.yhat)], [0 0], 'k--');
    
    % Add fit
    if ~strcmp(options,'none')
        % Sort values
        data = [est.res, est.yhat];
        data = sortrows(data, 2);
        
        switch options
            case 'linear'
                % Estimate linear
                fit = ols(data(:,1), data(:,2));
            case 'quad'
                fit = ols(data(:,1), [data(:,2), data(:,2).^2]);
            otherwise
                error('Invalid fit');
        end
        
        % Plot
        % Note: Only the first column of X to prevent multiple lines appearing
        % in a quadratic fit
        plot(fit.X(:,1), fit.yhat,'r-')
    end    
    
    xlabel('Fitted values');
    ylabel('Residuals');
    title('Residuals vs Fitted');
    hold off

end

