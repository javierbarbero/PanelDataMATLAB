function [  ] = estdisp( est )
%ESTDISP Display estiamtion results
%   ESTDISP( est ) Display estiamtion results stored in an 'estout'
%   structure.
%
%   Example
%     
%      estdisp(est);
%
%   See also ESTOUT, TESTDISP
%
%   Copyright 2013-2015 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 26, May, 2016
%

    % Fill missing variables names. In case the user has change some names
    % and not specify all of them
    est = defaultVarNames(est);

    fprintf('_________________________________________________________\n');
    
    % Estimation method
    fprintf('<strong>');

    switch(est.method)
        case 'OLS'
            fprintf('Ordinary Least Squares (OLS)');
        case 'IV2SLS'
            fprintf('Two Stage Least Squares (IV2SLS)');
        % Panel
        case 'PO'
            fprintf('Panel: Pooling estimation (PO)')
        case 'FE'
            fprintf('Panel: Fixed effects (within) (FE)');
        case 'BE'
            fprintf('Panel: Between estimation (BE)');
        case 'RE'
            fprintf('Panel: Random effects (RE)');
        % IV Panel
        case 'PO2SLS'
            fprintf('IV Panel: Pooling estimation (PO2SLS)');
        case 'FE2SLS'
            fprintf('IV Panel: Fixed effects two stage least squares (FE2SLS)')
        case 'BE2SLS'
            fprintf('IV Panel: Between estimation two stage least squares (BE2SLS)');
        case 'RE2SLS'
            fprintf('IV Panel: Random effects two stage least squares (RE2SLS)');
        case 'EC2SLS'
            fprintf('Panel: Baltagi''s error components two stage least squares (EC2SLS)');
        % Spatial
        case 'S2SLS'
            fprintf('Spatial Two Stage Least Squares (S2SLS)');
        case 'FES2SLS'
            fprintf('Spatial Panel: Fixed effects spatial two stage least squares (FES2SLS)')
        case 'BES2SLS'
            fprintf('Spatial Panel: Between estimation spatial two stage least squares (BES2SLS)');
        case 'RES2SLS'
            fprintf('Spatial Panel: Random effects spatial two stage least squares (RES2SLS)');
        case 'ECS2SLS'
            fprintf('Spatial Panel: Error components spatial two stage least squares (ECS2SLS)');
        otherwise
            error('Unknown method');
    end
    fprintf('</strong>\n');
    
    % Display equatio info
    pos = 1;
    for i=1:est.g

        disp(' ');
        
        % Display equation if Multi-Eq
        if est.isMultiEq ~= 0
            fprintf('%s equation: \n', char(est.ynames(i)));
        end

        % Observations
        fprintf('N = %d  ',est.N);
        
        % Panel group and time
        if est.isPanel
            fprintf('n = %d  ', est.n);
            
            if est.isBalanced
                fprintf('T = %d ', est.T);
                fprintf('(Balanced panel)')
            else
                fprintf('T = %d,...,%.2f,...,%d ',min(est.Tid),est.Tmean,max(est.Tid));
                fprintf('(Unbalaced panel)')
            end
        end
        fprintf('\n');

        % Goodness of fit
        if ~isnan(est.r2(i))
            fprintf('R-squared = %6.5f',est.r2(i));
        end
        if ~isnan(est.adjr2(i))
            fprintf('    Adj R-squared = %6.5f \n',est.adjr2(i));
        elseif ~isnan(est.r2(i));
            % To prevent a blank line if nor r2 nor adjr2 is available
            fprintf('\n');
        end
        
        % Joint significance test        
        if ~est.isMultiEq
            if est.hasConstant
                R = [eye(est.k-1) zeros(est.k-1,1)];
                r = zeros(est.k-1,1);
            else                
                R = eye(est.k);
                r = zeros(est.k,1);
            end
        else
            % Zeros of previos equations + eye of the equation + zero
            % (constant) + zeros(for the next equations)
        	R = [zeros(est.keq(i)-1,pos-1) eye(est.keq(i)-1) zeros(est.keq(i)-1,1) zeros(est.keq(i)-1, est.k-est.keq(i)-pos+1) ];
            r = zeros(est.keq(i)-1,1);
            
            if est.hasConstant == 0
                error('Not implemented')
            end
        end
   
        
        % Wald test of joint significance
        test = waldsigtest(est,R,r);
        
        if est.isAsymptotic
            % Is Chi-Square
            fprintf('Wald Chi2(%d) = %f ', test.df, test.value);
            fprintf('p-value = %5.4f \n', test.p);
        else
            % Is F
            fprintf('Wald F(%d, %d) = %f ', test.df(1), test.df(2), test.value);
            fprintf('p-value = %5.4f \n', test.p);
        end

        % RSS, ESS and TSS
        if ~est.isMultiEq
            % Only in single-equation
            if ~isnan(est.RSS) || ~isnan(est.ESS) || ~isnan(est.TSS)
                if ~isnan(est.RSS)
                    fprintf('RSS = %f ',est.RSS);
                end
                if ~isnan(est.ESS)
                    fprintf('ESS = %f ',est.ESS);
                end
                if  ~isnan(est.TSS)
                    fprintf('TSS = %f ',est.ESS);
                end
                fprintf('\n');
                %fprintf('RSS = %f ESS = %f TSS %f \n',est.RSS, est.ESS, est.TSS);
            end
        end
        
        % Increase positoin
        pos = pos + est.keq(i);
    end
    
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

    
    disp(' ');
    % Regression table
    pos = 1;
    for i=1:est.g
        fprintf('----------------------------------------------------------------------\n');
        fprintf('%15.15s |   Coefficient   %s      %s    p-value\n',char(est.ynames(i)),stderrname,statname);
        fprintf('----------------------------------------------------------------------\n');
        for j=1:1:est.k   
            fprintf('%15.15s | %13.6f  %12.6f  %10.4f  %7.3f ',char(est.xnames(pos)),est.coef(pos),est.stderr(pos),test.value(pos), test.p(pos));

            % p-value stars
            if test.p(pos) < 0.01
                fprintf('***');
            elseif test.p(pos) < 0.05
                fprintf('** ');
            elseif test.p(pos) < 0.10
                fprintf('*  ');
            end

            % New line
            fprintf('\n');
            
            % Increase position
            pos = pos + 1;
        end
        
    end
    fprintf('----------------------------------------------------------------------\n');
    
    % Spatial eror model
    if est.isSpatial && est.slagerror
        srho = est.srho;
        %stderrsrho = sqrt(est.sOmega(end,end));
        stderrsrho = est.sigma_srho;
        tsrho = srho / stderrsrho;
        psrho = (1 - normaldist(abs(tsrho)))' * 2;
        fprintf('            rho | %13.6f  %12.6f  %10.4f  %7.3f ',srho, stderrsrho, tsrho, psrho);
        
        % p-value stars
        if psrho < 0.01
            fprintf('***');
        elseif psrho < 0.05
            fprintf('** ');
        elseif psrho < 0.10
            fprintf('*  ');
        end

        % New line
        fprintf('\n');
            
        fprintf('----------------------------------------------------------------------\n');
    end
    
    % PANEL RE variance information
    if est.isPanel && (strcmpi(est.options.method,'re') || strcmp(est.options.method,'ec'))
        if ~est.isSpatial
            fprintf('sigma_mu = %f      rho_mu = %f \n',sqrt(est.sigma2_mu), est.rho_mu);
        end
        fprintf(' sigma_v = %f     sigma_1 = %f \n',sqrt(est.sigma2_v ), sqrt(est.sigma2_1));    
        
        if est.isBalanced
            fprintf('   theta = %f\n', est.theta);
        else
            fprintf('   theta = %f,...,%f,...,%f\n', min(est.theta), median(est.theta), max(est.theta));
        end
        
        fprintf('----------------------------------------------------------------------\n');
        
    end
    
    
    
    % IV information
    if est.isInstrumental 
        
        fprintf('Endogenous: ');
        for i=1:est.nendog
            fprintf('%s ',char(est.xnames(est.endog(i))));
        end
        
        fprintf('\n');
        
        % Display information about instruments only if estimation is not
        % spatial
        
        if ~est.isSpatial
        
            fprintf('Instruments (exogenous): ');
            for i=1:est.nexog
                fprintf('%s ',char(est.xnames(est.exog(i))));
            end

            fprintf('\n');

            fprintf('Instruments (new): ');
            for i=1:est.lnew
                fprintf('%s ',char(est.znames(i)));
            end
            fprintf('\n');
        end
        
        fprintf('----------------------------------------------------------------------\n');        
    end
    
    
    
    
    
    disp(' ');

end

