function [ est ] = panel( id, time, y, X, method, varargin )
%PANEL Panel data estimation
%   Computes panel data estimation for fixed effects, between, and random
%   effects.
%
%   est = PANEL(id, time, y, X, method) computes a panel data estimation of
%   y into X with the specified method, where id and time are the individual
%   and time identifiers. The function returns an estimation ouput structure,
%   estout.
%   est = PANEL(id, time, y, X, method, Name,Value) computes panel data 
%   estimation with additional properties using one or more Name,Value pair
%   arguments.
%
%   method is one of the following:
%   - 'po': pooling estimation
%   - 'fe': fixed effects (whitin) estimation
%   - 'be': between estimation
%   - 're': random effects GLS estimation
%
%   Additional properties:
%   - 'vartype': compute heteroskedasticity 'robust' or 'cluster' variance 
%   matrix estiamtion. Default 'homo' for homoskedasticity.
%   - 'clusterid': specify the cluster variable. Only if 'vartype' is set
%   to 'cluster'. Default equal to id.
%   - 'dfcorrection': 0 to supress degrees-of-freedom correction of the
%   variance. Default 1.
%
%   Example
%     
%      est = panel(id, time, y, X, 'fe');
%      est = panel(id, time, y, X, 'be');
%      est = panel(id, time, y, X, 're');
%
%   See also ESTOUT, OLS, IVPANEL, SPANEL
%
%   Copyright 2013-2017 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 9, August, 2017
%

    % Create output structure
    est = estout();
    
    % Check input
    if nargin < 3
        error('Dependent variable not specified');
    end
    if nargin < 4
        error('Independent variables not specified');
    end
    if size(y,2) ~= 1
        error('Y must be a column vector of the dependant variable');
    end
    if size(y,1) ~= size(X,1)
        error('Number of rows in Y must be equal to number of rows in X');
    end
    if size(id,1) ~= size(y,1)
        error('Number of rows in id must be equal to number of rows in Y');
    end
    if size(time,1) ~= size(y,1)
        error('Number of rows in time must be equal to number of rows in Y');
    end
    if nargin < 5
        error('Panel method not specified')
    end   
        
    % Parse Additional options
    p = inputParser;
    if verLessThan('matlab', '8.2')
        addPar = @(v1,v2,v3,v4) addParamValue(v1,v2,v3,v4);
    else
        addPar = @(v1,v2,v3,v4) addParameter(v1,v2,v3,v4);
    end
    addPar(p,'vartype','homo',...
                 @(x) any(validatestring(x,{'homo','robust','cluster'})));
    addPar(p,'dfcorrection',1,@(x) isnumeric(x));
    addPar(p,'clusterid',id,@(x) length(x) == length(y))
    p.parse(varargin{:})
    options = p.Results;
    
    % Extract table names and convert data to array
    id = extracttable(id);
    time = extracttable(time);
    [y, ynames] = extracttable(y);
    [X, xnames] = extracttable(X);  
    options.clusterid = extracttable(options.clusterid);

    % Error if NaN's in input data
    if any(isnan(id)) ||any(isnan(time)) || any(isnan(y)) || any(any(isnan(X))) || any(any(isnan(options.clusterid)))
        error('NaN values not allowed in input data. Remove all rows with NaN''s before using this function.');
    end
    
    % Check method
    if ~ismember(method,{'po','fe','be','re'})
        error('Incorrect method')
    end
             
    % Automatically change 'robust' to 'cluster'
    if strcmp(options.vartype,'robust') 
        options.vartype = 'cluster';        
    end
    
    % Check if robust
    if strcmp(options.vartype,'cluster')
        robust = 1;
        
        if strcmp(method,'be')
            error('Robust standard errors not available for between estimation');
        end
    else
        robust = 0;
    end
    
    % Get number of observations
    N = size(y,1);         
    
    % Get balanced and T variables
    [ isBalanced, idx, n, T, Tid, Tmean, Thmean] = isbalancedpanel( id, time );
    
    % Sort variables
    id = id(idx);
    time = time(idx);
    y = y(idx,:);
    X = X(idx,:);
    options.clusterid = options.clusterid(idx); % Sort clusterid
    
    % Store original data (sorted)
    est.id = id;
    est.uid = unique(id);
    est.time = time;
    est.idx = idx;
    est.y = y;
    est.X = X;
               
    % Number of varaibles and hasConstant
    k = size(X,2);
    if ~strcmp(method,'fe')
        k = k + 1;  % For the constant term in non fe estiamtions
        hasConstant = 1;
    else
        hasConstant = 0;
    end
    
    % Transform variables
    switch method
        case 'po'
            % Add constant term
            X = [X ones(N,1)];
            
            % No transformation por pool
            ytr = y;
            Xtr = X;
            
            resdf = N-k;
        case 'fe'      
            % Fixed effects transformation
            ytr = y - groupmeans(id,y,'replicate',1);
            Xtr = X - groupmeans(id,X,'replicate',1);
            
            if ~robust
                resdf = (N-n)-k; % (N-n) instead of (n*T) for compatibility with unbalanced panels
            else
                resdf = n-1;
            end

        case 'be'
            % Add constant term
            X = [X ones(N,1)];
            
            % Between transformation.
            ytr = groupmeans(id, y);
            Xtr = groupmeans(id, X);
            
            resdf = n-k;
        case 're'
            % Get time invariant variables to remove from the first 'fe'
            % and 'be' regressions
            tinvariant = istinvariant(id,X);
            
            % Get sigma2_v from 'fe' residuals
            estfe = panel(id,time,y,X(:,~tinvariant),'fe');
            sigma2_v = estfe.resvar;
            
            % Get sigma2_1 from 'be' residuals
            estbe = panel(id,time,y,X(:,~tinvariant),'be');
            sigma2_1 = estbe.resvar;
            
            % sigma2_u and rho_var
            if isBalanced
                % Calculate sigma2_mu
                sigma2_mu = (sigma2_1 - sigma2_v/T);
            else                                
                % Calculate sigma2_mu using the harmonic mean of T
                sigma2_mu = (sigma2_1 - sigma2_v/Thmean);
            end
            % correct sigma2_mu if < 0
            if sigma2_mu <= 0
                sigma2_mu = 0;
            end
            rho_mu = sigma2_mu / (sigma2_mu + sigma2_v);
            
            % Theta
            if isBalanced
                theta = 1 - sqrt(sigma2_v./(T*sigma2_mu + sigma2_v));
            else
                % Tid_g = repelem(Tid,Tid,1);
                Tid_g = replicate_elements(Tid,Tid);
                theta = 1 - sqrt(sigma2_v./(Tid_g.*sigma2_mu + sigma2_v));
            end
            
            % Add constant term
            X = [X ones(N,1)];
            
            % Random effects transformation
            if isBalanced
                ytr = y - theta * groupmeans(id,y,'replicate',1);            
                Xtr = X - theta * groupmeans(id,X,'replicate',1);
            else
                ytr = y - theta .* groupmeans(id,y,'replicate',1);            
                Xtr = X - repmat(theta,1,k) .* groupmeans(id,X,'replicate',1);
            end
            
            resdf = N-k;
                       
    end
    
    % Compute estimates
    coef = (Xtr'*Xtr)\Xtr'*ytr;
    
    % Fitted values
    yhattr = Xtr*coef;
    yhat = X*coef;
    
    % Residuals
    res = ytr - yhattr;
        
    % Inverse of Xtr'Xtr
    invXtrXtr = ((Xtr'*Xtr)\eye(k));
    
    % Residual variance
    resvar = (res'*res) ./ resdf;
    
    % Covariance matrix
    switch options.vartype
        case 'homo'
            varcoef = resvar * invXtrXtr;
        case {'robust','cluster'}
            % White-robust standard errors adjusted for cluster
            % Get unique ids
            id_clus = options.clusterid;
            uid_clus = unique(id_clus);
            n_clus = length(uid_clus);
            
            % Compute total sum
            total_sum = 0;
            for i=1:1:n_clus;
                total_sum = total_sum + Xtr(id_clus==uid_clus(i),:)'*res(id_clus==uid_clus(i))*res(id_clus==uid_clus(i))'*Xtr(id_clus==uid_clus(i),:);
            end
            varcoef = invXtrXtr * total_sum * invXtrXtr; 
              
            % Degrees of freedom correction
            if options.dfcorrection
                if strcmp(method,'fe')
                    varcoef = (n_clus/(n_clus-1)) * ((N)/(N-k)) .* varcoef;
                else
                    varcoef = (n_clus/(n_clus-1)) * ((N-1)/(N-k)) .* varcoef;
                end
            end
            
            % Change resdf to the correct value with cluster
            resdf = n_clus - 1;
    end

    % Standard errors
    stderr = sqrt(diag(varcoef));
    
    % Goodness of fit
    if strcmp(method,'be')
    %    M0 = eye(n) - 1/n * ones(n);
        adjr2_correction = (n - 1) ./ (resdf);
    else
    %    M0 = eye(N) - 1/N * ones(N);
        if strcmp(method,'fe') && robust
            adjr2_correction = (N - 1) ./ ((N-n)-k); 
        else
            adjr2_correction = (N - 1) ./ (resdf);
        end
    end
    RSS = res' * res;
    TSS = y' * y;
    ESS = TSS - RSS;
    %r2 = 1 - (res' * M0 * res) ./ (ytr' * M0 * ytr);
    r2 = 1 - RSS ./ sum((ytr - mean(ytr)).^2);
    adjr2 = 1 - adjr2_correction .* (1 - r2);

    % Compute correc r2 for 're'
    if strcmp(method,'re')
        temp = corrcoef(yhat,y);
        r2 = temp(1,2).^2;
        adjr2 = 1 - adjr2_correction .* (1 - r2);
    end
    
    % Save estimation
    est.method = upper(method);
    est.options = options;
    est.options.method = method;
    est.n = n;
    est.T = T;
    est.N = N;
    est.k = k;
    est.isPanel = 1;
    est.isBalanced = isBalanced;
    est.isMultiEq = 0;
    est.isLinear = 1;
    est.hasConstant = hasConstant;
    est.isRobust = robust;
    est.coef = coef;
    est.varcoef = varcoef;
    est.stderr = stderr;
    est.yhat = yhat;
    est.yhattr = yhattr;
    est.res = res;
    est.resvar = resvar;
    est.resdf = resdf;
    
    est.Tid = Tid;
    est.Tmean = Tmean;
    est.Thmean = Thmean;
    if strcmp(method,'re')
        est.theta = theta;
        est.sigma2_mu = sigma2_mu;
        est.sigma2_1 = sigma2_1;
        est.sigma2_v = sigma2_v;
        est.rho_mu = rho_mu;
    end
    
    est.RSS = RSS;
    est.ESS = ESS;
    est.TSS = TSS;
    est.r2 = r2;
    est.adjr2 = adjr2;    
    if strcmp(method,'re')
        est.isAsymptotic = 1;
    else
        est.isAsymptotic = 0;
    end
    est.vartype = options.vartype;
            
    % Set default var names
    est.ynames = ynames;
    est.xnames = xnames;
    est = defaultVarNames(est);
    
end

% OLD CODE: 
% Usefull matrices
    %{
    if isBalanced
        Jbar = ones(T)./T;
        P = kron(eye(n),Jbar);  % Group means and replicate
        Q = eye(N) - P;         % Difference from group means
        Tid_g = (T*ones(N,1));  % Time periods per id, repeated for all observations
    else
        P = zeros(N);
        Tid_g = zeros(N,1);
        lastt=1;
        for i=1:n
            Jbar = ones(Tid(i))./Tid(i);
            P(lastt:lastt+Tid(i)-1, lastt:lastt+Tid(i)-1) = Jbar;
            Tid_g(lastt:lastt+Tid(i)-1,1) = Tid(i);
            lastt = lastt + Tid(i);
        end
        Q = eye(N) - P;
    end  
    %}

        %{
        case 'fe':
            ytr = Q*y;
            Xtr = Q*X;
        %}

        %{
        case 're':
            if isBalanced
                Qre = eye(N) - theta * P;
            else
                Qre = eye(N) - repmat(theta,1,N) .* P;
            end
            ytr = Qre * y;
            Xtr = Qre * X;   

        %}
