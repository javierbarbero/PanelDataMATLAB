function [ est ] = ivpanel( id, time, y, X, Z, method, varargin )
%IVPANEL Instrumental Variables Panel data estimation
%   Computes instrumental variables panel data estimation for fixed effects,
%   between, and random effects.
%
%   est = IVPANEL(id, time, y, X, Z, method) computes instrumental variables
%   panel data estimation of y into X, using the additional instruments in Z,
%   with the specified method, where id and time are the individual and time
%   identifiers. The function returns an estimation ouput structure, estout.
%   est = IVPANEL(id, time, y, X, Z, method, Name,Value) computes panel data 
%   estimation  with additional properties using one or more Name,Value pair
%   arguments.
%
%   method is one of the following:
%   - 'po': pooling estimation
%   - 'fe': fixed effects (whitin) estimation
%   - 'be': between estimation
%   - 're': random effects GLS estimation
%   - 'ec': Baltagi's error components estimation.
%
%   Additional properties:
%   - 'constant': 0 to ommit the constan term (intercept) in the
%   estiamtion. Default 1.
%   - 'endog': a vector of indices of the endogenous variables in X. This
%   option is mandatory.
%   - 'useinstruments': use the specified instruments instead of the ones
%   specified in Z plus the exogenous X.
%
%   Example
%     
%      est = ivpanel(id, time, y, X, Z, 'fe');
%      est = ivpanel(id, time, y, X, Z, 'be');
%      est = ivpanel(id, time, y, X, Z, 're');
%      est = ivpanel(id, time, y, X, Z, 'ec');
%
%   See also ESTOUT, IV2SLS, PANEL, SPANEL
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
    if nargin < 6
        error('Panel method not specified')
    end      
    
    % Parse Additional options
    p = inputParser;
    if verLessThan('matlab', '8.2')
        addPar = @(v1,v2,v3,v4) addParamValue(v1,v2,v3,v4);
    else
        addPar = @(v1,v2,v3,v4) addParameter(v1,v2,v3,v4);
    end
    addPar(p,'constant',1,...
                 @(x) isnumeric(x));
    addPar(p,'vartype','homo',...
                 @(x) any(validatestring(x,{'homo'})));
    addPar(p,'endog',[],@(x) isnumeric(x) && size(x,2) > 0 && size(x,2) <= size(X,2));
    addPar(p,'useinstruments',[],@(x) isnumeric(x));
    p.parse(varargin{:})
    options = p.Results;

    % Extract table names and convert data to array
    id = extracttable(id);
    time = extracttable(time);
    [y, ynames] = extracttable(y);
    [X, xnames] = extracttable(X);  
    [Z, znames] = extracttable(Z);
            
    % Error if NaN's in input data
    if any(isnan(id)) ||any(isnan(time)) || any(isnan(y)) || any(any(isnan(X))) || any(any(isnan(Z))) || any(any(isnan(options.useinstruments)))
        error('NaN values not allowed in input data. Remove all rows with NaN''s before using this function.');
    end
        
    % Check method
    if ~ismember(method,{'po','fe','be','re','ec'})
        error('Incorrect method')
    end
    
    % Get endogenous and exogenous ariables 
    endog = unique(options.endog);
    if length(endog) < 1
        error('No endogenous variables specified. Use the ''endog'' option.')
    end
    exog = 1:size(X,2);
    exog(ismember(exog,endog)) = [];
    
    % Get and check number of instruments
    lself = size(X(:,exog),2);
    lnew = size(Z,2);
    l = lself + lnew;
    
    if l < size(X,2) && isempty(options.useinstruments)
        error('Number of instruments must be at least equal to the number of variables');
    end
                 
    % Get number of observations
    N = size(y,1);
    
    % Get balanced and T variables
    [ isBalanced, idx, n, T, Tid, Tmean, Thmean] = isbalancedpanel( id, time );
    
    % Sort panel
    id = id(idx);
    time = time(idx);
    y = y(idx,:);
    X = X(idx,:);
    if isempty(options.useinstruments)
        Z = Z(idx,:);
    else
        options.useinstruments = options.useinstruments(idx,:);
    end
    
    % Store original data (shorted)
    est.id = id;
    est.uid = unique(id);
    est.time = time;
    est.idx = idx;
    est.y = y;
    est.X = X;
    est.Z = Z;
    
        
    % Number of varaibles and hasConstant
    k = size(X,2);
    nendog = size(X(:,endog),2);
    nexog = size(X(:,exog),2);
    if ~strcmp(method,'fe')
        k = k + 1;  % For the constant term in non fe estiamtions
        nexog = nexog;
        hasConstant = 1;
    else
        hasConstant = 0;
    end
    
    % Transform variables
    switch method
        case 'po'
            % Add constant term
            X = [X ones(N,1)];
            
            % Build matix of instruments            
            Z = [X(:,exog) Z ones(N,1)];
                        
            % Replace instruments if specified in the useinstruments option
            if ~isempty(options.useinstruments)
                Z = options.useinstruments;
                lself = NaN;
                lnew = size(Z,2);
                l = lnew;
            end
            
            % No transformation por pool
            ytr = y;
            Xtr = X;
            Ztr = Z;
            
            resdf = N;
        case 'fe'      
            % Build matix of instruments
            Z = [X(:,exog) Z];

            % Replace instruments if specified in the useinstruments option
            if ~isempty(options.useinstruments)
                Z = options.useinstruments;
                lself = NaN;
                lnew = size(Z,2);
                l = lnew;
            end
            
            % Fixed effects transformation
            ytr = y - groupmeans(id,y,'replicate',1);
            Xtr = X - groupmeans(id,X,'replicate',1);
            Ztr = Z - groupmeans(id,Z,'replicate',1);
            
            resdf = (N-n)-k; % (N-n) instead of (n*T) for compatibility with unbalanced panels

        case 'be'
            % Add constant term
            X = [X ones(N,1)];
            
            % Build matix of instruments
            Z = [X(:,exog) Z ones(N,1)];
                       
            % Replace instruments if specified in the useinstruments option
            if ~isempty(options.useinstruments)
                Z = options.useinstruments;
                lself = NaN;
                lnew = size(Z,2);
                l = lnew;
            end
            
            % Between transformation.
            ytr = groupmeans(id, y);
            Xtr = groupmeans(id, X);
            Ztr = groupmeans(id, Z);
            
            resdf = n-k;
        case {'re','ec'}
            % Replace instruments if specified in the useinstruments option
            if ~isempty(options.useinstruments)
                Z = options.useinstruments;
                lself = NaN;
                lnew = size(Z,2);
                l = lnew;
            end
            
            % Get time invariant variables to remove from the first 'fe'
            % and 'be' regressions
            tinvariantX = istinvariant(id,X);
            tinvariantZ = istinvariant(id,Z);
            
            % Stop if all endogenous, endogenous or Z variables are time invariant
            if all(tinvariantX(:,endog) == 1)
                error('Not all endogenous variables can be time-invariant.')
            end
            if all(tinvariantX(:,exog) == 1)
                error('Not all exogenous variables can be time-invariant.')
            end
            if all(tinvariantZ == 1)
                error('Not all new instruments can be time-invariant.')
            end
            
            % Get sigma2_v from 'fe' residuals
            estfe = ivpanel(id,time,y,X(:,~tinvariantX),Z(:,~tinvariantZ),'fe','endog',endog,'useinstruments',options.useinstruments);
            sigma2_v = estfe.resvar;
            
            % Get sigma2_1 from 'be' residuals
            estbe = ivpanel(id,time,y,X(:,~tinvariantX),Z(:,~tinvariantZ),'be','endog',endog,'useinstruments',options.useinstruments);
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
                Tid_g = repelem(Tid,Tid,1);
                theta = 1 - sqrt(sigma2_v./(Tid_g.*sigma2_mu + sigma2_v));
            end
                                    
            % Build matix of instruments
            Z = [X(:,exog) Z];
            
            % Replace instruments if specified in the useinstruments option
            if ~isempty(options.useinstruments)
                Z = options.useinstruments;
                lself = NaN;
                lnew = size(Z,2);
                l = lnew;
            end
            
            % If Baltagi's error components EC2SLS
            if strcmp(method,'ec')
                Zmean = groupmeans(id,Z,'replicate',1);
                Z = [Z - Zmean, Zmean];
            end
            
            % Add constant term
            X = [X ones(N,1)];
            if isempty(options.useinstruments)
                Z = [Z ones(N,1)];       
            end
            
            % Random effects transformation
            if isBalanced
                ytr = y - theta * groupmeans(id,y,'replicate',1);            
                Xtr = X - theta * groupmeans(id,X,'replicate',1);
                Ztr = Z - theta * groupmeans(id,Z,'replicate',1);
            else
                ytr = y - theta .* groupmeans(id,y,'replicate',1);            
                Xtr = X - repmat(theta,1,k) .* groupmeans(id,X,'replicate',1);
                Ztr = Z - repmat(theta,1,k) .* groupmeans(id,Z,'replicate',1);
            end
            
            resdf = N-k;
                       
    end  
    
    % First Stage
    Xhat = Ztr*((Ztr'*Ztr)\Ztr'*Xtr);

    % Compute estimates: 2SLS
    coef = (Xhat'*Xtr)\Xhat'*ytr;
    
    % Fitted values
    yhattr = Xtr*coef;
    yhat = X*coef;
    
    % Residuals
    res = ytr - yhattr;
    
    % Inverse of Xtr'Xtr
    invXtrXtr = ((Xhat'*Xtr)\eye(k));
    
    % Residual variance
    resvar = (res'*res) ./ resdf;
    
    % Covariance matrix
    varcoef = resvar * invXtrXtr;
    
    % Standard errors
    stderr = sqrt(diag(varcoef));
    
    % Goodness of fit
    if strcmp(method,'be')
    %    M0 = eye(n) - 1/n * ones(n);
        adjr2_correction = (n - 1) ./ (resdf);
    else
    %    M0 = eye(N) - 1/N * ones(N);
        adjr2_correction = (N - 1) ./ (resdf);
    end
    RSS = res' * res;
    TSS = y' * y;
    ESS = TSS - RSS;
    % r2 = 1 - (res' * M0 * res) ./ (ytr' * M0 * ytr);
    r2 = 1 - RSS ./ sum((ytr - mean(ytr)).^2);
    adjr2 = 1 - adjr2_correction .* (1 - r2);
    
    % Compute correc r2 for 're'
    if strcmp(method,'re') || strcmp(method,'ec')
        temp = corrcoef(yhat,y);
        r2 = temp(1,2).^2;
        adjr2 = 1 - adjr2_correction .* (1 - r2);
    end

    
    % Save estimation
    est.method = strcat(upper(method),'2SLS');
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
    est.isRobust = 0;
    est.l = l;
    est.lself = lself;
    est.lnew = lnew;
    est.isInstrumental = 1;
    est.endog = endog;
    est.exog = exog;
    est.nendog = nendog;
    est.nexog = nexog;
    est.coef = coef;
    est.varcoef = varcoef;
    est.stderr = stderr;
    est.Xhat = Xhat;
    est.yhat = yhat;
    est.yhattr = yhattr;
    est.res = res;
    est.resvar = resvar;
    est.resdf = resdf;
    
    est.Tid = Tid;
    est.Tmean = Tmean;
    est.Thmean = Thmean;
    if strcmp(method,'re') || strcmp(method,'ec')
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
    est.isAsymptotic = 1;
    est.vartype = options.vartype;
        
    % Set default var names
    est.ynames = ynames;
    est.xnames = xnames;
    est.znames = znames;
    est = defaultVarNames(est);
    
    

end

