function [ est ] = ols( y, X, varargin )
%OLS Ordinary Least Squares estimation
%   Computes Ordinary Least Squares (OLS) estimation for the specified data
%
%   est = OLS(y, X) computes an OLS estimation of y into X. Returns an
%   estimation ouput structure, estout.
%   est = OLS(y, X, Name,Value) OLS estimation with additional properties 
%   using one or more Name,Value pair arguments.
%
%   Additional properties:
%   - 'constant': 0 to ommit the constan term (intercept) in the
%   estiamtion. Default 1.
%   - 'vartype': compute heteroskedasticity 'robust' or 'cluster' variance 
%   matrix estiamtion. Default 'homo' for homoskedasticity.
%   - 'clusterid': specify the cluster variable. Only if 'vartype' is set
%   to 'cluster'.
%   - 'dfcorrection': 0 to supress degrees-of-freedom correction of the
%   variance. Default 1.
%
%   Example
%     
%      est = ols(y, X);
%
%   See also ESTOUT, IV2SLS, S2SLS, PANEL
%
%   Copyright 2013-2015 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 16, June, 2016
%

    % Create output structure
    est = estout();
    
    % Check input
    if nargin < 1
        error('Dependent variable not specified');
    end
    if nargin < 2
        error('Independent variables not specified');
    end
    if size(y,2) ~= 1
        error('Y must be a column vector of the dependant variable');
    end
    if size(y,1) ~= size(X,1)
        error('Number of rows in Y must be equal to number of rows in X');
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
                 @(x) any(validatestring(x,{'homo','robust','cluster'})));
    addPar(p,'dfcorrection',1,@(x) isnumeric(x));
    addPar(p,'clusterid',[],@(x) length(x) == length(y))
    p.parse(varargin{:})
    options = p.Results;

    % Extract table names and convert data to array
    [y, ynames] = extracttable(y);
    [X, xnames] = extracttable(X);  
            
    % Error if NaN's in input data
    if any(isnan(y)) || any(any(isnan(X))) || any(any(isnan(options.clusterid)))
        error('NaN values not allowed in input data. Remove all rows with NaN''s before using this function.');
    end
        
    % Check if constant
    hasConstant = options.constant;
    
    % Check if robust
    if strcmp(options.vartype,'robust') || strcmp(options.vartype,'cluster')
        robust = 1;
    else
        robust = 0;
    end
    
    % Store original data
    est.y = y;
    est.X = X;
    
    % Get number of observations
    n_clus = size(y,1);
    N = n_clus;
    
    % Add constant term
    if hasConstant
        X = [X ones(N,1)];
    end
    
    % Number of variables
    k = size(X,2);
    
    % Compute estimates
    coef = (X'*X)\X'*y;
    
    % Fitted values
    yhat = X*coef;
    
    % Residuals
    res = y - yhat;
    
    % Residuasl degrees of freedom
    resdf = N - k;
    
    % Inverse of X'X
    invXX = ((X'*X)\eye(k));
    
    % Residual variance
    resvar = (res'*res);
    resvar = resvar ./ resdf;
        
    % Covariance matrix
    switch options.vartype
        case 'homo'
            varcoef = resvar * invXX;
        case 'robust'
            % White-robust standard errors
            varcoef = invXX * (X' * diag(diag(res*res')) * X) * invXX;            
            if options.dfcorrection
                varcoef = N/(N-k) * varcoef;
            end
        case 'cluster'
            % Get unique ids
            id_clus = options.clusterid;
            uid_clus = unique(id_clus);
            n_clus = length(uid_clus);
            
            % Compute total sum
            total_sum = 0;
            for i=1:1:n_clus;
                total_sum = total_sum + X(id_clus==uid_clus(i),:)'*res(id_clus==uid_clus(i))*res(id_clus==uid_clus(i))'*X(id_clus==uid_clus(i),:);
            end            
            varcoef = invXX * total_sum * invXX; 
            
            if options.dfcorrection
                varcoef = n_clus/(n_clus-1) * (N-1)/(N-k) * varcoef;
            end
            
            % Change resdf to the correct value with cluster
            resdf = n_clus - 1;
    end
    
    % Standard errors
    stderr = sqrt(diag(varcoef));
       
    
    % Goodness of fit
    % M0 = eye(N) - 1/N * ones(N);
    RSS = res' * res;
    TSS = y' * y;
    ESS = TSS - RSS;
    % r2 = 1 - (res' * M0 * res) ./ (y' * M0 * y);
    r2 = 1 - RSS ./ sum((y - mean(y)).^2);
    adjr2 = 1 - (N - 1) ./ (N - k) .* (1 - r2);
    
    % Save estimation
    est.method = 'OLS';
    est.options = options;
    est.n = n_clus;
    est.T = 1;
    est.N = N;
    est.k = k;
    est.isMultiEq = 0;
    est.isLinear = 1;
    est.hasConstant = hasConstant;
    est.isRobust = robust;
    est.coef = coef;
    est.varcoef = varcoef;
    est.stderr = stderr;
    est.yhat = yhat;
    est.res = res;
    est.resvar = resvar;
    est.resdf = resdf;
    est.RSS = RSS;
    est.ESS = ESS;
    est.TSS = TSS;
    est.r2 = r2;
    est.adjr2 = adjr2;
    est.isAsymptotic = 0;
    est.vartype = options.vartype;
    
    % Set default var names
    est.ynames = ynames;
    est.xnames = xnames;
    est = defaultVarNames(est);
    


end

