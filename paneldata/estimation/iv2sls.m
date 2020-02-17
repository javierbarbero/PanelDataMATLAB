function [ est ] = iv2sls( y, X, Z, varargin )
%IV2SLS Instrumental Variables Two Stage Least Squares estimation
%   Computes Instrumental Variables Two Stage Least Squares estimation (IV2SLS)
%   estimation for the specified data and instruments.
%
%   est = IV2SLS(y, X) computes an IV2SLS estimation of y into X, using the
%   additional instruments in Z. Returns an estimation ouput structure, estout.
%   est = IV2SLS(y, X, Z, Name,Value) IV2SLS estimation with additional 
%   properties using one or more Name,Value pair arguments.
%
%   Additional properties:
%   - 'constant': 0 to ommit the constan term (intercept) in the
%   estiamtion. Default 1.
%   - 'vartype': compute heteroskedasticity 'robust' variance matrix 
%   estiamtion. Default 'homo' for homoskedasticity.
%   - 'endog': a vector of indices of the endogenous variables in X. This
%   option is mandatory.
%   - 'useinstruments': use the specified instruments instead of the ones
%   specified in Z plus the exogenous X.
%
%   Example
%     
%      est = iv2sls(y, X, Z);
%
%   See also ESTOUT, OLS, S2SLS, IVPANEL
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
    if nargin < 3
        error('Independent variables not specified');
    end
    if size(y,2) ~= 1
        error('Y must be a column vector of the dependant variable');
    end
    if size(y,1) ~= size(X)
        error('Number of rows in Y must be equal to number of rows in X');
    end
    %{
    if size(y,1) ~= size(Z,1)
        error('Number of rows in Y must be equal to number of rows in Y and X');
    end
    %}    
        
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
                 @(x) any(validatestring(x,{'homo','robust'})));
    addPar(p,'endog',[],@(x) isnumeric(x) && size(x,2) > 0 && size(x,2) < size(X,2));
    addPar(p,'useinstruments',[],@(x) isnumeric(x));
    p.parse(varargin{:})
    options = p.Results;
    
    % Extract table names and convert data to array
    [y, ynames] = extracttable(y);
    [X, xnames] = extracttable(X);  
    [Z, znames] = extracttable(Z);
    
    % Error if NaN's in input data
    if any(isnan(y)) || any(any(isnan(X))) || any(any(isnan(Z))) || any(any(isnan(options.useinstruments)))
        error('NaN values not allowed in input data. Remove all rows with NaN''s before using this function.');
    end
    
    % Check if constant
    hasConstant = options.constant;
    
    % Check if robust
    if strcmp(options.vartype,'robust')
        robust = 1;
    else
        robust = 0;
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
    
    % Store original data
    est.y = y;
    est.X = X;
    est.Z = Z;
    
    % Get number of observations
    n = size(y,1);
    N = n;
    
    % Add constant term
    if hasConstant
        X = [X ones(N,1)];
    end
    
    % Build matix of instruments
    if hasConstant
        Z = [X(:,exog) Z ones(N,1)];
    else
        Z = [X(:,exog) Z];
    end
    
    
    % Number of variables
    k = size(X,2);
    nendog = size(X(:,endog),2);
    nexog = size(X(:,exog),2);
    
    % Replace instruments if specified in the useinstruments option
    if ~isempty(options.useinstruments)
        Z = options.useinstruments;
        lsef = NaN;
        lnew = size(Z,2);
        l = lnew;
    end
    
    % First Stage
    Xhat = Z*((Z'*Z)\Z'*X);

    % Compute estimates: 2SLS
    coef = (Xhat'*X)\Xhat'*y;
    
    % Fitted values
    yhat = X*coef;
    
    % Residuals
    res = y - yhat;
    
    % Residuasl degrees of freedom
    resdf = N;
    
    % Inverse of X'X
    invXhatX = ((Xhat'*X)\eye(k));
    
    % Residual variance
    resvar = (res'*res) ./ resdf;
    
    % Covariance matrix
    if ~robust
        varcoef = resvar * invXhatX;
    else
        % White-robust standard errors
        varcoef = invXhatX * (Xhat' * diag(diag(res*res')) * Xhat) * invXhatX;
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
    est.method = 'IV2SLS';
    est.n = n;
    est.T = 1;
    est.N = N;
    est.k = k;
    est.l = l;
    est.lself = lself;
    est.lnew = lnew;
    est.isMultiEq = 0;
    est.isLinear = 1;
    est.isInstrumental = 1;
    est.endog = endog;
    est.exog = exog;
    est.nendog = nendog;
    est.nexog = nexog;
    est.hasConstant = 1;
    est.isRobust = robust;
    est.coef = coef;
    est.varcoef = varcoef;
    est.stderr = stderr;
    est.Xhat = Xhat;
    est.yhat = yhat;
    est.res = res;
    est.resvar = resvar;
    est.resdf = resdf;
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

