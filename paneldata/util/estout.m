function [ est ] = estout(  )
%TESTOUT Generates a default estout structure to store estiamtion results
%   Generates a default estout structure to store estiamtion results
%
%   Estimation:
%   - method: estimation method.
%   - options: estimation options.
%
%   Model properties:
%   - hasConstant: 1 if the model has a constant term (intercep)
%   - isLinear: 1 if the model is a linear model
%   - isInstrumental: 1 if instrumental variables estimation
%   - isPanel: 1 if panel data estimation
%   - isRobust: 1 if heteroskedasticity-robust estimation
%   - isAsymptotic: 1 if asymptotic variance matrix
%   - isSpatial: 1 if spatial estimation
%   - isMultiEq: 1 if multi equation model (reserved for future use)
%   
%   Variable names:
%   - yvarnames: names of the dependent variables
%   - Xvarnames: names of the explanatory variables
%   - Zvarnaems: names of the instruments
%
%   Data matrices:
%   - y: dependent variable
%   - X: independent variable
%   - Z: matrix of instruments
%   - W: contiguity matrix of spatial models.
%   - Xhat: first stage estimated X's (for IV)
%
%   Observations, time and variables
%   - n: Number of units (for panel)
%   - T: Number of time periods (average if unbalanced)
%   - N: Total number of observations
%   - k: Number of regressors
%   
%   IV specific
%   - endog: indices of endogenous variables in X.
%   - exog: indices of exogenous variables in X.
%   - nendog: Number of endogenous variables
%   - nexog: Number of exogenous variables
%   - l: Number of total instuments
%   - lself: Number of self-instrumetns
%   - lnew: Number of new instruments.
%
%   Panel specific
%   - Tid: Number of time periods for each id
%   - id: ID's of units
%   - time: Time of units
%   - idx: Indexes of shorted data
%   - isBalanced: 1 if panel is balanced
%   - yhattr: Fitted values of the transformed model
%
%   Spatial specific
%   - slagy: 1 if the model includes a spatial lag in y
%   - slagerror: 1 if the model includes a spatial lag in the error structure.
%   - slagX: indices of the spatial lags of the X's included in the model
%   - srho: estimated coefficient of the spatial lag in the error structure.
%   - sigma_srho: standard error of the previous
%
%   Estiamtion results
%   - coef: Estimated coefficients
%   - yhat: Fitted values
%   - res: Residuals
%   - resvar: Residual variance
%   - resdf: Residuals degrees-of-freedom
%   - varcoef: Variance-covariance matrix
%   - stderr: Coef. Std. Errors
%   - RSS: Residual Sum of Squares
%   - ESS: Explained Sum of Squares
%   - TSS: Total Sum of Squares
%   - r2: R-squared
%   - adjr2: Adjusted R-squared
%
%
%   Example
%     
%      test = testout();
%
%   See also ESTDISP
%
%   Copyright 2013-2015 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 17, June, 2015
%

    % Estimation method
    est.method = NaN;
    est.options = NaN;

    % Model properties
    est.hasConstant = 0;       % The mode has constant term
    est.isLinear = 0;          % The model is a linear model
    est.isInstrumental = 0;    % Instrumental variables estimation
    est.isPanel = 0;           % Panel data estimation
    est.isRobust = 0;          % Robust estimation
    est.isAsymptotic = 0;    % Asymptotic covariance matrix
    est.isSpatial = 0;         % Spatial estimation
    est.isMultiEq = 0;         % Multi-equation model
    
    % Variable names
    est.ynames = NaN;       % y variables names
    est.xnames = NaN;       % X variable names
    est.znames = NaN;       % Z variables names
    
    % Data matrices
    est.y = NaN;               % Dependent variables
    est.X = NaN;               % Independent variables
    est.Z = NaN;               % Instrumental variables
    est.W = NaN;               % W matrix (spatial)
    est.Xhat = NaN;            % First-stage estiamted X's (for IV)
    
    % Observations, time and variables
    est.n = NaN;               % Number of units (for panel)
    est.T = NaN;               % Number of time periods (average if unbalanced)
    est.N = NaN;               % Total number of observations
    est.k = NaN;               % N. of regressors (including CONST)
    
    % IV specific
    est.nendog = NaN;           % Number of endogenous variables
    est.nexog = NaN;            % Number  of exogenous variables
    est.l = NaN;               % N. of total instuments
    est.lself = NaN;           % N. of self-instrumetns
    est.lnew = NaN;            % N. of new instruments.

    % Multiequation specifig
    est.g = 1;                 % Number of equations
    est.keq = NaN;             % Number of regressors in each equation 
 
    
    % Panel specific
    est.Tid = NaN;             % Number of time periods for each id
    est.id = NaN;              % ID's of units
    est.time = NaN;            % Time of units
    est.idx = NaN;             % Indexes of shorted data
    est.isBalanced = NaN;      % If panel is balanced
    est.yhattr = NaN;          % Fitted values of the transformed model
    
    % Spatial specific
    est.slagy = NaN;
    est.slagerror = NaN;
    est.slagX = NaN;
    est.srho = NaN;
    est.sigma_srho = NaN;
    
    % Estiamtion results
    est.coef = NaN;            % Estimated coefficients
    est.yhat = NaN;            % Fitted values
    est.res = NaN;             % Residuals
    est.resvar = NaN;          % Residual variance (homoskedasticity)
    est.resdf = NaN;           % Residuals degrees-of-freedom
    est.varcoef = NaN;         % Variance-covariance matrix
    est.stderr = NaN;          % Coef. Std. Errors
    est.RSS = NaN;             % Residual Sum of Squares
    est.ESS = NaN;             % Explained Sum of Squares
    est.TSS = NaN;             % Total Sum of Squares
    est.r2 = NaN;              % R-squared
    est.adjr2 = NaN;           % Adjusted R-squared

end

