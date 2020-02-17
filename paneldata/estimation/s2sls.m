function [ est ] = s2sls( y, X, W, varargin )
%S2SLS Spatial two stage least squares
%   Computes a Spatial Two Stage Least Squares(S2SLS) estimation.
%
%   est = S2SLS(y, X) computes a S2SLS estimation of y into X, using the
%   contiguity matrix W. Returns an estimation ouput structure, estout.
%   est = S2SLS(y, X, Z, Name,Value) S2SLS estimation with additional 
%   properties using one or more Name,Value pair arguments. Use this to
%   specify the spatial lags to include in the model
%
%   Additional properties:
%   - 'slagy': 1 to include a spatial lag of the dependent variables.
%   Default 1.
%   - 'slagerror': 1 to include a spatil lag of the error structure.
%   Default 0.
%   - 'slagX': a vector of indices to specify the spatial lags of the
%   variables in X to be added. Default [].
%   - 'endog': a vector of indices of endogenous variables in X.
%   - 'inst': instruments for the endogenous variables in X.
%
%   Example
%     
%      est = s2sls(y, X, W, 'slagy', 1);
%      est = s2sls(y, X, W, 'slagy', 1, 'slagerror', 1);
%
%   See also ESTOUT, OLS, IV2SLS, SPANEL
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
    addPar(p,'vartype','homo',...
                 @(x) any(validatestring(x,{'homo'})));
    addPar(p,'slagy',1,@(x) isnumeric(x));
    addPar(p,'slagX',[],@(x) isnumeric(X));
    addPar(p,'slagerror',0,@(x) isnumeric(X));
    addPar(p,'endog',[],@(x) isnumeric(x) && size(x,2) > 0 && size(x,2) < size(X,2));
    addPar(p,'inst',[],@(x) isnumeric(x));
    p.parse(varargin{:})
    options = p.Results;

    % Extract table names and convert data to array
    [y, ynames] = extracttable(y);
    [X, xnames] = extracttable(X);    
    options.inst = extracttable(options.inst);
            
    % Error if NaN's in input data
    if any(isnan(y)) || any(any(isnan(X))) || any(any(isnan(W))) || any(any(isnan(options.inst)))
        error('NaN values not allowed in input data. Remove all rows with NaN''s before using this function.');
    end
        
    % Get endogenous and exogenous ariables 
    endog = unique(options.endog);
    exog = 1:size(X,2);
    exog(ismember(exog,endog)) = [];   
    
    % Get and check number of instruments
    lself = size(X(:,exog),2);
    Z = options.inst;
    lnew = size(Z,2);
    l = lself + lnew;
    
    if l < size(X,2)
        error('Number of instruments must be at least equal to the number of variables');
    end
    
    % Build X matrix for the Spatial model
    Xspa = X;
    slagy = options.slagy;
    if slagy
        % Add spatial lag of y;
        Xspa = [Xspa W*y];
        endog = [endog, size(Xspa,2)];
    end
    
    slagX = unique(options.slagX);
    noslagX = 1:size(X,2);
    noslagX(ismember(noslagX,slagX)) = [];
    if ~all(isnan(slagX)) && ~isempty(slagX)
        % Add spatial lags of indicated X
        Xspa = [Xspa W*X(:,slagX)];
    else
        slagX = [];
    end
    
    % Build matrix of  instruments   
    exnoslagX = noslagX;
    exnoslagX(~ismember(exnoslagX,exog)) = [];

    N = size(y,1);
    
    
    % -------------------------------------
    % Step 1.1: Initial estimation
    % And unique if there is no spatial correlation in the error term
    % -------------------------------------
    % Compute regression regression
    if endog
        if slagy
            H = [W*X(:,exnoslagX), W^2*X(:,exog), W^3*X(:,slagX), Z];
        else
            H = Z;
        end
    
        estTemp = iv2sls(y, Xspa, H,'endog',endog);
    else
        estTemp = ols(y, Xspa);
    end
        
    
    % If spatial lag in the error term
    slagerror = options.slagerror;
    if slagerror
        
        % Build matrix of instruments (incluidng exogenous X)
        if slagy
            H = [X(:,exog) W*X(:,exnoslagX), W^2*X(:,exog), W^3*X(:,slagX), Z, ones(N,1)];
        else
            H = [X(:,exog), Z, ones(N,1)];
        end
        
        % Get spatial matrix for the errors
        M = W;
                
        % Get residuals
        u = estTemp.res;
        
        % Build u bar
        ubar = W * u;
        
        % -------------------------------------
        % Step 1.2: Get initial estimate of rho
        % -------------------------------------        
        [A1, A2] = computeA1A2(M);
        
        Gamma = N^-1 * [u'*(A1+A1')*ubar, - ubar'*A1*ubar;
                        u'*(A2+A2')*ubar, - ubar'*A2*ubar];
        gamma = N^-1 * [u'*A1*u; u'*A2*u];
        
        perrfnc = @(x) (Gamma*[x; x^2]-gamma)'*(Gamma*[x; x^2]-gamma);
        
        options = optimset('TolFun',1e-10,'TolX',1e-10);
        rho = fminsearch(perrfnc,0.2,options);
        
        % -------------------------------------
        % Step 2.1: GS2SLS estimation
        % -------------------------------------
        % Cochrane-Ocutt transformation
        % Add the constant here to X since it need to be transformed
        ytr = y - rho*M*y;
        Xspatr = [Xspa ones(N,1)] - rho*M*[Xspa ones(N,1)];
        
        % Estimate the transformed model
        if endog
            estTemp = iv2sls(ytr, Xspatr, [],'endog',endog,'constant',0,'useinstruments',H);
        else
            estTemp = ols(ytr, Xspatr,'constant',0);
        end
        
        % -------------------------------------
        % Step 2.2: Consistent and efficient GMM estimator of rho
        % -------------------------------------
        
        % Get residuals using untransformed y and X
        u = y - [Xspa ones(N,1)]*estTemp.coef;     

        % Build u bar
        ubar = W * u;

        A1 = (1 + ((N^-1) * trace(M'*M))^2)^-1 * (M'*M - N^-1 * trace(M'*M) * eye(N));
        A2 = M;
        
        Gamma = N^-1 * [u'*(A1+A1')*ubar, - ubar'*A1*ubar;
                        u'*(A2+A2')*ubar, - ubar'*A2*ubar];
        gamma = N^-1 * [u'*A1*u; u'*A2*u];       
               
        Z = Xspatr;        
        Psi = computePsi(rho, u, H, Z, W);
        
        perrfnc = @(x) (Gamma*[x; x^2]-gamma)' * ((Psi)\eye(size(Psi))) * (Gamma*[x; x^2]-gamma);
        
        rho = fminsearch(perrfnc,0.2,options);  
        
        % -------------------------------------
        % Correct Variance/Covariance computation
        % -------------------------------------
        
        % Transform Z again with the new rho
        Z = [Xspa ones(N,1)] - rho*M*[Xspa ones(N,1)];
        
        
        [Psi, sigma2, P, QHH, a1, a2, mu3] = computePsi(rho, u, H, Z, W);

        
        % Psi elements
        Psi_dd = sigma2 * QHH;
        Psi_dr = sigma2 * N^-1 * (H' * [a1, a2]) + mu3 * N^-1 * (H' * [diag(A1), diag(A2)]);
         
        % J marix
        J = Gamma * [1; 2*rho];
        
        % Omega elements
        invPsi = Psi \ eye(size(Psi));
        Omega_dd = P' * Psi_dd * P;
        Omega_rr = inv(J' * invPsi * J);
        Omega_rr = (Omega_rr\eye(size(Omega_rr))); % Invert it
        Omega_dr = P' * Psi_dr * invPsi * J * Omega_rr;
          
        % Omega matrix
        Omega = N^-1 * [Omega_dd, Omega_dr; Omega_dr', Omega_rr];

        % Substitute the correct variance in the estimation output
        varcoef = Omega(1:end-1,1:end-1);
        stderr = sqrt(diag(Omega(1:end-1,1:end-1)));
        % est.sOmega = Omega;
        sigma_srho = sqrt(Omega(end,end));

        
    else
        rho = NaN;
        sigma_srho = NaN;
        varcoef = estTemp.varcoef;
        stderr = sqrt(diag(varcoef));
    end
    
    
    % Number of regressors (+1 for the constant term)
    hasConstant = 1;
    k = size(Xspa,2) + 1;
    Xspa = [Xspa, ones(N,1)];
    
    % Store original data (shorted)
    est.y = y;
    est.X = X;
    est.W = W;
    
    % Coefficients, fitted values and residuals
    n = N;
    T = 1;
    coef = estTemp.coef;
    yhat = Xspa * coef;
    res = y - yhat;
    resdf = estTemp.resdf;
    RSS = res'*res;
    ESS = NaN;
    TSS = NaN;
    options.vartype = 'homo';
    isInstrumental = ~isempty(endog);
    nendog = length(endog);
    nexog = length(exog);
    
    temp = corrcoef(yhat,y);
    r2 = temp(1,2).^2;
    adjr2 = NaN;
    
    % Store original data (shorted)
    est.y = y;
    est.X = X;
    est.W = W;

    % Save estimation
    est.method = 'S2SLS';
    est.options = options;
    est.n = n;
    est.T = T;
    est.N = N;
    est.k = k;
    est.isPanel = 0;
    est.isMultiEq = 0;
    est.isLinear = 1;
    est.hasConstant = hasConstant;
    est.isRobust = 0;
    est.l = l;
    est.lself = lself;
    est.lnew = lnew;
    est.isInstrumental = isInstrumental;
    est.endog = endog;
    est.exog = exog;
    est.nendog = nendog;
    est.nexog = nexog;
    est.coef = coef;
    est.varcoef = varcoef;
    est.stderr = stderr;
    %est.Xhat = Xhat;
    est.yhat = yhat;
    %est.yhattr = yhattr;
    est.res = res;
    %est.resvar = resvar;
    est.resdf = resdf;
    
    est.RSS = RSS;
    est.ESS = ESS;
    est.TSS = TSS;
    est.r2 = r2;
    est.adjr2 = adjr2;    
    est.isAsymptotic = 1;
    est.vartype = options.vartype;
    
    
    % Spatial
    est.isSpatial = 1;
    est.slagy = slagy;
    est.slagX = slagX;
    est.slagerror = slagerror;
    est.srho = rho;
    est.sigma_srho = sigma_srho;

    % Set default var names
    est.ynames = ynames;
    est.xnames = xnames;
    est = defaultVarNames(est);
    

end

function [A1, A2] = computeA1A2(M)

    N = size(M,1);

    A1 = (1 + ((N^-1) * trace(M'*M))^2)^-1 * (M'*M - N^-1 * trace(M'*M) * eye(N));
    A2 = M;
        
end

function [Psi, sigma2, P, QHH, a1, a2, mu3] = computePsi(rho, u, H, Z, M)

    % Get A1 and A2
    [A1, A2] = computeA1A2(M);

% Get N
    N = size(u,1);

% Get epsilon by transforming resitudals
        epsilon = u - rho * (M * u);

% Q matrices
        QHH = (N^-1*(H'*H));
        iQHH = inv(QHH);
        QHZ = (N^-1*(H'*Z));   
        
        % P and T matrices
        P = iQHH * QHZ * inv(QHZ'*iQHH*QHZ);
        T = H*P;
        
        % alpha and a
        alpha1 = -N^-1 * (Z'*(A1+A1')*epsilon);
        a1 = T*alpha1;
        alpha2 = -N^-1 * (Z'*(A2+A2')*epsilon);
        a2 = T*alpha2;
        
        sigma2 = 1/N * (epsilon'*epsilon);
        mu3 = N^-1 * sum(epsilon.^3);
        mu4 = N^-1 * sum(epsilon.^4);
        
        % Psi
        Psi12 = sigma2^2 * (2*N)^-1 * trace( (A1 + A1')*(A2+A2') ) + ...
            sigma2 * N^-1 * (a1'*a2) + ...
            N^-1 * (mu4 - 3*sigma2^2) * diag(A1)' * diag(A2) + ...
            N^-1 * mu3 * (a1' * diag(A2) + a2' * diag(A1));
        
        Psi21 = sigma2^2 * (2*N)^-1 * trace( (A2 + A2')*(A1+A1') ) + ...
            sigma2 * N^-1 * (a2'*a1) + ...
            N^-1 * (mu4 - 3*sigma2^2) * diag(A2)' * diag(A1) + ...
            N^-1 * mu3 * (a2' * diag(A1) + a1' * diag(A2));
        
        Psi11 = sigma2^2 * (2*N)^-1 * trace( (A1 + A1')*(A1+A1') ) + ...
            sigma2 * N^-1 * (a1'*a1) + ...
            N^-1 * (mu4 - 3*sigma2^2) * diag(A1)' * diag(A1) + ...
            N^-1 * mu3 * (a1' * diag(A1) + a1' * diag(A1));
        
        Psi22 = sigma2^2 * (2*N)^-1 * trace( (A2 + A2')*(A2+A2') ) + ...
            sigma2 * N^-1 * (a2'*a2) + ...
            N^-1 * (mu4 - 3*sigma2^2) * diag(A2)' * diag(A2) + ...
            N^-1 * mu3 * (a2' * diag(A2) + a2' * diag(A2));
        
        Psi = [Psi11 Psi12; Psi21 Psi22];
    

end


% OLD CODE:
% KP1998
    %{
    if slagerror
        
        % Get residuals
        u = est.res;
        
        % Build u bar
        ubar = W * u;
        sum(ubar);
        ubarbar = W * ubar;
        sum(ubarbar);
        
        % Build A matrices
        N = est.N;
        
        A11 =  1/N * sum(ubar.^2);
        A12 = -2/N * sum(u .* ubar);
        A13 = -1;
        A21 =  1/N * sum(ubarbar.^2);
        A22 = -2/N * sum(ubar .* ubarbar);
        A23 = -1/N * trace(W'*W);
        A31 =  1/N * sum(ubar .* ubarbar);
        A32 = -1/N * (sum(ubar.^2) + sum(u .* ubarbar));
        A33 =  0;
        
        A1 = [A11 A12 A13; A21 A22 A23; A31 A32 A33];
        
        A2 = [-1/N*sum(u.^2); -1/N*sum(ubar.^2); -1/N*sum(u.*ubar)];
        
        rho = A1\A2;
        
        srho = rho(2)
        
        %  Transform the model
        y = (eye(N) - srho*W)*y;
        Xspa = (eye(N) - srho*W)*Xspa;
        
        % Estimate the transformed model
        if slagy
            est = iv2sls(y, Xspa, H,'endog',endog);
        else
            est = ols(y, Xspa);
        end
        
    end
    %}
    
