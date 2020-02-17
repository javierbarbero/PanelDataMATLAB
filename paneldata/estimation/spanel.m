function [ est ] = spanel( id, time, y, X, W, method, varargin )
%SPANEL Spatial Panel data estimation
%   Computes a spatial panel data estimation.
%
%   est = SPANEL(id, time, y, X, W, method) computes a spatial panel data 
%   estimation of y into X, using the contiguity matrix W,with the specified
%   method, where id and time are the individual and time identifiers. The 
%   function returns an estimation ouput structure, estout.
%   est = SPANEL(id, time, y, X, Z, method, Name,Value) computes spatial 
%   panel data estimation with additional properties using one or more 
%   Name,Value pair arguments.
%
%   method is one of the following:
%   - 'fe': fixed effects (whitin) estimation
%   - 're': random effects GLS estimation
%   - 'ec': Baltagi's error components estimation.
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
%      est = s2sls(id, time, y, X, W, 'fe', 'slagy', 1);
%      est = s2sls(id, time, y, X, W, 'be', 'slagy', 1, 'slagerror', 1);
%
%   See also ESTOUT, S2SLS, PANEL, IVPANEL 
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
    id = extracttable(id);
    time = extracttable(time);
    [y, ynames] = extracttable(y);
    [X, xnames] = extracttable(X);  
    options.inst = extracttable(options.inst);
            
    % Error if NaN's in input data
    if any(isnan(id)) ||any(isnan(time)) || any(isnan(y)) || any(any(isnan(X))) || any(any(isnan(W))) || any(any(isnan(options.inst)))
        error('NaN values not allowed in input data. Remove all rows with NaN''s before using this function.');
    end    
    
    % Check method
    if ~ismember(method,{'fe','re','ec'})
        error('Incorrect method')
    end
    
    % If pool
    if strcmpi(method,'po')
        error('For pooling spatial estimations use the ''s2sls'' function');
    end
    
    % Get balanced and T variables
    [isBalanced, idx, n, T, Tid, Tmean, Thmean] = isbalancedpanel( id, time );
    
    % Store original W matrix
    Worig = W;
    
    % Check W matrix
    if size(W,1) == n
        % Builg Big W matrix
        W = sparse(W);
        W = kron(W, eye(T));
    end
    
    % Sort panel
    id = id(idx);
    time = time(idx);
    y = y(idx,:);
    X = X(idx,:);
    % Instruments are sorted later
    
    % Get endogenous and exogenous ariables 
    endog = unique(options.endog);
    exog = 1:size(X,2);
    exog(ismember(exog,endog)) = [];   
    
    % Get and check number of instruments
    lself = size(X(:,exog),2);
    Z = options.inst;
    if ~isempty(Z)
        Z = Z(idx,:);   % Sort instruments
    end
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
    
    
    % Step 1.1: 
    % And unique if there is no spatial correlation in the error term
    % Compute regression regression
    % Compute ivpanel regression
    
    % Build matrix of instruments (exlucing exogenous X)
    if endog
        H = [W*X(:,exnoslagX), W^2*X(:,exog), W^3*X(:,slagX), Z];
        estTemp = ivpanel(id, time, y, Xspa, H, method,'endog',endog);
    else
        estTemp = panel(id, time, y, Xspa, method);
    end
    
    % If spatial lag in the error term
    slagerror = options.slagerror;
    if slagerror
        
        % Check if the W matrix has size n (and not n*T)
        if size(Worig,1) ~= n
            error('The W matrix must have size n if the model has a spatial lag error');
        end

        % Get spatial matrix for the errors
        M = W;
        Morig = Worig;
        
        % Get n and T
        n = estTemp.n;
        T = estTemp.T;
        
        % Build full matrix of instruments (including exogenous X)
        if slagy
            Hw1 = [X(:,exog) W*X(:,exnoslagX), W^2*X(:,exog), W^3*X(:,slagX), Z];
            H = [X(:,exog) W*X(:,exnoslagX), W^2*X(:,exog), W^3*X(:,slagX), Z, ones(N,1)];
        else
            Hw1 = [X(:,exog), Z];
            H = [X(:,exog), Z, ones(N,1)];
        end
        
        % SEM with Fixed effects and no additional exogenous covariates
        switch method
            case 'fe'
                % Panel (IV if endog)
                if endog
                    estTemp = ivpanel(id, time, y, Xspa, [], method,'endog',endog, 'constant', 0,'useinstruments',Hw1);
                else
                    estTemp = panel(id, time, y, Xspa, method);
                end
                res = estTemp.res;   

                % Builg G and g (Kapoor et al 2007) (For Fixed Effects)
                u = res;
                ub = M * res;
                ubb = M * ub;

                G11 =  2/(n*(T-1)) * (u'*ub);
                G12 = -1/(n*(T-1)) * (ub'*ub);
                G13 =  1;
                G21 =  2/(n*(T-1)) * (ubb' * ub);
                G22 = -1/(n*(T-1)) * (ubb' * ubb);
                G23 =  1/n * trace(Morig' * Morig);
                G31 =  1/(n*(T-1)) * (u' * ubb + ub' * ub);
                G32 = -1/(n*(T-1)) * (ub' * ubb);
                G33 = 0;

                G = [G11, G12, G13; G21, G22, G23; G31, G32, G33];

                g1 = 1/(n*(T-1)) * (u'*u);
                g2 = 1/(n*(T-1)) * (ub'*ub);
                g3 = 1/(n*(T-1)) * (u'*ub);

                g = [g1; g2; g3];
             

                % Initial values for the solver ("pars" option in R default to
                % TRUE)
                wres = M*res;

                x0(1) = (res'*res)\(res'*wres);
                x0(2) = (res'*res)/(n*T); % Different if endogenous

                % Objective function
                perrfnc = @(x) (G * [x(1); x(1)^2; x(2) ] - g)' * (G * [x(1); x(1)^2; x(2) ] - g);
                options = optimset('TolFun',1e-8,'TolX',1e-8);
                rhos = fminsearch(perrfnc,x0,options);

                % Chrocrane-Orcut transformation
                rho = rhos(1);
                sigma2_srho = rhos(2);

                ytr = y - rho * M * y;
                Xtr = Xspa - rho * M * Xspa;

                % Estimate the transformed model (IV if endog)
                % est = panel(id, time, ytr, Xtr, method);
                % Manuel within regression for correct standard errors
                if endog                    
                    if slagy
                        % Add te W*X and W^2*X instruments
                        Z = [Z W*X(:,exnoslagX), W^2*X(:,exog), W^3*X(:,slagX),];
                    end
                    
                    estTemp = iv2sls(ytr - groupmeans(id,ytr,'replicate',1), Xtr - groupmeans(id,Xtr,'replicate',1), [Z - groupmeans(id,Z,'replicate',1)],'constant',0,'endog',endog);
                    sigma_v = sigma2_srho;
                    varcoef = (1/(n*T) * (estTemp.Xhat'*estTemp.Xhat)/sigma_v)\eye(size(Xtr,2)) / (n*T);
                else
                    estTemp = ols(ytr - groupmeans(id,ytr,'replicate',1), Xtr - groupmeans(id,Xtr,'replicate',1),'constant',0);
                    varcoef = estTemp.varcoef;
                end
                
                stderr = sqrt(diag(varcoef));
                sigma_srho = sqrt(sigma2_srho);
                    
            case 're'
                % OLS (IV if endog)
                if endog
                    estTemp1 = ivpanel(id,time,y,Xspa,[],'fe','endog',endog','constant',0,'useinstruments',Hw1);
                    estTemp2 = ivpanel(id,time,y,Xspa,[],'be','endog',endog','constant',0,'useinstruments',H);
                    res = estTemp1.res;
                    res2 = estTemp2.res;
                else
                    estTemp = ols(y, Xspa);    
                    res = estTemp.res; 
                end
                

                % Builg G and g (Piras 2013) (For Random Effects)
                u = res;
                ub = M * res;
                ubb = M * ub;
                
                Qu = u - groupmeans(id,u,'replicate',1);
                Qub = ub - groupmeans(id,ub,'replicate',1);
                Qubb = ubb - groupmeans(id,ubb,'replicate',1);
                
                G11 =  2/(n*(T-1)) * (u'*Qub);
                G12 = -1/(n*(T-1)) * (ub'*Qub);
                G13 =  1;
                G21 =  2/(n*(T-1)) * (ubb' * Qub);
                G22 = -1/(n*(T-1)) * (ubb' * Qubb);
                G23 =  1/n * trace(Morig' * Morig);
                G31 =  1/(n*(T-1)) * (u' * Qubb + ub' * Qub);
                G32 = -1/(n*(T-1)) * (ub' * Qubb);
                G33 = 0;

                G = [G11, G12, G13; G21, G22, G23; G31, G32, G33];

                g1 = 1/(n*(T-1)) * (u'*Qu);
                g2 = 1/(n*(T-1)) * (ub'*Qub);
                g3 = 1/(n*(T-1)) * (u'*Qub);

                g = [g1; g2; g3];

                % Initial values for the solver ("pars" option in R default to
                % TRUE)
                wres = M*res;

                x0(1) = (res'*res)\(res'*wres);
                x0(2) = (res'*res)/(n*T); % Different if endogenous

                % Objective function
                perrfnc = @(x) (G * [x(1); x(1)^2; x(2) ] - g)' * (G * [x(1); x(1)^2; x(2) ] - g);
                options = optimset('TolFun',1e-8,'TolX',1e-8);
                rhos = fminsearch(perrfnc,x0,options);                

                % Chrocrane-Orcut transformation
                rho = rhos(1);
                sigma2_srho = rhos(2);

                ytr = y - rho * M * y;
                    % Add 1 for the constant term
                Xtr = [Xspa ones(n*T,1)] - rho * M * [Xspa ones(n*T,1)];
                
                % Compue theta (different if IV)
                if endog
                    % Apply the sqrt(T) adjustment to the residuals
                    res2 = res2 * sqrt(T);
                    sigma2_1 = n^-1 * (res2 - rho * Morig * res2)' *  (res2 - rho * Morig * res2);
                else                    
                    sigma2_1 = n^-1 * (res - rho * M* res)' *( groupmeans(id,u,'replicate',1) - rho * groupmeans(id,ub,'replicate',1));
                end
                sigma2_v =  sigma2_srho;
                theta = 1 - sqrt(sigma2_v/sigma2_1);
                
                % RE transformation
                ytr = ytr - theta * groupmeans(id,ytr,'replicate',1);
                Xtr = Xtr - theta * groupmeans(id,Xtr,'replicate',1);
                if endog
                    if slagy
                        % Add te W*X and W^2*X instruments
                        Z = [Z X(:,exog) W*X(:,exnoslagX), W^2*X(:,exog), W^3*X(:,slagX)];
                        Z = [groupmeans(id,Z,'replicate',1),Z-groupmeans(id,Z,'replicate',1),ones(n*T,1)];
                        estTemp = iv2sls(ytr,Xtr,[],'endog',endog,'constant',0,'useinstruments',Z);
                    else
                        estTemp = iv2sls(ytr,Xtr,[groupmeans(id,Z,'replicate',1),Z-groupmeans(id,Z,'replicate',1)],'endog',endog,'constant',0);
                    end
                    
                    
                    % Correct estTempimation of covariance matrix
                    varcoef = (1/(n*T) * (estTemp.Xhat'*estTemp.Xhat)/sigma2_v)\eye(size(Xtr,2)) / (n*T);
                else
                    estTemp = ols(ytr,Xtr,'constant',0);
                    
                    % Correct estimation of covariance matrix
                    varcoef = (1/(n*T) * (Xtr'*Xtr)/sigma2_v)\eye(size(Xtr,2)) / (n*T);
                end                
                
                % Correct standard errors
               stderr = sqrt(diag(varcoef));    
               sigma_srho = sqrt(sigma2_srho);
                
        end
                
        
    else
        rho = NaN;
        sigma_srho = NaN;
        stderr = estTemp.stderr;
        varcoef = estTemp.varcoef;
        if strcmp(method,'re') || strcmp(method,'ec')
            theta = estTemp.theta;
            % sigma2_mu = estTemp.sigma2_mu;
            sigma2_1 = estTemp.sigma2_1;
            sigma2_v = estTemp.sigma2_v;
            % rho_mu = estTemp.rho_mu;
        end
    end
        
    % Number of regressors (+1 for the constant if not FE)
    hasConstant = ~strcmpi(method,'fe');
    k = size(Xspa,2) + hasConstant;
    if hasConstant
        Xspa = [Xspa, ones(N,1)];
    end
    
    % Coefficients, fitted values and residuals
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
    est.id = id;
    est.uid = unique(id);
    est.time = time;
    est.idx = idx;
    est.y = y;
    est.X = X;
    est.W = Worig;
    
    % Save estimation
    est.method = strcat(upper(method),'S2SLS');
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
    
    est.Tid = Tid;
    est.Tmean = Tmean;
    est.Thmean = Thmean;
    if strcmp(method,'re') || strcmp(method,'ec')
        est.theta = theta;
        % est.sigma2_mu = sigma2_mu;
        est.sigma2_1 = sigma2_1;
        est.sigma2_v = sigma2_v;
        % est.rho_mu = rho_mu;
    end
    
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




% CODE from OLD full function
    
    %{
    % Get number of observations
    N = size(y,1);
    
    % Get uniques
    id_uniq = unique(id);
    time_uniq = unique(time);
    idt = [id time];   
    
    % Check no id and time is repeated
    if size(idt,1) ~= size(unique(idt,'rows'))
        error('Multiple observations with same id and time');
    end    
    
    % Sort panel
    [~, idx] = sortrows(idt,[1 2]);
    id = id(idx);
    time = time(idx);
    y = y(idx,:);
    X = X(idx,:);
    
    % Get number of units
    n = length(id_uniq);
    T = length(time_uniq);
    
    % Check balanced
    if N == (n * T)
        isBalanced = 1;
        Tid = (T*ones(n,1));
        Tmean = T;
        Thmean = T;
    else
        isBalanced = 0;
        % Compute number of time periods per ID
        Tid = nan(n,1);
        for i=1:n
            Tid(i,1) = sum(id==i);
        end
        % Compute mean and harmonic mean
        Tmean = mean(Tid);
        Thmean = n ./ sum(1./Tid);
    end
    
    % Store original data (shorted)
    est.id = id;
    est.uid = unique(id);
    est.time = time;
    est.idx = idx;
    est.y = y;
    est.X = X;
    est.W = W;
    
    % Usefull matrices
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
    
    % Build matrix of instruments
    Z = [X, W*X, W^2*X];
        
    % Number of varaibles and hasConstant
    k = size(X,2);
    if ~strcmp(method,'fe')
        k = k + 1;  % For the constant term in non fe estiamtions
        hasConstant = 1;
    else
        hasConstant = 0;
    end
    
    % Add the lagged Y to the number of variables
    k = k + 1;
    
    % Transform variables
    switch method
        case 'po'
            % Add spatial lag and constant term
            X = [X W*y ones(N,1)];
            Z = [Z ones(N,1)];
            
            % No transformation por pool
            ytr = y;
            Xtr = X;
            Ztr = Z;
            
            resdf = N;
        case 'fe'      
            % Add spatial lag
            X = [X W*y];
            
            % Fixed effects transformation
            ytr = Q*y;
            Xtr = Q*X;
            Ztr = Q*Z;
            
            resdf = (N-n)-k; % (N-n) instead of (n*T) for compatibility with unbalanced panels

        case 'be'
            % Add spatial lag and constant term
            X = [X W*y ones(N,1)];
            Z = [Z ones(N,1)];
            
            % Between transformation.
            ytr = groupmeans(id, y);
            Xtr = groupmeans(id, X);
            Ztr = groupmeans(id, Z);
            
            resdf = n-k;
        case {'re','ec'}
            % Get time invariant variables to remove from the first 'fe'
            % and 'be' regressions
            tinvariant = istinvariant(id,X);
            
            % Get sigma2_v from 'fe' residuals
            estfe = spanel(id,time,y,X(:,~tinvariant),W,'fe');
            sigma2_v = estfe.resvar;
            
            % Get sigma2_1 from 'be' residuals
            estbe = spanel(id,time,y,X(:,~tinvariant),W,'be');
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
                theta = 1 - sqrt(sigma2_v./(Tid_g.*sigma2_mu + sigma2_v));
            end
            
            % If Baltagi's error components EC2SLS
            if strcmp(method,'ec')
                Z = [Q*Z, P*Z];
            end
            
            % Add spatial lag and constant term
            X = [X W*y ones(N,1)];
            Z = [Z ones(N,1)];
            
            % Random effects transformation
            if isBalanced
                Qre = eye(N) - theta * P;
            else
                Qre = eye(N) - repmat(theta,1,N) .* P;
            end
            ytr = Qre * y;
            Xtr = Qre * X;   
            Ztr = Qre * Z;
            
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
        M0 = eye(n) - 1/n * ones(n);
        adjr2_correction = (n - 1) ./ (resdf);
    else
        M0 = eye(N) - 1/N * ones(N);
        adjr2_correction = (N - 1) ./ (resdf);
    end
    RSS = res' * res;
    TSS = y' * y;
    ESS = TSS - RSS;
    r2 = 1 - (res' * M0 * res) ./ (ytr' * M0 * ytr);
    adjr2 = 1 - adjr2_correction .* (1 - r2);
    
    
    % Save estimation
    est.method = upper(method);
    est.options = options;
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
    est.isInstrumental = 0;
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
    %}
    



