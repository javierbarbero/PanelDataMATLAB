function [ test ] = waldsigtest( est, varargin )
%WALDSIGTEST Wald significance test for linear hypothesis.
%   Computes the Wald significance test for linear hypothesis.
%
%   testo = WALDSIGTEST( est ) Computes the Wald significance test for 
%   linear hypothesis for the specified estimation output. Returns a test 
%   output structure, testout.
%   %   testo = WALDSIGTEST( est ) Computes the Wald significance test for 
%   linear hypothesis for the specified estimation output, uwing the R and 
%   r matrices. Returns a test output structure, testout.
%
%   Example
%
%      test = waldsigtest(est);
%      test = waldsigtest(est, R, r);
%
%   See also TESTOUT
%
%   Copyright 2013-2015 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 17, June, 2015
%

    % Create otuput structure
    test = testout();

    % Individual significance Wald test
    if nargin == 1
        % Test is not R and r
        Randr = 0;
                
        if est.isAsymptotic
            % Normal distribution
            Wdist = est.coef ./ est.stderr;
            df = NaN;
            p = (1 - normaldist(abs(Wdist)))' * 2;
        else
            % Student t
            Wdist = est.coef ./ est.stderr;
            df = est.resdf;
            p = (1 - tdist(abs(Wdist), df))' * 2;
        end   
        
    % R and r Wald test
    elseif nargin == 3
        % Test is R and r
        Randr = 1;
        
        R = varargin{1};
        r = varargin{2};
        % Check if size is correct
        if size(R,1) ~= size(r,1)
            error('Incorrect size of the R and r test matrices')
        end
        
        q = size(R,1);       
        
        % Get coefficients and var coef
        coef = est.coef;
        varcoef = est.varcoef;
        resdf = est.resdf;       
        
        if est.isAsymptotic
            % Chi2 distribution
            Wdist = (R*coef - r)' * ((R*varcoef*R')\eye(q)) *(R*coef - r);
            df = q;
            p = 1 - chi2dist(Wdist, df);
        else     
            % F distribution
            Wdist = (R*coef - r)' * ((R*varcoef*R')\eye(q)) *(R*coef - r) / q; 
            df = [q, resdf];
            p = 1 - fdist(Wdist, df(1), df(2));
        end
    else
        error('Incorrect number of arguments')
    end
    
    % Save test results
    test.test = 'WALD';
    test.value = Wdist;
    test.p = p;
    test.df = df;
    test.isAsymptotic = est.isAsymptotic;
    test.isRobust = est.isRobust;
    
    % waldsigtest specific
    test.Randr = Randr;
    if Randr == 1
        test.R = R;
        test.r = r;
    end

end

