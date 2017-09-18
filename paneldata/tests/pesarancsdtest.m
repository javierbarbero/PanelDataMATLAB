function [ test ] = pesarancsdtest( est )
%PESARANCSDTEST Pesaran cross-sectional dependence test
%   Computes the Pesaran cross-sectional dependence test
%
%   testo = PESARANCSDTEST( est ) Computes the Pesaran cross-sectional 
%   dependence test for the specified estimation output. Returns a test 
%   output structure, testout.
%
%   Example
%     
%      test = pesarancsdtest(est);
%
%   See also TESTOUT
%
%   Copyright Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 11, August, 2017
%

    if est.isMultiEq
        error('Pesaran''s CSD test not available for multi-equation models')
    end
    
    if strcmpi(est.options.method,'be')
        error('Pesaran''s CSD test can not be performed after between estimation');
    end
    
    if est.isInstrumental
        error('Pesaran''s CSD test not valid for instrumetnal estimation')
    end
    
    % Create otuput structure
    test = testout();
    
    % Get used data
    id = est.id;
    uid = est.uid;
    time = est.time;
    res = est.res;
    T = est.T;
    n = est.n;
    
    % Compute CD statistic
    if est.isBalanced
        total_sum = 0;
        for i=1:1:n-1
            for j=i+1:1:n
                temp = corrcoef(res(id == uid(i)),res(id == uid(j)) );
                total_sum = total_sum + temp(1,2);
            end
        end
        CD = sqrt(2*T/(n*(n-1))) * total_sum;
    else

        % Unbalanced
        total_sum = 0;
        for i=1:1:n-1
            for j=i+1:1:n
                % Get common values
                Tij = intersect(time(id == uid(i)), time(id == uid(j)) );
                
                % Display error message if less than 2 common time periods
                if length(Tij) < 2
                    error('Some units have 1 or less common time periods. Test cannot be performed.')
                end
                
                % Get common time residuals for both groups
                resi = res(id == uid(i) & ismember(time, Tij));
                resj = res(id == uid(j) & ismember(time, Tij));
                
                % Compute correlation coefficient
                temp = corrcoef(resi, resj);
                
                % Add to the total sum corrected by the number of common
                % observations
                total_sum = total_sum + sqrt(length(Tij)) * temp(1,2);
            end
        end
        CD = sqrt(2/(n*(n-1))) * total_sum;
    end
    df = NaN;
    p = 1 - normaldist(abs(CD));

    
    % Store results
    test.test = 'PESARANCSD';
    test.value = CD;
    test.p = p;
    test.df = df;
    test.isAsymptotic = 0;
    test.isRobust = 0;
    
end

