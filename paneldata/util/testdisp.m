function [  ] = testdisp( test )
%TESTDISP Display test results
%   TESTDISP( test ) Display test results stored in a 'testout' structure
%
%   Example
%     
%      testdisp(test);
%
%   See also TESTOUT, ESTDISP
%
%   Copyright 2013-2015 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 29, February, 2016
%

    switch(test.test)
        case 'WALD'
            testdisp_WALD( test )
        case 'HAUSMAN'
            testdisp_HAUSMAN( test )
        case 'HREGENDOG'
            testdisp_HREGENDOG( test )
        case 'WUENDOG'
            testdisp_WUENDOG( test )
        case 'SARGANOI'
            testdisp_SARGANOI( test )
        case 'RESET'
            testdisp_RESET( test )
        case 'BPHET'
            testdisp_BPHET( test ) 
        case 'WHITEHET'
            testdisp_WHITEHET( test )
        % Panel Specific
        case 'EFFECTSF'
            testdisp_EFFECTSF( test )
        case 'BPRE'
            testdisp_BPRE( test )
        case 'MUNDLAKVA'
            testdisp_MUNDLAKVA( test )
        case 'POOL'
            testdisp_POOL( test );
        case 'WOOLSERIAL'
            testdisp_WOOLSERIAL( test );
        case 'BLSERIAL'
            testdisp_BLSERIAL( test );
        case 'PESARANCSD'
            testdisp_PESARANCSD( test );
        % Spatial specific
        case 'BSJKSA'
            testdisp_BSJKSA( test );
        case 'SHAUSMAN'
            testdisp_SHAUSMAN( test )
        otherwise
            error('testdisp not implemented for this test');
    end

end

function [] = testdisp_WALD( test )

    if test.Randr == 0
        error('testdisp only available for R and r test waldsigtest. Use estdisp to display individual significance.');
    end
    
    % Display input
    fprintf('_________________________________________________________\n');
    disp('<strong>Wald joint significance test</strong>');
    %{
    disp('R = ');
    disp(test.R);
    disp('r = ');
    disp(test.r);
    %}
    disp(' ');
    
    % Display output
    if test.isAsymptotic
        fprintf('Chi2(%d) = %f \n', test.df, test.value);
        fprintf('p-value = %5.4f \n', test.p);
    else
        fprintf('F(%d, %d) = %f \n', test.df(1), test.df(2), test.value);
        fprintf('p-value = %5.4f \n', test.p);
    end
    
    disp(' ');

end

function [] = testdisp_HAUSMAN( test )
    
    ktest = test.df; % The number of coefficients tested are the df
    
    estAname = test.estA.method;
    estBname = test.estB.method;
    
    stderrdif = sqrt(diag(test.estA.varcoef(1:ktest,1:ktest) - test.estB.varcoef(1:ktest,1:ktest)));

    fprintf('_________________________________________________________\n');
    fprintf('<strong>Hausman''s test of specification</strong>\n')
    disp(' ');
    
    fprintf('------------------------------------------------------------------------\n');
    fprintf('        Varname | %13s %13s    Coef. Diff    S.E. Diff  \n', strcat('A:',estAname), strcat('B:',estBname));

    fprintf('------------------------------------------------------------------------\n');
    for i=1:1:ktest
        %fprintf('%7.7s: %11f (%10f)',char(varnames(i)),pestre.coef(i),sqrt(pestre.varcoef(i,i)));
        fprintf('%15.15s | %13.6f %13.6f %13.6f %12.6f\n',char(test.estA.xnames(i)),test.estA.coef(i), test.estB.coef(i),test.estA.coef(i)-test.estB.coef(i), stderrdif(i));
    end
    fprintf('------------------------------------------------------------------------\n');
    % disp(' ');
    fprintf('A is consistent under H0 and H1 (A = %s)\n',estAname)
    fprintf('B is consistent under H0        (B = %s)\n',estBname);
    fprintf('H0: coef(A) - coef(B)  = 0 \n');
    fprintf('H1: coef(A) - coef(B) != 0 \n');
    %disp(' ');
    fprintf('      H = %f ~ Chi2(%d)\n',test.value,test.df);
    if test.value > 0
        fprintf('p-value = %5.4f \n',test.p);
    else
        fprintf('WARNING: Invalid estimates. H <= 0');
    end
    
    disp(' ');

end

function [  ] = testdisp_HREGENDOG( test )
    
    fprintf('_________________________________________________________\n');
    fprintf('<strong>Hausman''s regression test of endogeneity</strong>\n')
    disp(' ');
    
    fprintf('Test robust to heteroskedasticity \n');
    
    fprintf('H0: Variables are exogenous \n');
    
    fprintf('F(%d,%d) = %f \n',test.df(1),test.df(2),test.value);
    fprintf('p-value = %5.4f \n',test.p);

    disp(' ');
    
end

function [  ] = testdisp_WUENDOG( test )
    
    fprintf('_________________________________________________________\n');
    fprintf('<strong>Wu''s variable addition test of endogeneity</strong>\n')
    disp(' ');
    
    fprintf('Test robust to heteroskedasticity \n');
    
    fprintf('H0: Variables are exogenous \n');
    
    fprintf('F(%d,%d) = %f \n',test.df(1),test.df(2),test.value);
    fprintf('p-value = %5.4f \n',test.p);

    disp(' ');
    
end


function [  ] = testdisp_SARGANOI( test )
    
    fprintf('_________________________________________________________\n');
    fprintf('<strong>Sargan''s test of overidentification</strong>\n')
    disp(' ');
    
    if test.isRobust
        %error('Not yet implemented. Wooldridge');
        disp('Wooldridge''s (1995) robust score test');
    end
    
    fprintf('H0: Instruments are uncorrelated with the error term \n');
    
    fprintf('  Score = %f ~ Chi2(%d)\n',test.value, test.df);
    fprintf('p-value = %5.4f \n',test.p);

    disp(' ');
    
end

function [ ] = testdisp_RESET( test )

    fprintf('_________________________________________________________\n');
    fprintf('<strong>Ramsey''s Regression Equation Specification Error Test (RESET)</strong>\n')
    disp(' ');
    
    fprintf('Regression of y on X ');    
    for i=2:test.powers
        fprintf('yhat.^%d ',i)
    end
    fprintf(' and 1 \n');
    
    fprintf('H0: Powers of yhat are not jointly significant \n');
    fprintf('H1: The model is mis-specified \n');
    
    fprintf('F(%d,%d) = %f \n',test.df(1),test.df(2),test.value);
    fprintf('p-value = %5.4f \n',test.p);

    disp(' ');

end

function [ ] = testdisp_BPHET( test )

    fprintf('_________________________________________________________\n');
    fprintf('<strong>Breusch-Pagan''s test of Heteroskedasticity</strong>\n')
    disp(' ');

    disp('Koenker''s (1981) version without the normality assumption');
    if ~test.isInstrumental
        disp('Regression of u.^2 on X and 1 '); 
    else
        disp('Regression of u.^2 on exogenous X, Z, and 1 '); 
    end
    
    fprintf('H0: Homoskedasticity\n');
    
    fprintf('  Score = %f ~ Chi2(%d)\n',test.value,test.df);
    fprintf('p-value = %5.4f \n',test.p);

    disp(' ');
    
end

function [] = testdisp_WHITEHET( test )

    fprintf('_________________________________________________________\n');
    fprintf('<strong>White''s test of Heteroskedasticity</strong>\n')
    disp(' ');
    
    if ~test.isInstrumental
        disp('Regression of u.^2 on X X.^2 cross-products and 1'); 
    else
        disp('Regression of u.^2 on exogenous X X.^2 Z Z.^2 cross-products and 1');
    end
    if test.duprm > 0
        fprintf('%d duplicates removed \n', test.duprm)
    end
    
    fprintf('H0: Homoskedasticity\n');
    
    fprintf('  Score = %f ~ Chi2(%d)\n',test.value,test.df);
    fprintf('p-value = %5.4f \n',test.p);

    disp(' ');

end

function [] = testdisp_EFFECTSF( test )


    fprintf('_________________________________________________________\n');
    fprintf('<strong>F test of individual effects</strong>\n')
    disp(' ');
    
    fprintf('H0: All mu_i = 0\n');
    
    fprintf('F(%d,%d) = %f \n',test.df(1), test.df(2), test.value);
    fprintf('p-value = %5.4f \n',test.p);

    disp(' ');
    
end


function [] = testdisp_BPRE( test )


    fprintf('_________________________________________________________\n');
    fprintf('<strong>Breusch-Pagan''s LM test for random effects</strong>\n')
    disp(' ');
    
    disp('Baltagi and Li (1990) version of the Breusch and Pagan (1980) test'); 
    
    fprintf('H0: sigma2_mu = 0\n');
    
    fprintf('  LM = %f ~ Chi2(%d)\n',test.value,test.df);
    fprintf('p-value = %5.4f \n',test.p);

    disp(' ');

end


function [  ] = testdisp_MUNDLAKVA( test )
    
    fprintf('_________________________________________________________\n');
    fprintf('<strong>Mundlak''s variable addition test for fixed or random effects</strong>\n')
    disp(' ');
    
    if test.isRobust
        fprintf('Test robust to heteroskedasticity \n');
    end
    
    fprintf('H0: Group means are zero. Random effects. \n');
    
    fprintf('Chi2(%d) = %f \n',test.df,test.value);
    fprintf(' p-value = %5.4f \n',test.p);

    disp(' ');
    
end

function [] = testdisp_POOL( test )


    fprintf('_________________________________________________________\n');
    fprintf('<strong>Test of poolability</strong>\n')
    disp(' ');
    
    fprintf('H0: Stability of coefficients\n');
    
    fprintf('F(%d,%d) = %f \n',test.df(1), test.df(2),test.value);
    fprintf('p-value = %5.4f \n',test.p);

    disp(' ');

end


function [] = testdisp_WOOLSERIAL( test )


    fprintf('_________________________________________________________\n');
    fprintf('<strong>Wooldridge''s test for serial correlation</strong>\n')
    disp(' ');
      
    fprintf('H0: Corr(res_{T-1}, res_T) = rho. No serial correlation\n');
    fprintf('rho = -1/(T-1) = %f \n',test.rho);
    
    fprintf('F(%d,%d) = %f \n',test.df(1), test.df(2),test.value);
    fprintf('p-value = %5.4f \n',test.p);

    disp(' ');

end

function [] = testdisp_BLSERIAL( test )


    fprintf('_________________________________________________________\n');
    fprintf('<strong>Baltagi and Li''s test for serial correlation and random effects</strong>\n')
    disp(' ');

    fprintf('H0: No random effects and no serial correlation.\n');
    fprintf('H1: Random effects or serial correlation.\n')

    fprintf('Chi2(%d) = %f \n',test.df, test.value);
    fprintf('p-value = %5.4f \n',test.p);

    disp(' ');

end

function [] = testdisp_PESARANCSD( test )


    fprintf('_________________________________________________________\n');
    fprintf('<strong>Pesaran''s test of cross sectional dependence</strong>\n')
    disp(' ');
    
    %disp('Baltagi and Li (1990) version of the Breusch and Pagan (1980) test'); 
    
    fprintf('H0: Corr(res_{it}, res_{jt}) = 0 for i != j\n');
    
    fprintf('     CD = %f \n',test.value);
    fprintf('p-value = %5.4f \n',test.p);

    disp(' ');

end

function [] = testdisp_BSJKSA( test )


    fprintf('_________________________________________________________\n');
    fprintf('<strong>Baltagi, Song, Jung and Koh''s test for serial correlation,\n  spatial autocorrelation and random effects</strong>\n')
    disp(' ');

    fprintf('H0: No spatial autocorrelation, no serial error correlation and no re.\n');
    fprintf('H1: Spatial autocorrelation or serial error correaltion or random effects.\n')

    fprintf('Chi2(%d) = %f \n',test.df, test.value);
    fprintf('p-value = %5.4f \n',test.p);

    disp(' ');

end


function [] = testdisp_SHAUSMAN( test )
    
    ktest = test.df; % The number of coefficients tested are the df
    
    estAname = test.estA.method;
    estBname = test.estB.method;
    
    stderrdif = sqrt(diag(test.estA.varcoef(1:ktest,1:ktest) - test.estB.varcoef(1:ktest,1:ktest)));

    fprintf('_________________________________________________________\n');
    fprintf('<strong>Spatial Hausman''s test</strong>\n')
    disp(' ');
    
    disp('Mutl and Pfaffermayr (200/');
    
    fprintf('------------------------------------------------------------------------\n');
    fprintf('        Varname | %13s %13s    Coef. Diff    S.E. Diff  \n', strcat('A:',estAname), strcat('B:',estBname));

    fprintf('------------------------------------------------------------------------\n');
    for i=1:1:ktest
        %fprintf('%7.7s: %11f (%10f)',char(varnames(i)),pestre.coef(i),sqrt(pestre.varcoef(i,i)));
        fprintf('%15.15s | %13.6f %13.6f %13.6f %12.6f\n',char(test.estA.xnames(i)),test.estA.coef(i), test.estB.coef(i),test.estA.coef(i)-test.estB.coef(i), stderrdif(i));
    end
    fprintf('------------------------------------------------------------------------\n');
    % disp(' ');
    fprintf('A is consistent under H0 and H1 (A = %s)\n',estAname)
    fprintf('B is consistent under H0        (B = %s)\n',estBname);
    fprintf('H0: coef(A) - coef(B)  = 0 \n');
    fprintf('H1: coef(A) - coef(B) != 0 \n');
    %disp(' ');
    fprintf('      H = %f ~ Chi2(%d)\n',test.value,test.df);
    if test.value > 0
        fprintf('p-value = %5.4f \n',test.p);
    else
        fprintf('WARNING: Invalid estimates. H <= 0');
    end
    
    disp(' ');

end