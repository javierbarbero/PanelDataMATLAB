clear all
clc

% Tolerance
tolCoef = 1e-5;
tolTestStrong = 1e-2;
tolTestLow = 1e-1;
tolTest = 1e-4;

% Load data
load '..\data\CigarData'

y = log(c);
X = [log(price), log(ndi), log(pimin)];
endog = 1;
Z = [log(ndi_1), log(pimin_1)];
id = state;

% Pool
iv = iv2sls(y,X,Z,'endog',endog);
regpo = ivpanel(id,year,y,X,Z,'po','endog',endog);

assert(all(abs(regpo.coef - iv.coef) <= tolCoef),'Estimated coefficients')
assert(all(abs(regpo.stderr - iv.stderr) <= tolCoef),'Standard errors')
assert(abs(regpo.r2 - iv.r2) <= tolCoef,'R2')
assert(abs(regpo.RSS -  iv.RSS) <= tolCoef,'RSS')
assert(abs(regpo.resdf - iv.resdf) <= tolCoef,'res DF')

% FE
regivfe = ivpanel(id,year,y,X,Z,'fe','endog',endog);

realcoeff = [-1.0163549; 0.53784801; 0.31237215];
assert(all(abs(regivfe.coef - realcoeff) <= tolCoef),'Estimated coefficients')

readstderr = [0.24919655; 0.02303342; 0.22839452];
assert(all(abs(regivfe.stderr - readstderr) <= tolCoef),'Standard errors')

assert(abs(regivfe.r2 - 0.6406418221500125) <= tolCoef,'R2')
assert(abs(regivfe.RSS -  7.731113593920569) <= tolCoef,'RSS')
assert(abs(regivfe.resdf - 1285) <= tolCoef,'res DF')

R = [eye(regivfe.k)];
r = zeros(regivfe.k,1);
wt = waldsigtest(regivfe,R,r);
%assert(abs(wt.value - 5094492.537879358) <= tolTestStrong,'Wald joint significance')

% Overll individual effect
[oeff, soeff] = ieffects(regivfe,'overall');
assert(all(abs([oeff, soeff] - [ 2.991411,  0.08110602 ]) <= tolTestLow),'Wald joint significance (Robust)')

% Sargan over id
soi = sarganoitest(regivfe);
assert(abs(soi.value - 25.52015409617647) <= tolTestStrong,'SARGANOIT value')
assert(abs(soi.p -  4.37785894954e-07) <= tolTest,'SARGANOIT p')
assert(abs(soi.df -  1) <= tolTest,'SARGANOIT df')


% BE
regivbe = ivpanel(id,year,y,X,Z,'be','endog',endog);

realcoeff = [-3.2752235; 0.83220216; 1.1810666; 6.1738967];
assert(all(abs(regivbe.coef - realcoeff) <= tolCoef),'Estimated coefficients')

readstderr = [2.613918; 0.40038775; 1.3237463; 3.2967346];
assert(all(abs(regivbe.stderr - readstderr) <= tolCoef),'Standard errors')

assert(abs(regivbe.r2 - 0.0285289946745334) <= tolCoef,'R2')
assert(abs(regivbe.RSS -   1.526528740750127) <= tolCoef,'RSS')
assert(abs(regivbe.resdf - 42) <= tolCoef,'res DF')

R = [eye(regivbe.k-1) zeros(regivbe.k-1,1)];
r = zeros(regivbe.k-1,1);
wt = waldsigtest(regivbe,R,r);
assert(abs(wt.value - 6.660387374883474) <= tolTestStrong,'Wald joint significance')

% RE
regivre = ivpanel(id,year,y,X,Z,'re','endog',endog);

realcoeff = [-1.0071126; 0.53747325; 0.30356681; 2.9921215];
assert(all(abs(regivre.coef - realcoeff) <= tolCoef),'Estimated coefficients')

readstderr = [0.24735393; 0.0230282; 0.22643016; 0.08566819];
assert(all(abs(regivre.stderr - readstderr) <= tolCoef),'Standard errors')

assert(abs(regivre.r2 - 0.4170237931223892) <= tolCoef,'R2')
%assert(abs(regre.adjr2 - 0.9933422638044184) <= tolCoef,'Adj R2')
%assert(abs(regre.RSS -  0.2977007950047585) <= tolCoef,'RSS')
%assert(abs(regre.resdf - 43) <= tolCoef,'res DF')

R = [eye(regivre.k-1) zeros(regivre.k-1,1)];
r = zeros(regivre.k-1,1);
wt = waldsigtest(regivre,R,r);
assert(abs(wt.value - 1820.426690406653) <= tolTestStrong,'Wald joint significance')

% Sargan over id
% NOT IN STATA


% Error Components EC2SLS
regivec = ivpanel(id,year,y,X,Z,'ec','endog',endog);

realcoeff = [-0.99268084; 0.53641048; 0.2903893; 2.9951236];
assert(all(abs(regivec.coef - realcoeff) <= tolCoef),'Estimated coefficients')

readstderr = [0.23586856; 0.02235608; 0.21596955; 0.08419795];
assert(all(abs(regivec.stderr - readstderr) <= tolCoef),'Standard errors')

assert(abs(regivec.r2 - 0.4168560064893315) <= tolCoef,'R2')
%assert(abs(regre.adjr2 - 0.9933422638044184) <= tolCoef,'Adj R2')
%assert(abs(regre.RSS -  0.2977007950047585) <= tolCoef,'RSS')
%assert(abs(regre.resdf - 43) <= tolCoef,'res DF')

R = [eye(regivec.k-1) zeros(regivec.k-1,1)];
r = zeros(regivec.k-1,1);
wt = waldsigtest(regivec,R,r);
assert(abs(wt.value - 1825.252840874937) <= tolTestStrong,'Wald joint significance')

