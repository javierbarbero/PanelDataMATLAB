clear all
clc

% Tolerance
tolCoef = 1e-6;
tolTestStrong = 1e-2;
tolTestLow = 1e-1;
tolTest = 1e-4;

% Read data
load '..\data\MunnellData'
y = log(gsp);
X = [log(pcap), log(pc), log(emp), unemp];
endog = 1;
Z = [log(hwy), log(water)];

% IV2SLS estimation
iv1 = iv2sls(y,X,Z,'endog',endog);

realcoeff = [0.19696384; 0.29667868; 0.56665739; -0.00691241; 1.5608487];
assert(all(abs(iv1.coef - realcoeff) <= tolCoef),'Estimated coefficients')

readstderr = [0.02100863; 0.01089453; 0.01585087; 0.00141818; 0.06234629];
assert(all(abs(iv1.stderr - readstderr) <= tolCoef),'Standard errors')

assert(abs(iv1.r2 - 0.9925388114834532) <= tolCoef,'R2')
assert(abs(iv1.adjr2 - 0.9925020115400918) <= tolCoef,'Adj R2')
assert(abs(iv1.RSS - 6.340584258257985) <= tolCoef,'RSS')
assert(abs(iv1.resdf - 816) <= tolCoef,'res DF')

R = [eye(iv1.k-1) zeros(iv1.k-1,1)];
r = zeros(iv1.k-1,1);
wt = waldsigtest(iv1,R,r);
assert(abs(wt.value - 108562.2700678448) <= tolTestStrong,'Wald joint significance')

% Sargen over id
soi = sarganoitest(iv1);
assert(abs(soi.value - 43.50743715314226) <= tolTestStrong,'SARGANOIT value')
assert(abs(soi.p -  4.22351188425e-11) <= tolTest,'SARGANOIT p')
assert(abs(soi.df -  1) <= tolTest,'SARGANOIT df')

% Breush-Pagan Heteroskedasticity test
bpt = bphettest(iv1);
assert(abs(bpt.value - 91.8186454163849) <= tolTestStrong,'BPHET value')
assert(abs(bpt.p -  2.78715114613e-18) <= tolTest,'BPHET p')
assert(all(abs(bpt.df -  5) <= tolTest),'BPHET df')

wht = whitehettest(iv1);
assert(abs(wht.value - 284.6907936700811) <= tolTest,'WHIEHET value')
assert(abs(wht.p -  1.06820968342e-48) <= tolTest,'WHIEHET p')
assert(all(abs(wht.df -  20) <= tolTest),'WHIEHET df')

% Hausman test
ols1 = ols(y,X);
ht = hausmantest(iv1,ols1);
assert(abs(ht.value - 11.96631447455107) <= tolTest,'HAUSMANTEST value')
assert(abs(ht.p -  .0176035252884972) <= tolTest,'HAUSMANTEST p')
assert(all(abs(ht.df -  4) <= tolTest),'HAUSMANTEST df')

% Wu variable addition test
wu = wuendogtest(iv1);
assert(abs(wu.value - 12.17513391342789) <= tolTest,'WUENDOG value')
assert(abs(wu.p -  .0005104520294718) <= tolTest,'WUENDOG p')
assert(all(abs(wu.df -  [1, 810]) <= tolTest),'WUENDOG df')

% IV2SLS robust estimation
iv1r = iv2sls(y,X,Z,'endog',endog,'vartype','robust');

readstderr = [0.02141369;0.01293101;0.02084542;0.00134565;0.07360416];
assert(all(abs(iv1r.stderr - readstderr) <= tolCoef),'Standard errors (Robust)')

R = [eye(iv1r.k-1) zeros(iv1r.k-1,1)];
r = zeros(iv1r.k-1,1);
wtr = waldsigtest(iv1r,R,r);
assert(abs(wtr.value -  134804.1371646649) <= tolTestLow,'Wald joint significance (Robust)')

% Sargen over id (Robust)
soir = sarganoitest(iv1r);
assert(abs(soir.value -  58.488295601461) <= tolTestStrong,'SARGANOIT value (Robust)')
assert(abs(soir.p - 2.04503945461e-14) <= tolTest,'SARGANOIT p (Robust)')
assert(abs(soir.df -  1) <= tolTest,'SARGANOIT df (Robust)')

% Wu variable addition test (Robust)
wur = wuendogtest(iv1r);
assert(abs(wur.value - 17.35936415554554) <= tolTest,'WUENDOG value')
assert(abs(wur.p -  .0000342623535047) <= tolTest,'WUENDOG p')
assert(all(abs(wur.df -  [1, 810]) <= tolTest),'WUENDOG df')

