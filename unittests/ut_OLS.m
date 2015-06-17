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

% OLS estimation
ols1 = ols(y,X);

realcoeff = [0.15500703; 0.30919016; 0.59393489 ; -0.00673298; 1.6433022];
assert(all(abs(ols1.coef - realcoeff) <= tolCoef),'Estimated coefficients')

readstderr = [0.01715377; 0.01027199; 0.01374746; 0.00141638; 0.05758725];
assert(all(abs(ols1.stderr - readstderr) <= tolCoef),'Standard errors')

assert(abs(ols1.r2 - 0.9925934477502721) <= tolCoef,'R2')
assert(abs(ols1.adjr2 - 0.9925569172829492) <= tolCoef,'Adj R2')
assert(abs(ols1.RSS - 6.294153873533901) <= tolCoef,'RSS')
assert(abs(ols1.resdf - 811) <= tolCoef,'res DF')

R = [eye(ols1.k-1) zeros(ols1.k-1,1)];
r = zeros(ols1.k-1,1);
wt = waldsigtest(ols1,R,r);
assert(abs(wt.value - 27171.66027401781) <= tolTestStrong,'Wald joint significance')

% Reset test
rt = resettest(ols1);
assert(abs(rt.value - 17.97404570550421) <= tolTest,'RESET value')
assert(abs(rt.p - 2.67377131614e-11) <= tolTest,'RESET p')
assert(all(abs(rt.df -  [3, 808]) <= tolTest),'RESET df')

% Breush-Pagan Heteroskedasticity test
bpt = bphettest(ols1);
assert(abs(bpt.value - 80.03252229533157) <= tolTestStrong,'BPHET value')
assert(abs(bpt.p -  1.71440990325e-16) <= tolTest,'BPHET p')
assert(all(abs(bpt.df -  4) <= tolTest),'BPHET df')

wt = whitehettest(ols1);
assert(abs(wt.value - 295.2051752610454) <= tolTest,'WHIEHET value')
assert(abs(wt.p -  1.18069053179e-54) <= tolTest,'WHIEHET p')
assert(all(abs(wt.df -  14) <= tolTest),'WHIEHET df')

% OLS robust estimation
ols1r = ols(y,X,'vartype','robust');

readstderr = [0.01857351;0.01251743;0.01959449;0.00134067;0.07098896];
assert(all(abs(ols1r.stderr - readstderr) <= tolCoef),'Standard errors (Robust)')

R = [eye(ols1r.k-1) zeros(ols1r.k-1,1)];
r = zeros(ols1r.k-1,1);
wtr = waldsigtest(ols1r,R,r);
assert(abs(wtr.value - 32796.11022272207) <= tolTestLow,'Wald joint significance (Robust)')

