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

% Pool
reg = ols(y,X);
regpo = panel(id,year,y,X,'po');
assert(all(abs(regpo.coef - reg.coef) <= tolCoef),'Estimated coefficients')
assert(all(abs(regpo.stderr - reg.stderr) <= tolCoef),'Standard errors')
assert(abs(regpo.r2 - reg.r2) <= tolCoef,'R2')
assert(abs(regpo.adjr2 - reg.adjr2) <= tolCoef,'Adj R2')
assert(abs(regpo.RSS -  reg.RSS) <= tolCoef,'RSS')
assert(abs(regpo.resdf - reg.resdf) <= tolCoef,'res DF')

% FE
regfe = panel(id,year,y,X,'fe');

realcoeff = [0-.0261493; 0.29200668; 0.7681595;  -0.00529775];
assert(all(abs(regfe.coef - realcoeff) <= tolCoef),'Estimated coefficients')

readstderr = [0.02900156; 0.02511967; 0.03009174; 0.00098873];
assert(all(abs(regfe.stderr - readstderr) <= tolCoef),'Standard errors')

assert(abs(regfe.r2 - 0.9413356320514598) <= tolCoef,'R2')
assert(abs(regfe.adjr2 - 0.9374195551334289) <= tolCoef,'Adj R2')
assert(abs(regfe.RSS -  1.111188158537215) <= tolCoef,'RSS')
assert(abs(regfe.resdf - 764) <= tolCoef,'res DF')

R = [eye(regfe.k)];
r = zeros(regfe.k,1);
wt = waldsigtest(regfe,R,r);
assert(abs(wt.value - 3064.809389569204) <= tolTestStrong,'Wald joint significance')

% Overll individual effect
[oeff, soeff] = ieffects(regfe,'overall');
assert(all(abs([oeff, soeff] - [ 2.3528981, 0.17481312 ]) <= tolTestLow),'Wald joint significance (Robust)')

% Test of individual effects
effF = effectsftest(regfe);

assert(abs(effF.value - 75.82042805979164) <= tolTest,'EFFECTSF value')
assert(abs(effF.p -  0.0000) <= tolTest,'EFFECTSF p')
assert(all(abs(effF.df -  [47, 764]) <= tolTest),'EFFECTSF df')

% BE
regbe = panel(id,year,y,X,'be');

realcoeff = [0.17936513; 0.30195423; 0.57612738;  -0.00389029; 1.5894443];
assert(all(abs(regbe.coef - realcoeff) <= tolCoef),'Estimated coefficients')

readstderr = [0.07197194; 0.04182148; 0.05637459; 0.00990835; 0.23297958];
assert(all(abs(regbe.stderr - readstderr) <= tolCoef),'Standard errors')

assert(abs(regbe.r2 - 0.9939088796508508) <= tolCoef,'R2')
assert(abs(regbe.adjr2 - 0.9933422638044184) <= tolCoef,'Adj R2')
assert(abs(regbe.RSS -  0.2977007950047585) <= tolCoef,'RSS')
assert(abs(regbe.resdf - 43) <= tolCoef,'res DF')

R = [eye(regbe.k-1) zeros(regbe.k-1,1)];
r = zeros(regbe.k-1,1);
wt = waldsigtest(regbe,R,r);
assert(abs(wt.value - 1754.114160252832) <= tolTestStrong,'Wald joint significance')

% RE
regre = panel(id,year,y,X,'re');

realcoeff = [0.0044388; 0.31054829; 0.72967051; -0.00617248; 2.1354106];
assert(all(abs(regre.coef - realcoeff) <= tolCoef),'Estimated coefficients')

readstderr = [0.02341731; 0.01980475; 0.02492022; 0.00090728; 0.13346148];
assert(all(abs(regre.stderr - readstderr) <= tolCoef),'Standard errors')

assert(abs(regre.r2 - 0.9916667034632465) <= tolCoef,'R2')
%assert(abs(regre.adjr2 - 0.9933422638044184) <= tolCoef,'Adj R2')
%assert(abs(regre.RSS -  0.2977007950047585) <= tolCoef,'RSS')
%assert(abs(regre.resdf - 43) <= tolCoef,'res DF')

R = [eye(regre.k-1) zeros(regre.k-1,1)];
r = zeros(regre.k-1,1);
wt = waldsigtest(regre,R,r);
assert(abs(wt.value - 19131.08908355508) <= tolTestStrong,'Wald joint significance')

% Hausman test
ht = hausmantest(regfe, regre);

assert(abs(ht.value - 9.525397763992501) <= tolTest,'HAUSMAN value')
assert(abs(ht.p -  0.0492279877199629) <= tolTest,'HAUSMAN p')
assert(all(abs(ht.df -  4) <= tolTest),'HAUSMAN df')

% Mundlak test
mu = mundlakvatest(regfe);

assert(abs(mu.value - 9.718078857971815) <= tolTest,'MUNDLAK value')
assert(abs(mu.p -  0.0454540315736643) <= tolTest,'MUNDLAK p')
assert(all(abs(mu.df -  4) <= tolTest),'MUNDLAK df')

% BP test for effects
bpef = bpretest(regre);

assert(abs(bpef.value -  4134.96140397660) <= tolTestLow,'BPEF value')
assert(abs(bpef.p -  0) <= tolTest,'BPEF p')
assert(all(abs(bpef.df -  1) <= tolTest),'BPEF df')

% Wooldridge serial correlation test
wos = woolserialtest(regfe);

assert(abs(wos.value -   680.2993384282088) <= tolTestLow,'WOOLSERIAL value')
assert(abs(wos.p -  1.32286995196e-29) <= tolTest,'WOOLSERIAL p')
assert(all(abs(wos.df -  [1, 47]) <= tolTest),'WOOLSERIAL df')

% Baltagi Li (1991) (1995 book) serial correaltion and random effects
blse = blserialtest(regre);

assert(abs(blse.value -   4187.597268866712) <= tolTestLow,'BLSERIAL value')
assert(abs(blse.p -  0) <= tolTest,'BLSERIAL p')
assert(all(abs(blse.df -  2) <= tolTest),'BLSERIAL df')


% Pesaran CSD test
pecsdfe = pesarancsdtest(regfe);
assert(abs(pecsdfe.value -   30.36845266767461) <= tolTest,'PESARANCSD value')

pecsdre = pesarancsdtest(regre);
assert(abs(pecsdre.value -   29.07877056772282) <= tolTest,'PESARANCSD value')

% -----------------
% ROBUST
% -----------------

% FE Robust
regfer = panel(id,year,y,X,'fe','vartype','robust');

readstderr = [0.06111477; 0.06254953; 0.08273258; 0.00252846];
assert(all(abs(regfer.stderr - readstderr) <= tolCoef),'Standard errors')

R = [eye(regfer.k)];
r = zeros(regfer.k,1);
wt = waldsigtest(regfer,R,r);
assert(abs(wt.value - 395.6090100801947) <= tolTestStrong,'Wald joint significance')

% Overll individual effect Robust
[oeff, soeff] = ieffects(regfer,'overall');
assert(all(abs([oeff, soeff] - [ 2.3528981, 0.31459401 ]) <= tolTestLow),'Wald joint significance (Robust)')

% RE Robust
regrer = panel(id,year,y,X,'re','vartype','robust');

readstderr = [0.05531069; 0.04416204; 0.07088248; 0.00236311; 0.24178717];
assert(all(abs(regrer.stderr - readstderr) <= tolCoef),'Standard errors')

R = [eye(regrer.k-1) zeros(regrer.k-1,1)];
r = zeros(regrer.k-1,1);
wt = waldsigtest(regrer,R,r);
assert(abs(wt.value - 4408.636408888779) <= tolTestStrong,'Wald joint significance')

% Mundlak test Robust
mur = mundlakvatest(regfer);

assert(abs(mur.value - 19.33303861494128) <= tolTest,'MUNDLAK value')
assert(abs(mur.p -  .0006759377929905) <= tolTest,'MUNDLAK p')
assert(all(abs(mur.df -  4) <= tolTest),'MUNDLAK df')




