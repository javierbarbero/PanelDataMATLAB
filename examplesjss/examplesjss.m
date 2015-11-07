clear
clc

% -----------------------
% BASIC PANEL DATA MODELS
% -----------------------

% Load data
load('MunnellData')

y = log(gsp);
X = [log(pcap), log(pc), log(emp), unemp];

ynames = {'lgsp'};
xnames = {'lpcap', 'lpc', 'lemp', 'unemp'};

% Fixed effects
fe = panel(id, year, y, X, 'fe');
fe.ynames = ynames;
fe.xnames = xnames;
estdisp(fe);

% Individual effects
ieff = ieffects(fe);
ieffectsdisp(fe);

% 'overall' constant term
ieffOver = ieffects(fe, 'overall');
ieffectsdisp(fe, 'overall');

% Between estimation
be = panel(id, year, y, X, 'be');
be.ynames = ynames;
be.xnames = xnames;
estdisp(be);

% Random effects
re = panel(id, year, y, X, 're');
re.ynames = ynames;
re.xnames = xnames;
estdisp(re);

% Confidence interavls
estcidisp(re);

% Robust standard errors
fer = panel(id, year, y, X, 'fe', 'vartype', 'robust');
fer.ynames = ynames;
fer.xnames = xnames;
estdisp(fer);

% -----------------------
% INSTRUMENTAL PANELS
% -----------------------

% Load data
load('CigarData')

y = log(c);
X = [log(price), log(ndi), log(pimin)];
Z = [log(ndi_1), log(pimin_1)];

ynames = {'lc'};
xnames = {'lprice', 'lndi', 'lpimin'};
znames = {'lndi_1', 'lpimin_1'};

% IV Fixed effects
ivfe = ivpanel(state, year, y, X, Z, 'fe', 'endog', 1);
ivfe.ynames = ynames;
ivfe.xnames = xnames;
ivfe.znames = znames;
estdisp(ivfe);

% EC2SLS
ec2sls = ivpanel(state, year, y, X, Z, 'ec', 'endog', 1);
ec2sls.ynames = ynames;
ec2sls.xnames = xnames;
ec2sls.znames = znames;
estdisp(ec2sls);

% -----------------------
% SPATIAL PANELS
% -----------------------

% Load data
load('MunnellData')
load('MunnellW')

y = log(gsp);
X = [log(pcap), log(pc), log(emp), unemp];

ynames = {'lgsp'};
xnames = {'lpcap', 'lpc', 'lemp', 'unemp'};

% SARAR
sarar = spanel(id, year, y, X, W, 're', 'slagy', 1, 'slagerror', 1);
sarar.ynames = ynames;
sarar.xnames = xnames;
estdisp(sarar);

% SARAR with endogenous
Z = [log(hwy), log(water)];

sarfe = spanel(id, year, y, X, W, 'fe', 'slagy', 1, 'slagerror', 1,...
                'endog', 1, 'inst', Z);
sarfe.ynames = ynames;
sarfe.xnames = xnames;
estdisp(sarfe);

% -----------------------
% TESTS
% -----------------------

% Testing Liner Hypothesis
% - Wald joint significance test
R = [1 0 0 0 0; 0 1 0 0 0];
r = [0; 0];

wald = waldsigtest(re, R, r);
testdisp(wald);

% Testing Poolability
pool = pooltest(re);
testdisp(pool);

% Testing individual effects
% - Chow F test
effF = effectsftest(fe);
testdisp(effF)

% - BP random effects
bpre = bpretest(re);
testdisp(bpre);

% Texting fixed vs. random effects
% - Hausman test
hausman = hausmantest(fe, re);
testdisp(hausman);

% - Mundlak
mundlak = mundlakvatest(fe);
testdisp(mundlak);

% Testing serial correlation
% - Wooldridge FE
woolfe = woolserialtest(fe);
testdisp(woolfe);

% - Baltagi Li RE
blre = blserialtest(re);
testdisp(blre);

% Testing cross-sectional dependence
% - Pesaran
pesaran = pesarancsdtest(fe);
testdisp(pesaran);

% Testing overidentification
% - Sargan
sargan = sarganoitest(ivfe);
testdisp(sargan);

% Testing spatial autocorrelation
% - BSJK
bsjk = bsjksatest(sarar);
testdisp(bsjk);



