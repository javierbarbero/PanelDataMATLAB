% Example with Munnell (1990) data. Spatial
clear
clc

% Load data
load('..\data\MunnellData.mat')

y = log(gsp);
X = [log(pcap), log(pc), log(emp), unemp];
Z = [log(hwy), log(water)];
T = length(unique(year));

ynames = {'lgsp'};
xnames = {'lpcap', 'lpc', 'lemp', 'unemp'};

load('..\data\MunnellW.mat');

Wbig = kron(W, eye(T));

% S2SLS
regs2sls = s2sls(y,X, Wbig);
regs2sls.ynames = ynames;
regs2sls.xnames = xnames;
estdisp(regs2sls);

% Durbin S2SLS
regs2slsdur = s2sls(y,X, Wbig, 'slagX',[1 2 3 4]);
regs2slsdur.ynames = ynames;
regs2slsdur.xnames = xnames;
estdisp(regs2slsdur);

% IV S2SLS
reg2slsiv = s2sls(y, X, Wbig, 'endog', 1, 'inst', Z);
reg2slsiv.ynames = ynames;
reg2slsiv.xnames = xnames;
estdisp(reg2slsiv);

% SEM S2SLS
reg2slssem = s2sls(y, X, Wbig, 'slagy', 0, 'slagerror', 1);
reg2slssem.ynames = ynames;
reg2slssem.xnames = xnames;
estdisp(reg2slssem);

% SARAR S2SLS
reg2slssem = s2sls(y, X, Wbig, 'slagy', 1, 'slagerror', 1);
reg2slssem.ynames = ynames;
reg2slssem.xnames = xnames;
estdisp(reg2slssem);

% SARAR Endogneous S2SLS
reg2slssem = s2sls(y, X, Wbig, 'slagy', 1, 'slagerror', 1, 'endog', 1, 'inst',Z);
reg2slssem.ynames = ynames;
reg2slssem.xnames = xnames;
estdisp(reg2slssem);

% Spatial Panel FE
regsfe = spanel(id, year, y, X, W, 'fe');
regsfe.ynames = ynames;
regsfe.xnames = xnames;
estdisp(regsfe);

% Spatial Panel RE
regsre = spanel(id, year, y, X, W, 're');
regsre.ynames = ynames;
regsre.xnames = xnames;
estdisp(regsre);

% Spatial Panel EC
regsec = spanel(id, year, y, X, W, 'ec');
regsec.ynames = ynames;
regsec.xnames = xnames;
estdisp(regsec);

% Spatial Panel Durbinb FE
regsfedur = spanel(id, year, y, X, W, 'fe','slagX',[1 2 3 4]);
regsfedur.ynames = ynames;
regsfedur.xnames = xnames;
estdisp(regsfedur);

% IV Spatial Panel
regsfeiv = spanel(id, year, y, X, W, 'fe','endog',1,'inst',Z);
regsfeiv.ynames = ynames;
regsfeiv.xnames = xnames;
estdisp(regsfeiv);

% Baltagi, Song, Jung and Koh tests
bsjk = bsjksatest(regsfe);
testdisp(bsjk)


% SEM FE
regsfe = spanel(id, year, y, X, W, 'fe','slagy',0,'slagerror',1);
regsfe.ynames = ynames;
regsfe.xnames = xnames;
estdisp(regsfe);

% SEM FE Endgoneous
regsfe = spanel(id, year, y, X, W, 'fe','slagy',0,'slagerror',1,'endog',1,'inst',Z);
regsfe.ynames = ynames;
regsfe.xnames = xnames;
estdisp(regsfe);

% SARAR FE Endgoneous
regsararefe = spanel(id, year, y, X, W, 'fe','slagy',1,'slagerror',1,'endog',1,'inst',Z);
regsararefe.ynames = ynames;
regsararefe.xnames = xnames;
estdisp(regsararefe);

% SEM RE
regsfe = spanel(id, year, y, X, W, 're','slagy',0,'slagerror',1);
regsfe.ynames = ynames;
regsfe.xnames = xnames;
estdisp(regsfe);

% SEM RE Endogenous
regsfe = spanel(id, year, y, X, W, 're','slagy',0,'slagerror',1,'endog',1,'inst',Z);
regsfe.ynames = ynames;
regsfe.xnames = xnames;
estdisp(regsfe);


% SARAR RE Endgoneous
regsararere = spanel(id, year, y, X, W, 're','slagy',1,'slagerror',1,'endog',1,'inst',Z);
regsararere.ynames = ynames;
regsararere.xnames = xnames;
estdisp(regsararere);

% Spatial Hausman Tests
ht = hausmantest(regsararefe, regsararere);
testdisp(ht);

