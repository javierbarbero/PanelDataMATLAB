clear all
clc

% Example with Unbalanced Munnell (1990) data
disp('-----------------------------------------');
disp('            PANEL DATA MODELS            ');
disp('-----------------------------------------');

% Load data
% Exported from Stata after:
% drop if public > 135000
load('MunnellUnbalanced');

y = gsp;
X = [lpcap, private, emp, unemp];
id = state;

% Sort data
[~, idx] = sortrows(X,[2]);
id = id(idx);
year = year(idx);
y = y(idx,:);
X = X(idx,:);

ynames = {'lgsp'};
xnames = {'lpcap', 'lpc', 'lemp', 'unemp','CONST'};

% OLS
regols = ols(y,X);
regols.ynames = ynames;
regols.xnames = xnames;
estdisp(regols);

% Panel FE
regfe = panel(id,year,y, X, 'fe');
regfe.ynames = ynames;
regfe.xnames = xnames;
estdisp(regfe);

% Panel BE
regbe = panel(id,year,y, X, 'be');
regbe.ynames = ynames;
regbe.xnames = xnames;
estdisp(regbe);

% Panel RE
regre = panel(id,year,y, X, 're');
regre.ynames = ynames;
regre.xnames = xnames;
estdisp(regre);

% F test of inividual effects
effF = effectsftest(regfe);
testdisp(effF);

% Hausman test
h = hausmantest(regfe,regre);
testdisp(h);

% BP test for effects
bpre = bpretest(regre);
testdisp(bpre);

% Mundlak test
mu = mundlakvatest(regfe);
testdisp(mu);

% Pool test
% po = pooltest(regfe);
% testdisp(po);

% Wooldridge serial test
wo = woolserialtest(regfe);
testdisp(wo);

wo = woolserialtest(regfe,'dfcorrection',0);
testdisp(wo);

% Pesaran CSD
pecsdfe = pesarancsdtest(regfe);
testdisp(pecsdfe);

pecsdre = pesarancsdtest(regre);
testdisp(pecsdre);

% Panel FE Robust
regfer = panel(id,year,y, X, 'fe', 'vartype','robust');
regfer.ynames = ynames;
regfer.xnames = xnames;
estdisp(regfer);

% Panel RE Robust
regrer = panel(id,year,y, X, 're', 'vartype','robust');
regrer.ynames = ynames;
regrer.xnames = xnames;
estdisp(regrer);

% Mundlak test
mur = mundlakvatest(regfer);
testdisp(mur);

% Individual effects
ieffectsdisp(regfer);

% Individual effects
ieffectsdisp(regfer,'overall');


