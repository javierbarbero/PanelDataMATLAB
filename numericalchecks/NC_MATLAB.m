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
xnames = {'lpcap','lpc','lemp','unemp'};

% Fixed effects
fe = panel(id, year, y, X, 'fe');
fe.ynames = ynames;
fe.xnames = xnames;
estdisp(fe);

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

% -----------------------
% INSTRUMENTAL PANELS
% -----------------------

% Load data
load('CigarData')

y = log(c);
X = [log(price), log(ndi), log(pimin)];
Z = [log(ndi_1), log(pimin_1)];

ynames = {'lc'};
xnames = {'lprice','lndi','lpimin'};
znames = {'lndi_1','lpimin_1'};

% IV Fixed effects
ivfe = ivpanel(state, year, y, X, Z, 'fe', 'endog', 1);
ivfe.ynames = ynames;
ivfe.xnames = xnames;
ivfe.znames = znames;
estdisp(ivfe);

% IV Random
ivre = ivpanel(state, year, y, X, Z, 're', 'endog', 1);
ivre.ynames = ynames;
ivre.xnames = xnames;
ivre.znames = znames;
estdisp(ivre);

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
xnames = {'lpcap','lpc','lemp','unemp'};

% SARAR Fe
sararfe = spanel(id, year, y, X, W, 'fe', 'slagy', 1, 'slagerror', 1);
sararfe.ynames = ynames;
sararfe.xnames = xnames;
estdisp(sararfe);

% SARAR Re
sararre = spanel(id, year, y, X, W, 're', 'slagy', 1, 'slagerror', 1);
sararre.ynames = ynames;
sararre.xnames = xnames;
estdisp(sararre);

