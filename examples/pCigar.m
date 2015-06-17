% Example with Cigar data. Baltagi and Levin (1992) and Baltagi, Griffin and Xiong (2000)
clear all
clc

% Load data
load('..\data\CigarData.mat')

y = log(c);
X = [log(price) log(ndi), log(pimin)];
Z = [log(ndi_1), log(pimin_1)];

ynames = {'lc'};
xnames = {'lprice', 'lndi', 'lpimin','CONST'};
znames = {'lndi_i','lpimin_1'};

regfe = panel(state, year, y, X, 'fe');
regfe.ynames = ynames;
regfe.xnames = xnames;
estdisp(regfe);

% IV Panel FE
regivfe = ivpanel(state, year, y, X, Z, 'fe', 'endog', 1);
regivfe.ynames = ynames;
regivfe.xnames = xnames;
regivfe.znames = znames;
estdisp(regivfe);

st = sarganoitest(regivfe);
testdisp(st)

% IV Panel BE
regivbe = ivpanel(state, year, y, X, Z, 'be', 'endog', 1);
regivbe.ynames = ynames;
regivbe.xnames = xnames;
regivbe.znames = znames;
estdisp(regivbe);

% IV Panel RE
regivre = ivpanel(state, year, y, X, Z, 're', 'endog', 1);
regivre.ynames = ynames;
regivre.xnames = xnames;
regivre.znames = znames;
estdisp(regivre);

st = sarganoitest(regivre);
testdisp(st)

% IV Panel Error components
regivec = ivpanel(state, year, y, X, Z, 'ec', 'endog', 1);
regivec.ynames = ynames;
regivec.xnames = xnames;
regivec.znames = znames;
estdisp(regivec);

% Regression with mlprice (endogenous) instead of lprice to test if they are correctly
% removed in the internal 're' regression
%{
mlprice = groupmeans(state,log(price),'rep');
X1test = [mlprice];
regivretest = ivpanel(state, year, y, X1test, X2, Z, 're');
regivretest.ynames = ynames;
regivretest.xnames = {'Mlprice', 'lndi', 'lpimin','CONST'};
regivretest.znames = znames;
estdisp(regivretest);
%}

%{
% Regression with mlndi (exogenous) instead of lmndi to test if they are correctly
% removed in the internal 're' regression
mlndi = groupmeans(state,log(ndi),'rep');
X2test = [mlndi, log(pimin)];
regivretest = ivpanel(state, year, y, X1, X2test, Z, 're');
regivretest.ynames = ynames;
regivretest.xnames = {'lprice', 'Mlndi', 'lpimin','CONST'};
regivretest.znames = znames;
estdisp(regivretest);

% Regression with mlndi (exogenous) and mlndi_i (instrument) instead of lmndi to test if they are correctly
% removed in the internal 're' regression
mlndi_1 = groupmeans(state,log(ndi_1),'rep');
X2test = [mlndi, log(pimin)];
Ztest = [mlndi_1, log(pimin_1)];
regivretest = ivpanel(state, year, y, X1, X2test, Ztest, 're');
regivretest.ynames = ynames;
regivretest.xnames = {'lprice', 'Mlndi', 'lpimin','CONST'};
regivretest.znames = znames;
estdisp(regivretest);
%}



