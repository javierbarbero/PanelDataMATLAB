clear all
clc

% Tolerance
tolCoef = 1e-6;
tolCoefLow = 1e-2;
tolTest = 1e-2;


% Read data
load '..\data\MunnellData'
y = log(gsp);
X = [log(pcap), log(pc), log(emp), unemp];
Z = [log(hwy), log(water)];

load '..\data\MunnellW'
Wbig = kron(W,eye(17));

% POOL unit test of Pool SAR are with respect to Sata results (std errors differe
% betwee Stata ad R)
% POOL SAR
poSAR = s2sls(y,X,Wbig);

assert(all(abs(poSAR.coef - [0.1474823; 0.3092149; 0.6026597; -0.0061726; -0.0092512; 1.7486408]) <= tolCoef),'Estimated coefficients')
assert(all(abs(poSAR.stderr - [0.0178043; 0.0102487; 0.0148493; 0.0014596  ;  0.0060324; 0.0895507]) <= tolCoefLow),'Standard errors')

% POOL SAR Endog
poSARen = s2sls(y,X,Wbig,'endog',1,'inst',Z);

assert(all(abs(poSARen.coef - [0.1673200; 0.3039439; 0.5886492; -0.0064108;  -0.0065639;  1.6833522]) <= tolCoef),'Estimated coefficients')
assert(all(abs(poSARen.stderr - [0.0219589; 0.0108646; 0.0173264; 0.0014708; 0.0061921; 0.0975046]) <= tolCoefLow),'Standard errors')

% POOL SEM
poSEM = s2sls(y,X,Wbig,'slagy',0,'slagerror',1);

assert(all(abs(poSEM.coef - [0.1444635; 0.3545111; 0.5683849; -0.0080987; 1.456315 ]) <= tolCoefLow),'Estimated coefficients')
assert(all(abs(poSEM.stderr - [0.0164729 ; 0.0109119; 0.0143316;  0.001697; 0.0578171]) <= tolCoefLow),'Standard errors')

% POOL SEM Endog (Differences from R and Stata)
poSEMen = s2sls(y,X,Wbig,'slagy',0,'slagerror',1,'endog',1,'inst',Z);

%assert(all(abs(poSEMen.coef - [0.2068948; 0.3362050; 0.5262000; -0.0086131; 1.3430265]) <= tolCoefLow),'Estimated coefficients')
%assert(all(abs(poSEMen.stderr - [0.0218915 ; 0.0117472; 0.0174271;  0.001697; 0.0017161]) <= tolCoefLow),'Standard errors')

% ----------------
% FIXED EFFECTS
% ----------------

% FE SAR
feSAR = spanel(id,year,y,X,W,'fe');

assert(all(abs(feSAR.coef - [-0.04040614; 0.21904067; 0.66833361; -0.00472828; 0.19166263]) <= tolCoef),'Estimated coefficients')
assert(all(abs(feSAR.stderr - [0.02666502; 0.02509769; 0.03077822; 0.00090997; 0.02617774]) <= tolCoef),'Standard errors')

% FE SAR Endog
feSARen = spanel(id,year,y,X,W,'fe','endog',1,'inst',Z);

assert(all(abs(feSARen.coef - [0.02112857; 0.21205733;  0.64455161; -0.00557698 ;  0.18486869]) <= tolCoef),'Estimated coefficients')
assert(all(abs(feSARen.stderr - [0.03465079; 0.02529328; 0.03183561; 0.00096473; 0.02664021]) <= tolCoef),'Standard errors')

% FE SEM
feSEM = spanel(id,year,y,X,W,'fe','slagy',0,'slagerror',1);

assert(all(abs(feSEM.coef - [0.0043026; 0.2144604; 0.7830897; -0.0025609]) <= tolCoef),'Estimated coefficients')
assert(all(abs(feSEM.stderr - [0.0253425; 0.0232533; 0.0279794; 0.0010546]) <= tolCoef),'Standard errors')

% FE SEM Endog
feSEMen = spanel(id,year,y,X,W,'fe','endog',1,'inst',Z,'slagy',0,'slagerror',1);

assert(all(abs(feSEMen.coef - [0.0738375; 0.2048004;  0.7522485; -0.0034398]) <= tolCoef),'Estimated coefficients')
assert(all(abs(feSEMen.stderr - [0.0364745; 0.0246330; 0.0314305; 0.0011619]) <= tolCoef),'Standard errors')

% FE SARAR
feSARAR = spanel(id,year,y,X,W,'fe','slagerror',1);

assert(all(abs(feSARAR.coef - [-0.0205827; 0.1936870; 0.7291745; -0.0037004; 0.1327087]) <= tolCoef),'Estimated coefficients')
assert(all(abs(feSARAR.stderr - [0.0268688; 0.0255383; 0.0303750; 0.0010235; 0.0245926]) <= tolCoef),'Standard errors')

% FE SARAR Endog
feSARARen = spanel(id,year,y,X,W,'fe','endog',1,'inst',Z,'slagerror',1);

assert(all(abs(feSARARen.coef - [0.0264323; 0.1885953; 0.7131347; -0.0042630; 0.1244800]) <= tolCoef),'Estimated coefficients')
assert(all(abs(feSARARen.stderr - [0.0352005; 0.0256520; 0.0315721; 0.0010737; 0.0249189]) <= tolCoef),'Standard errors')

% ----------------
% RANDOM EFFECTS
% ----------------

% RE SAR
reSAR = spanel(id,year,y,X,W,'re');

assert(all(abs(reSAR.coef - [0.02098103; 0.29001525; 0.71011141; -0.00641009; 0.03974083; 1.91197495]) <= tolCoef),'Estimated coefficients')
assert(all(abs(reSAR.stderr - [0.02475292; 0.02117623; 0.02677249; 0.00090826; 0.01508851; 0.16545521]) <= tolCoef),'Standard errors')

% RE SAR Endog
reSARen = spanel(id,year,y,X,W,'re','endog',1,'inst',Z);

assert(all(abs(reSARen.coef - [0.09117218; 0.27813078; 0.66319341; -0.00734469; 0.05160749; 1.56678660]) <= tolCoef),'Estimated coefficients')
assert(all(abs(reSARen.stderr - [0.03036366; 0.02137782; 0.02920288; 0.00093865;  0.01531029; 0.18610899]) <= tolCoef),'Standard errors')

% RE SEM
reSEM = spanel(id,year,y,X,W,'re','slagy',0,'slagerror',1);

assert(all(abs(reSEM.coef - [0.0533878; 0.2587524; 0.7268627; -0.0039258; 2.2178061]) <= tolCoef),'Estimated coefficients')
assert(all(abs(reSEM.stderr - [0.0221395; 0.0210013; 0.0253709; 0.0011000;  0.1352650]) <= tolCoef),'Standard errors')

% RE SEM Endog
reSEMen = spanel(id,year,y,X,W,'re','endog',1,'inst',Z,'slagy',0,'slagerror',1);

assert(all(abs(reSEMen.coef - [0.0890410; 0.2532223; 0.7036796; -0.0044968; 2.0962931]) <= tolCoef),'Estimated coefficients')
assert(all(abs(reSEMen.stderr - [0.0332734; 0.0209522; 0.0302675; 0.0011174; 0.1607192]) <= tolCoef),'Standard errors')

% RE SARAR
reSARAR = spanel(id,year,y,X,W,'re','slagerror',1);

assert(all(abs(reSARAR.coef - [0.04632588; 0.26797169; 0.72014854; -0.00523286; 0.02230657; 2.00687954]) <= tolCoef),'Estimated coefficients')
assert(all(abs(reSARAR.stderr - [0.02268646; 0.02047296; 0.02493860; 0.00097817; 0.01354214; 0.16835095]) <= tolCoef),'Standard errors')

% RE SARAR Endog
reSARARen = spanel(id,year,y,X,W,'re','endog',1,'inst',Z,'slagerror',1);

assert(all(abs(reSARARen.coef - [0.1066485; 0.2567983; 0.6806213; -0.0059865; 0.0311608; 1.7287076]) <= tolCoef),'Estimated coefficients')
assert(all(abs(reSARARen.stderr - [0.0276903; 0.0207423; 0.0271849; 0.0010105; 0.0138498; 0.1858177]) <= tolCoef),'Standard errors')

% ----------------
% ERROR COMPONENTS 
% ----------------

% EC SAR
ecSAR = spanel(id,year,y,X,W,'ec');

assert(all(abs(ecSAR.coef - [0.02244269; 0.28871841; 0.70835224; -0.00643463; 0.04256929; 1.89419519]) <= tolCoef),'Estimated coefficients')
assert(all(abs(ecSAR.stderr - [0.02472283; 0.02114987; 0.02673746; 0.00090743; 0.01501675; 0.16510513]) <= tolCoef),'Standard errors')

% RE SAR Endog
ecSARen = spanel(id,year,y,X,W,'ec','endog',1,'inst',Z);

assert(all(abs(ecSARen.coef - [0.09292665; 0.27706012; 0.66144585; -0.00737182; 0.05380668; 1.55034237]) <= tolCoef),'Estimated coefficients')
assert(all(abs(ecSARen.stderr - [0.02987237; 0.02134726; 0.02897316; 0.00093558; 0.01524689; 0.18417731]) <= tolCoef),'Standard errors')


% ----------------
% TESTS 
% ----------------

% BSJK test
bsjk = bsjksatest(reSARAR);

assert(abs(bsjk.value - 4290.422) <= tolTest,'BSJK value')
assert(abs(bsjk.p -  0) <= tolTest,'BSJK p')
assert(all(abs(bsjk.df -  3) <= tolTest),'BSJK df')


% Spatial hausman
haus = hausmantest(feSARARen, reSARARen);

assert(abs(haus.value - 33.13798) <= tolTest,'HAUSMAN value')
assert(abs(haus.p -  3.5333e-06) <= tolTest,'HAUSMAN p')
assert(all(abs(haus.df -  5) <= tolTest),'HAUSMAN df')



