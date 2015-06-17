* IV Models Checks

clear all
set more off

* Load Munnell Data
webuse productivity

* Generate log(public)
gen lpcap = log(public)

* IV2SLS
ivregress 2sls gsp (lpcap = hwy water) private emp unemp
matrix list e(b)
disp _se[lpcap]  _se[private] _se[emp] _se[unemp] _se[_cons]

ereturn list

test lpcap private emp unemp
return list

* Sargan test
estat overid
return list

* Hausman test
quietly reg gsp lpcap private emp unemp
estimates store ols
quietly ivregress 2sls gsp (lpcap = hwy water) private emp unemp
estimates store iv
hausman iv ols
return list

* Wu variable addition test test
estat endog
return list

* Bresuch-Pagan heteroscedasticity test
ivreg2 gsp (lpcap = hwy water) private emp unemp

ivhettest, nr2
return list

* White heteroskedasticity test
ivhettest, ivcp nr2
return list

* IV2SLS Robust
ivregress 2sls gsp (lpcap = hwy water) private emp unemp, robust
disp _se[lpcap]  _se[private] _se[emp] _se[unemp] _se[_cons]

test lpcap private emp unemp
return list

* Sargan test Robust
estat overid
return list

* Wu variable addition test test Robust
estat endog
return list
