* Linear Models Checks

clear all
set more off

* Load Munnell Data
webuse productivity

* Generate log(public)
gen lpcap = log(public)

* OLS
reg gsp lpcap private emp unemp
matrix list e(b)
disp _se[lpcap]  _se[private] _se[emp] _se[unemp] _se[_cons]

ereturn list

test lpcap private emp unemp
return list

* Ramsey RESET
estat ovtest 
return list

* Bresuch-Pagan heteroscedasticity test
estat hettest, rhs iid
return list

* White heteroskedasticity test
estat imtest, white
return list

* OLS Robust
reg gsp lpcap private emp unemp, robust
disp _se[lpcap]  _se[private] _se[emp] _se[unemp] _se[_cons]

test lpcap private emp unemp
return list
