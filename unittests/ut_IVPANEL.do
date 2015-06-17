* IV Models Checks

clear all
set more off

* Load Munnell Data
insheet using "..\data\CigarData_1.csv", comma

* Generate log of variables
gen lc = log(c)
gen lprice = log(price)
gen lndi = log(ndi)
gen lpimin = log(pimin)
gen lndi_1 = log(ndi_1)
gen lpimin_1 = log(pimin_1)

* Set the panel
xtset state year

* Fixed effects
xtivreg lc (lprice = lndi_1 lpimin_1) lndi lpimin, fe
matrix list e(b)
disp _se[lprice]  _se[lndi] _se[lpimin]  _se[_cons]

ereturn list

* Sargan test
xtivreg2 lc (lprice = lndi_1 lpimin_1) lndi lpimin, fe
ereturn list

* Between estimation
xtivreg lc (lprice = lndi_1 lpimin_1) lndi lpimin, be
matrix list e(b)
disp _se[lprice]  _se[lndi] _se[lpimin]  _se[_cons]

ereturn list

* Random effects
xtivreg lc (lprice = lndi_1 lpimin_1) lndi lpimin, re
matrix list e(b)
disp _se[lprice]  _se[lndi] _se[lpimin]  _se[_cons]

ereturn list

* Error components
xtivreg lc (lprice = lndi_1 lpimin_1) lndi lpimin, ec
matrix list e(b)
disp _se[lprice]  _se[lndi] _se[lpimin]  _se[_cons]

ereturn list
