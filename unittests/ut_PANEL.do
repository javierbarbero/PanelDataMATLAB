* Linear Models Checks

clear all
set more off

* Load Munnell Data
webuse productivity

* Generate log(public)
gen lpcap = log(public)

* Set the panel
xtset state year

* FE
xtreg gsp lpcap private emp unemp, fe
estimates store fixed
matrix list e(b)
disp _se[lpcap]  _se[private] _se[emp] _se[unemp] _se[_cons]

ereturn list

* BE
xtreg gsp lpcap private emp unemp, be
matrix list e(b)
disp _se[lpcap]  _se[private] _se[emp] _se[unemp] _se[_cons]

ereturn list

* RE
xtreg gsp lpcap private emp unemp, re
estimates store random
matrix list e(b)
disp _se[lpcap]  _se[private] _se[emp] _se[unemp] _se[_cons]

ereturn list

* BP test
xttest0
return list

* Hausman test
hausman fixed random 
return list

* Mundlak test
sort state
by state: egen mlpcap = mean(lpcap)
by state: egen mprivate = mean(private)
by state: egen memp = mean(emp)
by state: egen munemp = mean(unemp)

xtreg gsp lpcap private emp unemp mlpcap mprivate memp munemp, re

test mlpcap mprivate memp munemp
return list

* Woolridge serial correlation test
xtreg gsp lpcap private emp unemp, fe
predict res, e

sort state year
reg res l.res, cluster(state)
test l.res = -0.062500 
return list

* Baltagi Li serial correaltion test
xtreg gsp lpcap private emp unemp, re
xttest1, unadjust
return list

* Pesaran CSD test
xtreg gsp lpcap private emp unemp, fe
xtcsd, pesaran
return list

xtreg gsp lpcap private emp unemp, re
xtcsd, pesaran
return list

* ----------
* ROBUST
* ----------

* FE Robust
xtreg gsp lpcap private emp unemp, fe robust
matrix list e(b)
disp _se[lpcap]  _se[private] _se[emp] _se[unemp] _se[_cons]

ereturn list

* RE Robust
xtreg gsp lpcap private emp unemp, re robust
matrix list e(b)
disp _se[lpcap]  _se[private] _se[emp] _se[unemp] _se[_cons]

ereturn list

* Mundlak test Robust
xtreg gsp lpcap private emp unemp mlpcap mprivate memp munemp, re robust

test mlpcap mprivate memp munemp
return list
