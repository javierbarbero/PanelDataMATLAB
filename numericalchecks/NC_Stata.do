clear all

* -----------------------
* BASIC PANEL DATA MODELS
* -----------------------

insheet using MunnellData.csv

xtset id year

gen lgsp = log(gsp)
gen lpcap = log(pcap)
gen lpc = log(pc)
gen lemp = log(emp)

xtreg lgsp lpcap lpc lemp unemp, fe

xtreg lgsp lpcap lpc lemp unemp, be

xtreg lgsp lpcap lpc lemp unemp, re

* -----------------------
* INSTRUMENTAL PANELS
* -----------------------

clear all

insheet using CigarData.csv

xtset state year

gen c_1 = l.c
gen price_1 = l.price
gen ndi_1 = l.ndi
gen pimin_1 = l.pimin

gen lc = log(c)
gen lc_1 = log(c_1)
gen lprice = log(price)
gen lprice_1 = log(price_1)
gen lndi = log(ndi)
gen lndi_1 = log(ndi_1)
gen lpimin = log(pimin)
gen lpimin_1 = log(pimin_1)

drop if lc_1 >= .

xtivreg lc lndi lpimin (lprice = lndi_1 lpimin_1), fe

xtivreg lc lndi lpimin (lprice = lndi_1 lpimin_1), re

xtivreg lc lndi lpimin (lprice = lndi_1 lpimin_1), ec


