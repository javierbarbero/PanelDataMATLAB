# Libraries
library(splm)
library(sphet)
library(spdep)

# Load Data
data = read.csv("..\\data\\MunnellData.csv")
W = read.table("..\\data\\MunnellW.txt")
W = as.matrix(W)
lw = mat2listw(W)

fm <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
fme <- log(gsp) ~ log(pc) + log(emp) + unemp

# Build big W por Pool
bigWm <- kronecker(W,diag(17))
bigW <- mat2listw(bigWm)

# Export to geoda .GWT file (to import later into Stata)
write.sn2gwt(listw2sn(mat2listw(bigWm)), '..\\data\\bigW.GWT')

# Pool SAR
poSAR <- spreg(fm, data=data, listw=bigW, model="lag")
summary(poSAR)

# Pool SAR Endog
poSARen <- spreg(fme, data=data, listw=bigW, model="lag", endog = ~log(pcap), instruments = ~log(hwy) + log(water))
summary(poSARen)

# Pool SEM
poSEM <- spreg(fm, data=data, listw=bigW, model="error")
summary(poSEM)

# Pool SEM Endog
poSEMen <- spreg(fme, data=data, listw=bigW, model="error", endog = ~log(pcap), instruments = ~log(hwy) + log(water))
summary(poSEMen)

# ---------------
# FIXED EFFECTS
# ---------------

# FE SAR
feSAR <- spgm(fm, data=data, listw=lw, lag=TRUE, spatial.error=FALSE, model="within", method="w2sls")
summary(feSAR)

# FE SAR Endog
feSARen <- spgm(fme, data=data, listw=lw, lag=TRUE, spatial.error=FALSE, model="within", method="w2sls", endog= ~log(pcap), instruments= ~log(hwy) + log(water))
summary(feSARen)

# FE SEM
feSEM <- spgm(fm, data=data, listw=lw, lag=FALSE, spatial.error=TRUE, model="within", method="w2sls")
summary(feSEM)

# FE SEM Endog
feSEMen <- spgm(fme, data=data, listw=lw, lag=FALSE, spatial.error=TRUE, model="within", method="w2sls", endog= ~log(pcap), instruments= ~log(hwy) + log(water))
summary(feSEMen)

# FE SARAR
feSARAR <- spgm(fm, data=data, listw=lw, lag=TRUE, spatial.error=TRUE, model="within", method="w2sls")
summary(feSARAR)

# FE SARAR Endog
feSARARen <- spgm(fme, data=data, listw=lw, lag=TRUE, spatial.error=TRUE, model="within", method="w2sls", endog= ~log(pcap), instruments= ~log(hwy) + log(water))
summary(feSARAR)

# ---------------
# RANDOM EFFECTS
# ---------------

# RE SAR
reSAR <- spgm(fm, data=data, listw=lw, lag=TRUE, spatial.error=FALSE, model="random", method="g2sls")
summary(reSAR)

# RE SAR Endog
reSARen <- spgm(fme, data=data, listw=lw, lag=TRUE, spatial.error=FALSE, model="random", method="g2sls", endog= ~log(pcap), instruments= ~log(hwy) + log(water))
summary(reSARen)

# RE SEM
reSEM <- spgm(fm, data=data, listw=lw, lag=FALSE, spatial.error=TRUE, model="random", method="g2sls")
summary(reSEM)

# FE SEM Endog
reSEMen <- spgm(fme, data=data, listw=lw, lag=FALSE, spatial.error=TRUE, model="random", method="g2sls", endog= ~log(pcap), instruments= ~log(hwy) + log(water))
summary(reSEMen)

# RE SARAR
reSARAR <- spgm(fm, data=data, listw=lw, lag=TRUE, spatial.error=TRUE, model="random", method="g2sls")
summary(reSARAR)

# RE SARAR Endog
reSARARen <- spgm(fme, data=data, listw=lw, lag=TRUE, spatial.error=TRUE, model="random", method="g2sls", endog= ~log(pcap), instruments= ~log(hwy) + log(water))
summary(reSARAR)

# ---------------
# ERROR COMPONENTS EFFECTS
# ---------------

# EC SAR
ecSAR <- spgm(fm, data=data, listw=lw, lag=TRUE, spatial.error=FALSE, model="random", method="ec2sls")
summary(ecSAR)

# EC SAR Endog
ecSARen <- spgm(fme, data=data, listw=lw, lag=TRUE, spatial.error=FALSE, model="random", method="ec2sls", endog= ~log(pcap), instruments= ~log(hwy) + log(water))
summary(ecSARen)

# ---------------
# TESTS EFFECTS
# ---------------

# BSJK
bsjktest(fm, data=data, listw=lw, test="J")

# Spatial Hausman Tests
haus <- sphtest(reSARARen, feSARARen)
