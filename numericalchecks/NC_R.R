# Libraries
library(plm)
library(splm)

# -----------------------
# BASIC PANEL DATA MODELS
# -----------------------

# Load Data
dataMunnell = read.csv("MunnellData.csv")

fm <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp

# Panel FE
fe <- plm(fm, data=dataMunnell, model="within")
summary(fe)

# Panel BE
be <- plm(fm, data=dataMunnell, model="between")
summary(be)

# Panel RE
re <- plm(fm, data=dataMunnell, model="random")
summary(re)

# -----------------------
# INSTRUMENTAL PANELS
# -----------------------

# Load Data
dataCigar = read.csv("CigarData_1.csv")

fm <- log(c) ~ log(price) + log(ndi) + log(pimin) | log(ndi_1) + log(pimin_1) + log(ndi) + log(pimin)

# Panel IV FE
ivfe <- plm(fm, data=dataCigar, model="within")
summary(ivfe)

# Panel IV RE
ivre <- plm(fm, data=dataCigar, model="random")
summary(ivre)

# Panel IV EC
ivec <- plm(fm, data=dataCigar, model="random", inst.method="baltagi")
summary(ivec)

# -----------------------
# SPATIAL PANELS
# -----------------------

# Load Data
W = read.table("MunnellW.txt")
W = as.matrix(W)
lw = mat2listw(W)

fm <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp

# SARAR FE
sararfe <- spgm(fm, data=dataMunnell, listw=lw, lag=TRUE, spatial.error=TRUE, model="within", method="w2sls")
summary(sararfe)

# spatial coefficient (rho) and t statistic
rhofe <- sararfe$rho[1]
rhofe_t <- sararfe$rho[1] / sqrt(sararfe$rho[2]) 
print(rhofe)
print(rhofe_t)

# SARAR RE
sararre <- spgm(fm, data=dataMunnell, listw=lw, lag=TRUE, spatial.error=TRUE, model="random", method="g2sls")
summary(sararre)

# spatial coefficient (rho) and t statistic
rhore <- sararre$rho[1]
rhore_t <- sararre$rho[1] / sqrt(sararre$rho[2]) 
print(rhore)
print(rhore_t)


