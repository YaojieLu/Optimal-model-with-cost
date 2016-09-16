
options(digits=20)
library(deSolve)
source("Functions 1.r")
source("Functions 2.r")

ca <- 400

parms <- c(LAI=1,
           Vcmax=50, cp=30, Km=703, Rd=1,
           a=1.6, nZ=0.5, p=43200, l=1.8e-5, VPD=0.02,
           pe=-1.58*10^-3, b=4.38, kxmax=5, c=2.64, d=3.54)

k <- 0.05
MAP <- 1000
h3 <- 10

#averBf(0.187957156405029, parms)#0.187957156405029
optwL <- optimize(wLf, c(0.18, 0.21), maximum=T)
