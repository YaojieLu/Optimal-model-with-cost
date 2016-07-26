
options(digits=20)
library(deSolve)
source("Functions.r")

ca <- 400
parms <- c(LAI=1,
           Vcmax=50, cp=30, Km=703, Rd=1,
           a=1.6, nZ=0.5, p=43200, l=1.8e-5, VPD=0.02,
           pe=-1.58*10^-3, b=4.38, c=2.64, d=3.54, kxmax=5, h3=10)

k <- 0.05
MAP <- 1000
optwL1 <- optimize(wLf, c(0.18, 0.5), maximum=T)

k <- 0.1
MAP <- 1000
optwL2 <- optimize(wLf, c(0.18, 0.5), maximum=T)

k <- 0.025
MAP <- 1000
optwL4 <- optimize(wLf, c(0.18, 0.5), maximum=T)

k <- 0.05
MAP <- 2000
optwL3 <- optimize(wLf, c(0.18, 0.5), maximum=T)

k <- 0.05
MAP <- 500
optwL5 <- optimize(wLf, c(0.18, 0.5), maximum=T)
