
options(digits=20)
library(deSolve)
source("Functions.r")

ca <- 400
c <- 2.64
d <- 3.54

parms <- c(LAI=1,
           Vcmax=50, cp=30, Km=703, Rd=1,
           a=1.6, nZ=0.5, p=43200, l=1.8e-5, VPD=0.02,
           pe=-1.58*10^-3, b=4.38, kxmax=5)

k <- 0.05
MAP <- 1000
h3 <- 10

optwL <- optimize(wLf, c(0.18, 0.5), maximum=T)
