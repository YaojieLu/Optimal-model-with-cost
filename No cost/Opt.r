
######################### No cost #########################

options(digits=20)
library(deSolve)
source("No cost/Functions.r")

ca <- 400
k <- 0.05
MAP <- 3650
parms <- c(LAI=1,
           Vcmax=50, cp=30, Km=703, Rd=1,
           a=1.6, nZ=0.5, p=43200, l=1.8e-5, VPD=0.02)

muopt <- uniroot(muf, c(-15, -10), tol=.Machine$double.eps)
muopt
f <- Vectorize(function(w)gswf(w, muopt$root))
f(1)