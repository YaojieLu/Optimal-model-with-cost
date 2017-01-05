######################### No cost #########################

options(digits=20)
library(deSolve)
source("No cost/tf.r")

ca <- 400
k <- 0.05
MAP <- 365
parms <- c(LAI=1,
           Vcmax=50, cp=30, Km=703, Rd=1,
           a=1.6, nZ=0.5, p=43200, l=1.8e-5, VPD=0.02)

w0 <- 1e-3
gsw0 <- 1e-3
muf(-12, parms)
#muopt <- optimize(muf, c(-8, 0), tol=.Machine$double.eps, parms=parms)
#muopt
