
options(digits=20)
library(deSolve)
source("tf.r")

ca <- 400
k <- 0.025
MAP <- 1000
parms <- c(LAI=1,
           Vcmax=50, cp=30, Km=703, Rd=1,
           a=1.6, nZ=0.5, p=43200, l=1.8e-5, VPD=0.02,
           pe=-1.58*10^-3, b=4.38, c=2.64, d=3.54, kxmax=5, h3=10)

wL <- 0.2
#optimize(muf, c(-20, 0), tol=.Machine$double.eps, wL=wL)$minimum
mu <- optimize(muf, c(-20, 0), tol=.Machine$double.eps, wL=wL)$minimum
gsBCf(wL, parms, mu, ca)
gsBCf(1, parms, mu, ca)
gswf1 <- Vectorize(function(w)gswf(w, mu, wL))
curve(gswf1, wL, 1)