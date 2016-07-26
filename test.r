
options(digits=20)
library(deSolve)
source("tf.r")

wL <- 0.19
mul <- 1

ca <- 400
k <- 0.05
MAP <- 1000
parms <- c(LAI=1,
           Vcmax=50, cp=30, Km=703, Rd=1,
           a=1.6, nZ=0.5, p=43200, l=1.8e-5, VPD=0.02,
           pe=-1.58*10^-3, b=4.38, c=2.64, d=3.54, kxmax=5, h3=15)

#muf1 <- Vectorize(function(mu)muf(mu, wL))
#curve(muf1, -13, -11)
mu <- optimize(muf, c(-20, 0), tol=.Machine$double.eps, wL=wL)
mu

gswf1 <- Vectorize(function(w)gswf(w, mu$minimum, wL)*mul)
curve(gswf1, wL, 1, ylim=c(0, 0.2))
#f1(wL, parms, mu$minimum, ca)
#gswf1(wL)
#f1(1, parms, mu$minimum, ca)
#gswf1(1)

integralfnoc <- integralfnocf(wL, parms)
cPDF <- cPDFf(integralfnoc$value, wL, parms)
averB <- averBf(wL, cPDF, parms)

integralfnoc
cPDF
averB
