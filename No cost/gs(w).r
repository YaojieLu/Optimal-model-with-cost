
######################### No cost #########################

options(digits=20)
library(deSolve)
source("No cost/Functions.r")

ca <- 400
k <- 0.05
MAP <- 1000
parms <- c(LAI=1,
           Vcmax=50, cp=30, Km=703, Rd=1,
           a=1.6, nZ=0.5, p=43200, l=1.8e-5, VPD=0.02)
w0 <- 1e-5

#muopt <- uniroot(muf, c(-20, 0), tol=.Machine$double.eps)
#muopt
#f <- Vectorize(function(w)gswf(w, muopt$root))
#
## Figure
#windows(8, 6)
#par(mgp = c(2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(3.5, 4, 1, 2))
#curve(f, 0, 1,
#     xlab = expression(italic(w)~"(%)"), 
#     ylab = expression(italic(g[s])~(mol~m^-2~s^-1)),
#     xlim = c(0, 1), ylim = c(0, 0.4),
#     cex.lab = 1.3,
#     col = "red"
#)