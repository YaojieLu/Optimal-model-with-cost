
options(digits=20)
library(deSolve)
library(optimx)
source("Functions.r")

ca <- 400
k <- 0.025
MAP <- 1000
h3 <- 10
c <- 2.64
d <- 3.54

wL <- 0.199779892528234
gswL <- 0.010014785901787745

# Curve fit - gs(w) as a double power function
int <- c(0.07, 4, 0.08, 0.5)
#CFf_wrapper(int)
pars <- optimx(int, CFf_wrapper, itnmax=5000, method="BFGS", control=list(maximize=T))
gswf <- function(w){abs(pars$p1)*(w-wL)^pars$p2+abs(pars$p3)*(w-wL)^pars$p4+gswL}

# Figure
windows(8, 6)
par(mgp = c(2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(3.5, 4, 1, 2))
curve(gswf, wL, 1,
      xlim=c(0, 1), ylim=c(0, 0.2),
      xlab = expression(italic(w)), 
      ylab = expression(italic(g[s])~(mol~m^-2~s^-1)),
      cex.lab = 1.3)
