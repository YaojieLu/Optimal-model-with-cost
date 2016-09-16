
# Curve fit - gs(w) as a double power function
options(digits=20)
library(deSolve)
library(optimx)
source("Functions 1.r")

# Parameters
ca <- 400
k <- 0.05
MAP <- 1000
h3 <- 10
c <- 2.64
d <- 3.54

# Inputs
wL <- 0.24
gswL <- 0.0094749006994713589

# Optimization
int <- c(0.07, 4, 0.08, 0.5)
#CFf(int, wL, gswL)
pars <- optimx(int, CFf, itnmax=5000, method="BFGS", control=list(maximize=T), wL=wL, gswL=gswL)
pars
gswf <- function(w){abs(pars$p1)*(w-wL)^pars$p2+abs(pars$p3)*(w-wL)^pars$p4+gswL}

# Figure
#windows(8, 6)
par(mgp = c(2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(3.5, 4, 1, 2))
curve(gswf, wL, 1,
      xlim=c(0, 1), ylim=c(0, 0.2),
      xlab = expression(italic(w)), 
      ylab = expression(italic(g[s])~(mol~m^-2~s^-1)),
      cex.lab = 1.3)
