
options(digits=20)
library(deSolve)
source("Functions.r")

# Initialization
ca <- 400
c <- 2.64
d <- 3.54
parms <- c(LAI=1,
           Vcmax=50, cp=30, Km=703, Rd=1,
           a=1.6, nZ=0.5, p=43200, l=1.8e-5, VPD=0.02,
           pe=-1.58*10^-3, b=4.38, kxmax=5)

# Figure
windows(8, 6)
par(mgp = c(2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(3.5, 4, 1, 2))
plot(0, type="n",
     xlab = expression(italic(w)), 
     ylab = expression(italic(g[s])~(mol~m^-2~s^-1)),
     xlim = c(0, 1), ylim = c(0, 0.2),
     cex.lab = 1.3)

# Results
k <- 0.025
MAP <- 1000
h3 <- 10

wL <- 0.199779892528234
mu <- muf(wL, parms)
gswf1 <- Vectorize(function(w)gswf(w, wL, mu, parms))
curve(gswf1, wL, 1, add=T)

# Results
k <- 0.05
MAP <- 1000
h3 <- 10

wL <- 0.202228472749335
mu <- muf(wL, parms)
gswf1 <- Vectorize(function(w)gswf(w, wL, mu, parms))
curve(gswf1, wL, 1, add=T, col="red")

# Results
k <- 0.1
MAP <- 1000
h3 <- 10

wL <- 0.203561175249867
mu <- muf(wL, parms)
gswf1 <- Vectorize(function(w)gswf(w, wL, mu, parms))
curve(gswf1, wL, 1, add=T, col="blue")
