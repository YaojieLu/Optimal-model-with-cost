
source("Functions.r")

# w0opt & w0ESS
ca1 <- 400
w0opt <- w0optf(ca1)
w0ESS <- optimize(w0ESSf, c(0.113, 1), ca=ca1, tol=.Machine$double.eps)$minimum
ESSBf2 <- Vectorize(function(w)ESSBf1(w, ca1))

# Figure
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 2.5, 1), mfrow=c(1,1))
curve(ESSBf2, w0ESS, 1, col="black", type="l",
      xlab=expression(italic(w)),
      ylab=expression(italic(B)~(mu*mol~m^-2~s^-1)),
      xlim=c(0, 1), ylim=c(0, 20),
      cex.lab=1.3)
#abline(v=w0ESS, col="black")
#abline(v=w0opt, col="red")
#legend("topright", c("ESS", "Optimal"), lty=c(1, 1), col=c("black", "red"))

# Scenarios
####################################################################
k1 <- 0.025
####################################################################
MAP1 <- 365
mu <- muf(ca1, k1, MAP1)
optBf2 <- Vectorize(function(w)optBf1(w, ca1, k1, MAP1))
curve(optBf2, w0opt, 1, lty=1, col="red", add=T)
####################################################################
MAP1 <- 3650
mu <- muf(ca1, k1, MAP1)
optBf2 <- Vectorize(function(w)optBf1(w, ca1, k1, MAP1))
curve(optBf2, w0opt, 1, lty=2, col="red", add=T)
####################################################################
k1 <- 0.05
####################################################################
MAP1 <- 365
mu <- muf(ca1, k1, MAP1)
optBf2 <- Vectorize(function(w)optBf1(w, ca1, k1, MAP1))
curve(optBf2, w0opt, 1, lty=1, col="blue", add=T)
####################################################################
MAP1 <- 3650
mu <- muf(ca1, k1, MAP1)
optBf2 <- Vectorize(function(w)optBf1(w, ca1, k1, MAP1))
curve(optBf2, w0opt, 1, lty=2, col="blue", add=T)
####################################################################
k1 <- 0.1
####################################################################
MAP1 <- 365
mu <- muf(ca1, k1, MAP1)
optBf2 <- Vectorize(function(w)optBf1(w, ca1, k1, MAP1))
curve(optBf2, w0opt, 1, lty=1, col="darkgreen", add=T)
####################################################################
MAP1 <- 3650
mu <- muf(ca1, k1, MAP1)
optBf2 <- Vectorize(function(w)optBf1(w, ca1, k1, MAP1))
curve(optBf2, w0opt, 1, lty=2, col="darkgreen", add=T)
####################################################################
