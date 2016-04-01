
source("test (functions).r")

ca1 <- 400
k1 <- 0.025
MAP1 <- 365

w0opt <- w0optf(ca1)
w0ESS <- optimize(w0ESSf, c(0.113, 1), ca=ca1, tol=.Machine$double.eps)$minimum
mu <- muf(ca1, k1, MAP1)

disf2 <- Vectorize(function(w)disf(w, ca1, k1, MAP1))
xf <- function(w)disf2(w)-ESSf2(w)
curve(xf, 0.16, 0.19)

optf2 <- Vectorize(function(w)optf1(w, ca1, k1, MAP1))
curve(optf2, w0opt, 1)
optf2(0.18)
curve(optf2, w0opt, 0.17)
curve(optf2, 0.187, 1)
curve(optf2, 0.171, 0.186)
