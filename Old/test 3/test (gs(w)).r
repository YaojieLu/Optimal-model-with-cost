
options(digits=20)
source("test (Functions).r")

ca1 <- 400
k1 <- 0.025
MAP1 <- 365

wL <- 0.25
gsL <- B0f(wL, ca1)
mu1 <- muf(ca1, k1, MAP1)
mu <- mu1$minimum

source("test (Functions).r")
ESSBf2 <- Vectorize(function(w)ESSBf1(w, ca1, k1, MAP1))
optf2 <- Vectorize(function(w)optf1(w, ca1, k1, MAP1))
optBf2 <- Vectorize(function(w)optBf1(w, ca1, k1, MAP1))
optWUEf2 <- Vectorize(function(w)optWUEf1(w, ca1, k1, MAP1))

optf2(wL)
optBf2(wL)

curve(optf2, wL, 1)
curve(optWUEf2, wL, 1)
