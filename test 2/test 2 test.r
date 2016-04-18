
options(digits=20)
source("test 2/test 2 (functions).r")

ca1 <- 400
k1 <- 0.025
MAP1 <- 365

LHSopt <- LHSoptf(ca1)
wLopt <- LHSopt[1]
gsL <- LHSopt[2]
mu <- muf(ca1, k1, MAP1)

source("test 2/test 2 (functions).r")
ESSBf2 <- Vectorize(function(w)ESSBf1(w, ca1, k1, MAP1))
optf2 <- Vectorize(function(w)optf1(w, ca1, k1, MAP1))
optBf2 <- Vectorize(function(w)optBf1(w, ca1, k1, MAP1))
optWUEf2 <- Vectorize(function(w)optWUEf1(w, ca1, k1, MAP1))

optf2(wLopt)
optBf2(wLopt)

#curve(optf2, wLopt, 1)
#curve(optBf2, wLopt, 1)
#curve(optWUEf2, wLopt, 1)
curve(optf2, wLopt, 0.2)
curve(optBf2, wLopt, 0.2)
curve(optWUEf2, wLopt, 0.2)
