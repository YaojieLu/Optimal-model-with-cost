
# Curve fit - gs(w) as a double power function
options(digits=20)
library(deSolve)
library(optimx)
source("Functions 1.r")

ca <- 400
k <- 0.05
MAP <- 1000
h3 <- 10
c <- 2.64
d <- 3.54

gswLf <- function(gswL){
  int <- c(0.07, 4, 0.08, 0.5)
  res <- optimx(int, CFf, itnmax=5000, method="BFGS", control=list(maximize=T), wL=wL, gswL=gswL)
  message(gswL, " ", res$value)
  return(res$value)
}

wL <- 0.24
gswL <- optimize(gswLf, c(0, 0.02), tol=.Machine$double.eps, maximum=T)
gswL
