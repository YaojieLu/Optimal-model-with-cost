
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

wLf <- function(wL){
  gswLf <- function(gswL){
    int <- c(0.07, 4, 0.08, 0.5)
    f1 <- function(pars)-CFf(pars, wL=wL, gswL=gswL)
    res <- optimx(int, f1, itnmax=5000, method="BFGS")
    return(res$value)
  }
  
  res <- optimize(gswLf, c(0.001, 0.02), tol=.Machine$double.eps)
  message(wL, " ", res[1], " ",res[2])
  return(res$objective)
}

wL <- optimize(wLf, c(0.14, 1), tol=.Machine$double.eps)
