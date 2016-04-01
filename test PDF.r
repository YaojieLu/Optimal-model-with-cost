
source("test (functions).r")

ca1 <- 400
k1 <- 0.025
MAP1 <- 365

w0opt <- w0optf(ca1)
mu <- muf(ca1, k1, MAP1)

# Functions
PDFf <- function(ca, k, MAP,
                 Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
                 a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                 pe=-1.58*10^-3, b=4.38, h2=h/1000, xkmax=5, c=5.71, d=10.05, h3=10,
                 MDP=MAP/365, gamma=1/(MDP/k/1000)*nZ){
  
  optf2 <- Vectorize(function(w)optf1(w, ca, k, MAP))
  Ev <- function(w){h*VPD*optf2(w)}
  rEv <- function(w){1/Ev(w)}
  integralrEv <- Vectorize(function(w){integrate(rEv, w, 1, rel.tol=.Machine$double.eps^0.5, subdivisions = 10000)$value})
  fnoc <- function(w){1/Ev(w)*exp(-gamma*w-k*integralrEv(w))}
  #browser()
  integralfnoc <- integrate(fnoc, w0opt, 1, rel.tol=.Machine$double.eps^0.5)$value
  c <- 1/integralfnoc
  message(c)
  f1 <- function(w){
    res <- c*fnoc(w)
    message(w)
    return(res)}
  
  x <- seq(w0opt, 1, by=(1-w0opt)/101)
  res <- f1(x)
  return(res)
}
y <- PDFf(ca1, k1, MAP1)

# Figure
#windows(8, 6)
#par(mgp = c(2.2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(4, 4, 2.5, 4), mfrow=c(1,1))
