
# gsmax
gsmaxf <- function(w,
                   a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                   pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=5.71, d=10.05){
  
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)-(ps-x)*h2*kxf(x)
    res <- optimize(f1, c(-200,0), tol=.Machine$double.eps)$minimum
    return(res)
  }
  
  ps <- pe*w^(-b)
  pxmin <- pxminf(w)
  res <- (ps-pxmin)*h2*kxf(pxmin)/(h*VPD)
  return(res)
}
# m(w, gs)
mf <- function(w, gs,
               ca=400, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
               a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
               pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=5.71, d=10.05, h3=10){
  
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)-(ps-x)*h2*kxf(x)
    res <- optimize(f1, c(-200,0), tol=.Machine$double.eps)$minimum
    return(res)
  }
  # xylem water potential function
  pxf <- function(w, gs){
    ps <- pe*w^(-b)
    pxmin <- pxminf(w)
    f1 <- function(x)((ps-x)*h2*kxf(x)-h*VPD*gs)^2
    res <- optimize(f1, c(pxmin, ps), tol=.Machine$double.eps)$minimum
    return(res)
  }
  
  px <- pxf(w, gs)
  kx <- kxf(px)
  PLC <- 1-kx/kxmax
  res <- h3*PLC
  #browser()
  return(res)
}
# A(gs)
Af <- function(gs, ca=400, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-sqrt(((ca-Km)*gs+Rd-Vcmax)^2+4*gs*((ca*gs+Rd)*Km+cp*Vcmax)))
# B(w, gs)
Bf <- function(w, gs)Af(gs)-mf(w, gs)
# B/gs
Bgsf <- function(w, gs)Bf(w, gs)/gs
forgsf1 <- function(gs)Bgsf(0.2, gs)
forwf1 <- function(w)Bgsf(w, 0.01)

# Figures
Bf1 <- Vectorize(function(gs)Bf(1, gs))
curve(Bf1, 0.9, 1.1)
