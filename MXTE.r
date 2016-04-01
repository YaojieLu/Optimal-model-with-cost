
f1 <- function(gs,
               ca=400, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
               a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
               pe=-1.58*10^-3, b=4.38, h2=h/1000, xkmax=5, c=5.71, d=10.05){
  
  # photosynthesis rate function
  Af <- function(gs){LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-sqrt(((ca-Km)*gs+Rd-Vcmax)^2+4*gs*((ca*gs+Rd)*Km+cp*Vcmax)))}
  # xylem conductance function
  xkf <- function(x)xkmax*exp(-(-x/d)^c)
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)-(ps-x)*h2*xkf(x)
    res <- optimize(f1, c(-20,0), tol=.Machine$double.eps)$minimum
    return(res)
  }
  # xylem water potential function
  pxf <- function(w, gs){
    ps <- pe*w^(-b)
    pxmin <- pxminf(w)
    f1 <- function(x)((ps-x)*h2*xkf(x)-h*VPD*gs)^2
    res <- optimize(f1, c(pxmin, ps), tol=.Machine$double.eps)$minimum
    return(res)
  }
  # maximum gs function for given w
  gsmaxf <- function(w){
    ps <- pe*w^(-b)
    pxmin <- pxminf(w)
    gsmax <- (ps-pxmin)*h2*xkf(pxmin)/(h*VPD)
    return(gsmax)
  }
  # dA/dgs
  dAdgsf <- function(gs)(1/2)*LAI*(ca+Km-(4*ca*gs*Km+2*(ca-Km)*(gs*(ca-Km)+Rd-Vcmax)+4*(Km*(ca*gs+Rd)+cp*Vcmax))/(2*sqrt((gs*(ca-Km)+Rd-Vcmax)^2+4*gs*(Km*(ca*gs+Rd)+cp*Vcmax))))
  # dpx/dgs
  dpxdgsf <- function(w, gs){
    px <- pxf(w, gs)
    res <- h*VPD/(h*VPD*gs/d*c*(-px/d)^(c-1)-h2*xkmax*exp(-(-px/d)^c))
    return(res)
  }
  f1 <- function(gs)dAdgsf(gs)/-dpxdgsf(wtest, gs)
  #browser()
  res <- f1(gs)
  return(res)
}
f2 <- Vectorize(f1)

wtest <- 0.1
curve(f2, 0, 1, ylim=c(-100, 100))
wtest <- 0.2
curve(f2, 0, 1, col="darkgreen", add=T)
wtest <- 0.3
curve(f2, 0, 1, col="blue", add=T)
wtest <- 0.4
curve(f2, 0, 1, col="red", add=T)
