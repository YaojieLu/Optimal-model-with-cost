
# gsmax
gsmaxf1 <- function(w,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                    pe=-1.58*10^-3, b=4.38, h2=h/1000, xkmax=5, c=5.71, d=10.05){
  
  # xylem conductance function
  xkf <- function(x)xkmax*exp(-(-x/d)^c)
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)-(ps-x)*h2*xkf(x)
    res <- optimize(f1, c(-200,0), tol=.Machine$double.eps)$minimum
    return(res)
  }
  
  ps <- pe*w^(-b)
  pxmin <- pxminf(w)
  res <- (ps-pxmin)*h2*xkf(pxmin)/(h*VPD)
  return(res)
}
gsmaxf2 <- Vectorize(gsmaxf1)
# m(w, gs)
mf1 <- function(w, gs,
                ca=400, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
                a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                pe=-1.58*10^-3, b=4.38, h2=h/1000, xkmax=5, c=5.71, d=10.05, h3=10){
  
  # xylem conductance function
  xkf <- function(x)xkmax*exp(-(-x/d)^c)
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)-(ps-x)*h2*xkf(x)
    res <- optimize(f1, c(-200,0), tol=.Machine$double.eps)$minimum
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
  
  px <- pxf(w, gs)
  xk <- xkf(px)
  PLC <- 1-xk/xkmax
  res <- h3*PLC
  #browser()
  return(res)
}
# m(w, 0)
mgs0f1 <- function(w, pe=-1.58*10^-3, b=4.38, xkmax=5, c=5.71, d=10.05, h3=10){
  
  # xylem conductance function
  xkf <- function(x)xkmax*exp(-(-x/d)^c)
  
  px <- pe*w^(-b)
  xk <- xkf(px)
  PLC <- 1-xk/xkmax
  res <- h3*PLC
  #browser()
  return(res)
}
# A(gs)
Af <- function(gs, ca=400, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-sqrt(((ca-Km)*gs+Rd-Vcmax)^2+4*gs*((ca*gs+Rd)*Km+cp*Vcmax)))
# B(w, gs)
Bf1 <- function(w, gs)Af(gs)-mf1(w, gs)
# wtest
wtest <- 0.5
Bf2 <- Vectorize(function(gs)Bf1(wtest, gs))
mf2 <- Vectorize(function(gs)mf1(wtest, gs))
mf3 <- Vectorize(function(w)mf1(w, 0))
mgs0f2 <- function(w)(mf3(w)-mgs0f1(w))
gsmax <- gsmaxf2(wtest)
# Figure
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 2.5, 1), mfrow=c(1,1))
curve(gsmaxf2, 1e-5, 1,
      main=expression(italic(g[smax](w))),
      xlab=expression(italic(w)), ylab=expression(italic(g[smax])~(mol~m^-2~s^-1)),
      type="l", xlim=c(0, 1), ylim=c(0, 2), cex.lab=1.3)
curve(mgs0f1, 1e-5, 1,
      main=expression(italic(m(w,0))),
      xlab=expression(italic(w)), ylab=expression(italic(m)~(mu*mol~m^-2~s^-1)),
      type="l", xlim=c(0, 1), ylim=c(0, 12), cex.lab=1.3)
abline(v=0.1350649598204598, col="red")
curve(mf2, 1e-5, gsmax, type="l", xlim=c(0, gsmax), cex.lab=1.3)
curve(Bf2, 1e-5, gsmax, type="l", xlim=c(0, gsmax), cex.lab=1.3)
