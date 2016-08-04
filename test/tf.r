
# psf(w)
psf <- function(w, pe=-1.58*10^-3, b=4.38)pe*w^(-b)

# PLC(px)
PLCf <- function(px)1-exp(-(-px/d)^c)

# xylem conductance function
kxf <- function(px, kxmax=5)kxmax*(1-PLCf(px))

# minimum xylem water potential function at given w
pxminf <- function(w, a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, h2=l*LAI/nZ*p/1000){
  ps <- psf(w)
  f1 <- function(px)(ps-px)*h2*kxf(px)
  res <- optimize(f1, c(-20, 0), tol=.Machine$double.eps, maximum=T)$maximum
  return(res)
}

# gsmaxf(w)
gsmaxf <- Vectorize(function(w, a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02, h2=l*LAI/nZ*p/1000){
  ps <- psf(w)
  pxmin <- pxminf(w)
  res <- (ps-pxmin)*h2*kxf(pxmin)/(h*VPD)
  return(res)
})

# xylem water potential function
pxf <- function(w, gs, a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02, h2=l*LAI/nZ*p/1000){
  ps <- psf(w)
  pxmin <- pxminf(w)
  f1 <- function(px)((ps-px)*h2*kxf(px)-h*VPD*gs)^2
  res <- optimize(f1, c(pxmin, ps), tol=.Machine$double.eps)$minimum
  return(res)
}

# Af(gs, ca)
Af <- function(gs, ca, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))

# PLCwgsf(w, gs)
PLCwgsf <- function(w, gs){
  px <- pxf(w, gs)
  res <- PLCf(px)
  return(res)
}

# mf(w, gs)
mf <- function(w, gs)h3*PLCwgsf(w, gs)

# B(w, gs)
Bf <- function(w, gs, ca)Af(gs, ca)-mf(w, gs)

# Optimize wL
wLf <- Vectorize(function(wL){
  averB <- averBf(wL, parms, 1)
  message(wL, " ", averB)
  return(averB)
})
