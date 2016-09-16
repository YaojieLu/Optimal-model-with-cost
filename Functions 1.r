
# psf(w)
psf <- function(w, pe=-1.58*10^-3, b=4.38)pe*w^(-b)

# PLC(px)
PLCf <- function(px, c=2.64, d=3.54)1-exp(-(-px/d)^c)

# xylem conductance function
kxf <- function(px, kxmax=5)kxmax*(1-PLCf(px))

# minimum xylem water potential function at given w
pxminf <- function(w){
  ps <- psf(w)
  f1 <- function(px)(ps-px)*kxf(px)
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
  res <- ifelse(gs==0, ps, optimize(f1, c(pxmin, ps), tol=.Machine$double.eps)$minimum)
  return(res)
}

# Af(gs)
Af <- function(gs, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))

# PLCwgsf(w, gs)
PLCwgsf <- function(w, gs){
  px <- pxf(w, gs)
  res <- PLCf(px)
  return(res)
}

# mf(w, gs)
mf <- function(w, gs)h3*PLCwgsf(w, gs)

# B(w, gs)
Bf <- function(w, gs)Af(gs)-mf(w, gs)

# Optimize wL
wLf <- Vectorize(function(wL){
  message(wL)
  averB <- averBf(wL, parms)
  message(wL, " ", averB)
  return(averB)
})

# Curve fit
CFf <- function(pars,
                wL, gswL,
                LAI=1,
                Vcmax=50, cp=30, Km=703, Rd=1,
                a=1.6, nZ=0.5, p=43200, l=1.8e-5, VPD=0.02,
                pe=-1.58*10^-3, b=4.38, kxmax=5,
                h=l*a*LAI/nZ*p, h2=l*LAI/nZ*p/1000, gamma=1/((MAP/365/k)/1000)*nZ){
  
  p1 <- abs(pars[1]) # we require p1 be nonnegativ
  p2 <- pars[2]
  p3 <- abs(pars[3]) # we require p3 be nonnegativ
  p4 <- pars[4]
  gswf <- function(w)p1*(w-wL)^p2+p3*(w-wL)^p4+gswL
  
  Ef <- function(w)h*VPD*gswf(w)
  rEf <- function(w)1/Ef(w)
  integralrEf <- Vectorize(function(w)integrate(rEf, wL, w, rel.tol=.Machine$double.eps^0.25)$value)
  fnoc <- function(w)1/Ef(w)*exp(-gamma*w+k*integralrEf(w))
  
  f <- function(w)cPDF*fnoc(w)
  Bf1 <- function(w)Bf(w, gswf(w))
  fB <- Vectorize(function(w)f(w)*Bf1(w))
  
  integralfnoc <- integrate(fnoc, wL, 1, rel.tol=.Machine$double.eps^0.25)
  cPDF <- 1/(integralfnoc$value+1/k*exp(-gamma*wL))
  averB <- integrate(fB, wL, 1, rel.tol=.Machine$double.eps^0.25)
  return(averB$value)
}
