

# psf(w)
psf <- function(w, pe=-1.58*10^-3, b=4.38)pe*w^(-b)

# wf(ps)
wf <- function(ps, pe=-1.58*10^-3, b=4.38)(ps/pe)^(-1/b)

# PLC(px)
PLCf <- function(px, c=2.64, d=3.54)1-exp(-(-px/d)^c)

# xylem conductance function
kxf <- function(px, kxmax=5, c=2.64, d=3.54)kxmax*exp(-(-px/d)^c)

# minimum xylem water potential function at given w
pxminf <- function(w, a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, h2=l*LAI/nZ*p/1000){
  ps <- psf(w)
  f1 <- function(px)(ps-px)*h2*kxf(px)
  res <- optimize(f1, c(-20, 0), tol=.Machine$double.eps, maximum=T)$maximum
  return(res)
}

# Water balance in plants
gswpxf <- function(w, px, a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, h2=l*LAI/nZ*p/1000, VPD=0.02){
  ps <- psf(w)
  res <- (ps-px)*h2*kxf(px)/(h*VPD)
  return(res)
}

# Inverse pxmin(w)
inversepxminf <- function(px)optimize(function(w)(pxminf(w)-px)^2, c(0, 1), tol=.Machine$double.eps)$minimum

# Water balance in plants with pxmin(w)
gspxminf <- Vectorize(function(px)gswpxf(inversepxminf(px), px))

# gsmaxf(w)
gsmaxf <- function(w, a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02, h2=l*LAI/nZ*p/1000){
  ps <- psf(w)
  pxmin <- pxminf(w)
  res <- (ps-pxmin)*h2*kxf(pxmin)/(h*VPD)
  return(res)
}

# Inverse gsmaxf(w)
inversegsmaxf <- function(gs)optimize(function(w)(gsmaxf(w)-gs)^2, c(0, 1), tol=.Machine$double.eps)$minimum

# PLC with gsmaxf(w)
PLCgsmaxf <- Vectorize(function(gs)PLCwgsf(inversegsmaxf(gs), gs)*100)

# xylem water potential function
pxf <- function(w, gs, a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02, h2=l*LAI/nZ*p/1000){
  ps <- psf(w)
  pxmin <- pxminf(w)
  f1 <- function(px)((ps-px)*h2*kxf(px)-h*VPD*gs)^2
  res <- ifelse(pxmin<ps, optimize(f1, c(pxmin, ps), tol=.Machine$double.eps)$minimum, ps)
  return(res)
}

# Af(gs, ca)
Af <- function(gs, ca, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1){
  LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))
}

# PLCwgsf(w, gs)
PLCwgsf <- function(w, gs){
  px <- pxf(w, gs)
  res <- PLCf(px)
  return(res)
}

# PLCwgsf(w, gs) px < pxmin
PLCwgsf2 <- function(w, gs,
                     a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02, h2=l*LAI/nZ*p/1000){
  
  # xylem water potential function
  pxf1 <- function(w, gs){
    ps <- psf(w)
    pxmin <- pxminf(w)
    f1 <- function(x)((ps-x)*h2*kxf(x)-h*VPD*gs)^2
    res <- optimize(f1, c(-1000, pxmin), tol=.Machine$double.eps)$minimum
    return(res)
  }
  
  px <- pxf1(w, gs)
  res <- PLCf(px)
  return(res)
}

# mf(w, gs)
mf <- function(w, gs, h3=10)h3*PLCwgsf(w, gs)

# mf(w, gs) px < pxmin
mf2 <- function(w, gs, h3=10)h3*PLCwgsf2(w, gs)

# m with gsmaxf(w)
mgsmaxf <- Vectorize(function(gs)mf(inversegsmaxf(gs), gs))

# B(w, gs)
Bf <- function(w, gs, ca)Af(gs, ca)-mf(w, gs)

# ESS gs(w)
ESSf <- function(w, ca){
  f1 <- function(gs)Bf(w, gs, ca)
  res <- optimize(f1, c(0, gsmaxf(w)), tol=.Machine$double.eps, maximum=T)
  return(res$maximum)
}

# ESS gs(ps)
ESSpsf <- function(ps, ca){
  w <- wf(ps)
  res <- ESSf(w, ca)
  return(res)
}

# ESS A(w)
ESSAf <- function(w, ca)Af(ESSf(w, ca), ca)

# ESS m(w)
ESSmf <- function(w, ca)mf(w, ESSf(w, ca))

# ESS B(w)
ESSBf <- function(w, ca)ESSAf(w, ca)-ESSmf(w, ca)

# ESS B(ps)
ESSBpsf <- function(ps, ca){
  w <- wf(ps)
  res <- Bf(w, ESSf(w, ca), ca)
  return(res)
}

# ESS PLC(w)
ESSPLCf <- function(w, ca){
  px <- pxf(w, ESSf(w, ca))
  res <- PLCf(px)
  return(res)
}

# ESS PLC(ps)
ESSPLCpsf <- function(ps, ca)ESSPLCf(wf(ps), ca)

# integralfnoc of PDF
integralfnocf <- function(ca, k, MAP, wL,
                          LAI=1, a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                          gamma=1/((MAP/365/k)/1000)*nZ){
  
  ESSf1 <- Vectorize(function(w)ESSf(w, ca))
  Ef <- function(w){h*VPD*ESSf1(w)}
  rEf <- function(w){1/Ef(w)}
  integralrEf <- Vectorize(function(w){integrate(rEf, wL, w, rel.tol=.Machine$double.eps^0.3)$value})#
  fnoc <- function(w){1/Ef(w)*exp(-gamma*(w-wL)/(1-wL)+k*integralrEf(w)*1/(1-wL))*1/(1-wL)}
  res <- integrate(fnoc, wL, 1, rel.tol=.Machine$double.eps^0.3)
  return(res)
}

# PDF of w
PDFf <- function(w, ca, k, MAP, wL, cPDF,
                 LAI=1, a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                 gamma=1/((MAP/365/k)/1000)*nZ){
  
  ESSf1 <- Vectorize(function(w)ESSf(w, ca))
  Ef <- function(w){h*VPD*ESSf1(w)}
  rEf <- function(w){1/Ef(w)}
  integralrEf <- Vectorize(function(w){integrate(rEf, wL, w, rel.tol=.Machine$double.eps^0.3)$value})#
  res <- cPDF/Ef(w)*exp(-gamma*(w-wL)/(1-wL)+k*integralrEf(w)*1/(1-wL))*1/(1-wL)
  return(res)
}

# Average A
averAf <- function(ca, k, MAP, wL, cPDF){
  ESSAf1 <- Vectorize(function(w)ESSAf(w, ca))
  f1 <- function(w)ESSAf1(w)*PDFf(w, ca, k, MAP, wL, cPDF)
  res <- integrate(f1, wL, 1, rel.tol=.Machine$double.eps^0.3)
  return(res)
}

# Average m
avermf <- function(ca, k, MAP, wL, cPDF){
  ESSmf1 <- Vectorize(function(w)ESSmf(w, ca))
  f1 <- function(w)ESSmf1(w)*PDFf(w, ca, k, MAP, wL, cPDF)
  res <- integrate(f1, wL, 1, rel.tol=.Machine$double.eps^0.3)
  return(res)
}

# Average E
averEf <- function(ca, k, MAP, wL, cPDF,
                   LAI=1, a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02){
  ESSf1 <- Vectorize(function(w)ESSf(w, ca))
  Ef <- function(w)h*VPD*ESSf1(w)
  f1 <- function(w)Ef(w)*PDFf(w, ca, k, MAP, wL, cPDF)
  res <- integrate(f1, wL, 1, rel.tol=.Machine$double.eps^0.3)
  return(res)
}

# Average w
averwp1f <- function(ca, k, MAP, wL, cPDF){
  f1 <- function(w)w*PDFf(w, ca, k, MAP, wL, cPDF)
  res <- integrate(f1, wL, 1, rel.tol=.Machine$double.eps^0.3)
  return(res)
}

# Average cica
avercicaf <- function(ca, k, MAP, wL, cPDF){
  ESSAf1 <- Vectorize(function(w)ESSAf(w, ca))
  ESSf1 <- Vectorize(function(w)ESSf(w, ca))
  f1 <- function(w)1-ESSAf1(w)/ESSf1(w)/ca
  f2 <- function(w)f1(w)*PDFf(w, ca, k, MAP, wL, cPDF)
  
  res <- integrate(f2, wL, 1, rel.tol=.Machine$double.eps^0.3)
  return(res)
}

# Initializing
ca <- 400
ESSf1 <- Vectorize(function(w)ESSf(w, ca))
ESSBf1 <- Vectorize(function(w)ESSBf(w, ca))

wL <- uniroot(ESSBf1, c(0.12, 1), tol=.Machine$double.eps)$root

curve(ESSf1, wL, 1, add=T, col="blue")
