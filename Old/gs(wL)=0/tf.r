
# psf(w)
psf <- function(w, pe=-1.58*10^-3, b=4.38)pe*w^(-b)

# PLC(px)
PLCf <- function(px, c=2.64, d=3.54)1-exp(-(-px/d)^c)

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

# mf(w, gs)
mf <- function(w, gs, h3=10)h3*PLCwgsf(w, gs)

# B(w, gs)
Bf <- function(w, gs, ca)Af(gs, ca)-mf(w, gs)

# dgs/dw
dgsdwf <- function(w, y, parms, mu, wL){
  with(as.list(c(y, parms, mu, wL)), {
    
    h <- l*a*LAI/nZ*p
    gamma <- 1/((MAP/365/k)/1000)*nZ
    h2 <- l*LAI/nZ*p/1000
    
    px <- pxf(w, gs)
    res <- ((2*b*(-1+c)*c*gs*h*h3*pe*(-(px/d))^c*VPD*w^(-1+b))/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)-(4*b*c^2*gs*h*h3*pe*(-(px/d))^(2*c)*VPD*w^(-1+b))/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(2*b*c*h3*pe*(-(px/d))^c)/(exp((-(px/d))^c)*(w*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))))-(2*gamma*((-1+exp(-(-(px/d))^c))*h3+mu-(1/2)*LAI*((-ca)*gs-gs*Km+Rd-Vcmax+sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))))/(1-wL)-(gs*(-gamma+k/(gs*h*VPD))*(LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))+(2*c*h*h3*(-(px/d))^c*VPD*w^b)/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))))/(1-wL))/(gs*h^2*VPD^2*((LAI*((2*(ca+Km))/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)^2/((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))^(3/2)+(-3*ca^2*gs-3*gs*Km^2-2*ca*(3*gs*Km+Rd-Vcmax)-2*(2*cp*Vcmax+Km*(Rd+Vcmax)))/(gs*sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))))/(h^2*VPD^2)-(2*(-1+c)*c*exp((-(px/d))^c)*h3*(-(px/d))^c*w^(2*b))/(h2^2*kxmax^2*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(2*c^2*exp((-(px/d))^c)*h3*(-(px/d))^(2*c)*w^(2*b))/(h2^2*kxmax^2*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(4*c*h3*(-(px/d))^c*w^b)/(gs*h*h2*kxmax*VPD*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))-(2*(LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))+(2*c*h*h3*(-(px/d))^c*VPD*w^b)/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))))/(gs*h^2*VPD^2)))
    return(list(res))
  })
}


# Boundary condition at w=1
gs1f <- function(w, parms, mu, ca){
  with(as.list(c(w, parms, mu)), {
    
    h <- l*a*LAI/nZ*p
    h2=l*LAI/nZ*p/1000
    
    f1 <- function(gs){
      px <- pxf(w, gs)
      res <- (1/2)*gs*h*VPD*((-gs)*LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))+2*((-1+exp(-(-(px/d))^c))*h3+mu-(1/2)*LAI*((-ca)*gs-gs*Km+Rd-Vcmax+sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))))-(2*c*gs*h*h3*(-(px/d))^c*VPD*w^b)/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))))
      return(res^2)
    }
    
    gsmax <- gsmaxf(w)
    res <- optimize(f1, c(0.02, gsmax*0.99), tol=.Machine$double.eps)
    return(res)
  })
}

# optimize mu
muf <- function(mu, wL){
  gswL <- pp
  out <- ode(y=c(gs=gswL), times=c(wL, 1), func=dgsdwf, parms=parms, mu=mu, wL=wL)
  gs11 <- out[2, 2]
  gs12 <- gs1f(1, parms, mu, ca)
  res <- gs11-gs12$minimum
  return(res^2)
}

# gs(w)
gswf <- function(w, mu, wL){
  gswL <- pp
  res <- ode(y=c(gs=gswL), times=c(wL, w), func=dgsdwf, parms=parms, mu=mu, wL=wL)[2, 2]
  return(res)
}

# integralfnoc of PDF
integralfnocf <- function(wL, parms){
  with(as.list(c(wL, parms)), {
    
    h <- l*a*LAI/nZ*p
    gamma <- 1/((MAP/365/k)/1000)*nZ
    
    Ef <- function(w){h*VPD*gswf1(w)}
    rEf <- function(w){1/Ef(w)}
    integralrEf <- Vectorize(function(w){integrate(rEf, w, 1, rel.tol=.Machine$double.eps^0.25)$value})
    fnoc <- function(w){1/Ef(w)*exp(-gamma*(w-wL)/(1-wL)-k*integralrEf(w)*1/(1-wL))*1/(1-wL)}
    
    res <- integrate(fnoc, wL, 1, rel.tol=.Machine$double.eps^0.25)
    return(res)
  })
}

# cPDF
cPDFf <- function(integralfnoc, wL, parms){
  with(as.list(c(integralfnoc, wL, parms)), {
    
    h <- l*a*LAI/nZ*p
    gamma <- 1/((MAP/365/k)/1000)*nZ
    
    Ef <- function(w){h*VPD*gswf1(w)}
    rEf <- function(w){1/Ef(w)}
    res <- 1/(integralfnoc+1/k*exp(-k/(1-wL)*integrate(rEf, wL, 1, rel.tol=.Machine$double.eps^0.25)$value))
    return(res)
  })
}

# PDF of w
PDFf <- function(w, wL, cPDF, parms){
  with(as.list(c(wL, cPDF, parms)), {
    
    h <- l*a*LAI/nZ*p
    gamma <- 1/((MAP/365/k)/1000)*nZ
    
    Ef <- function(w){h*VPD*gswf1(w)}
    rEf <- function(w){1/Ef(w)}
    integralrEf <- Vectorize(function(w){integrate(rEf, w, 1, rel.tol=.Machine$double.eps^0.25)$value})
    res <- cPDF/Ef(w)*exp(-gamma*(w-wL)/(1-wL)-k*integralrEf(w)*1/(1-wL))*1/(1-wL)
    return(res)
  })
}

# Average B
averBf <- function(wL, cPDF, parms){
  with(as.list(c(wL, cPDF, parms)), {
    
    h <- l*a*LAI/nZ*p
    gamma <- 1/((MAP/365/k)/1000)*nZ
    
    Bf1 <- Vectorize(function(w)Bf(w, gswf1(w), ca))
    f1 <- function(w)Bf1(w)*PDFf(w, wL, cPDF, parms)
    
    res <- integrate(f1, wL, 1, rel.tol=.Machine$double.eps^0.25)
    return(res)
  })
}

# Average E/MAP
averEMAPf <- function(wL, cPDF, parms){
  with(as.list(c(wL, cPDF, parms)), {
    
    h <- l*a*LAI/nZ*p
    gamma <- 1/((MAP/365/k)/1000)*nZ
    
    Ef <- function(w){h*VPD*gswf1(w)}
    f1 <- function(w)Ef(w)*PDFf(w, wL, cPDF, parms)
    
    res <- integrate(f1, wL, 1, rel.tol=.Machine$double.eps^0.25)$value*nZ*1000*365/MAP
    return(res)
  })
}

# Average w
averwf <- function(wL, cPDF, parms, integralfnoc){
  with(as.list(c(wL, cPDF, parms)), {
    
    h <- l*a*LAI/nZ*p
    gamma <- 1/((MAP/365/k)/1000)*nZ
    
    f1 <- function(w)w*PDFf(w, wL, cPDF, parms)
    res <- integrate(f1, wL, 1, rel.tol=.Machine$double.eps^0.25)$value+wL*(1-cPDF*integralfnoc)
    return(res)
  })
}
