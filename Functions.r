
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

# Maximum averB
averBf <- function(wL, parms, mul){
  with(as.list(c(wL, parms, mul)), {
    h <- l*a*LAI/nZ*p
    h2 <- l*LAI/nZ*p/1000
    gamma <- 1/((MAP/365/k)/1000)*nZ
    
    # Functions
    # Boundary conditions
    gsBCf <- function(w, mu, wL){
      f1 <- function(gs){
        px <- pxf(w, gs)
        res <- (1/2)*gs*h*VPD*((-gs)*LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))+2*((-1+exp(-(-(px/d))^c))*h3+mu-(1/2)*LAI*((-ca)*gs-gs*Km+Rd-Vcmax+sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))))-(2*c*gs*h*h3*(-(px/d))^c*VPD*w^b)/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))))
        return(res^2)
      }
      
      gsmax <- gsmaxf(w)
      res <- optimize(f1, c(gsmax*0.1, gsmax), tol=.Machine$double.eps)$minimum
      return(res)
    }
    
    # dgs/dw
    dgsdwf <- function(w, y, parms, mu, wL){
      with(as.list(c(y, parms, mu, wL)), {
        px <- pxf(w, gs)
        res <- -((-((2*b*(-1+c)*c*gs*h*h3*pe*(-(px/d))^c*VPD*w^(-1+b))/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2))+(2*b*c^2*gs*h*h3*pe*(-(px/d))^(2*c)*VPD*w^(-1+b))/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)-(2*b*c*h3*pe*(-(px/d))^c)/(exp((-(px/d))^c)*(w*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))))+(2*b*c^2*gs*h*h3*pe*(-(px/d))^(2*c)*VPD*w^(-1+b)*(pe-c*pe+c*px*w^b))/(h2*kxmax*((-px)*w^b+c*(-(px/d))^c*(-pe+px*w^b))^3)+(2*gamma*((-1+exp(-(-(px/d))^c))*h3+mu-(1/2)*LAI*((-ca)*gs-gs*Km+Rd-Vcmax+sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))))/(1-wL)+(gs*(-gamma+k/(gs*h*VPD))*(LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))+(2*c*h*h3*(-(px/d))^c*VPD*w^b)/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))))/(1-wL))/(gs*h^2*VPD^2*((LAI*((2*(ca+Km))/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)^2/((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))^(3/2)+(-3*ca^2*gs-3*gs*Km^2-2*ca*(3*gs*Km+Rd-Vcmax)-2*(2*cp*Vcmax+Km*(Rd+Vcmax)))/(gs*sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))))/(h^2*VPD^2)-(2*(-1+c)*c*exp((-(px/d))^c)*h3*(-(px/d))^c*w^(2*b))/(h2^2*kxmax^2*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(2*c^2*exp((-(px/d))^c)*h3*(-(px/d))^(2*c)*w^(2*b))/(h2^2*kxmax^2*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)-(2*(LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))+(2*c*h*h3*(-(px/d))^c*VPD*w^b)/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))))/(gs*h^2*VPD^2)+(2*c*h3*(-(px/d))^c*w^b*(2*h2*kxmax*px^2*w^(2*b)+c^2*(-(px/d))^c*(pe-px*w^b)*((-exp((-(px/d))^c))*gs*h*(-1+(-(px/d))^c)*VPD*w^b+2*h2*kxmax*(-(px/d))^c*(pe-px*w^b))-c*(-(px/d))^c*w^b*(4*h2*kxmax*px*(-pe+px*w^b)+exp((-(px/d))^c)*gs*h*VPD*(pe+px*w^b))))/(gs*h*h2^2*kxmax^2*VPD*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^3))))
        return(list(res))
      })
    }
    
    # Optimize mu
    muf <- function(mu, wL){
      gswL <- gsBCf(wL, mu, wL)
      gs11 <- as.numeric(ode(y=c(gs=gswL), times=c(wL, 1), func=dgsdwf, parms=NULL, mu=mu, wL=wL, maxsteps=100000)[2, 2])#, method=c("ode45"), maxsteps=100000
      gs12 <- gsBCf(1, mu, wL)
      res <- gs11-gs12
      return(res^2)
    }
    
    # gs(w)
    gswf <- function(w, mu, wL){
      gswL <- gsBCf(wL, mu, wL)
      res <- ifelse(w==wL, gswL, as.numeric(ode(y=c(gs=gswL), times=c(wL, w), func=dgsdwf, parms=NULL, mu=mu, wL=wL, maxsteps=100000)[2, 2]))#, method=c("ode45"), maxsteps=100000
      return(res)
    }
    
    # integralfnoc of PDF
    integralfnocf <- function(wL){
      Ef <- function(w)h*VPD*gswf1(w)
      rEf <- function(w)1/Ef(w)
      integralrEf <- Vectorize(function(w)integrate(rEf, wL, w, rel.tol=.Machine$double.eps^0.25)$value)
      fnoc <- function(w)1/Ef(w)*exp(-gamma*(w-wL)/(1-wL)+k*1/(1-wL)*integralrEf(w))*1/(1-wL)
      
      res <- integrate(fnoc, wL, 1, rel.tol=.Machine$double.eps^0.25)
      return(res)
    }
    
    # PDF of w
    PDFf <- function(w, wL, cPDF){
      Ef <- function(w)h*VPD*gswf1(w)
      rEf <- function(w)1/Ef(w)
      integralrEf <- Vectorize(function(w)integrate(rEf, wL, w, rel.tol=.Machine$double.eps^0.25)$value)
      
      res <- cPDF/Ef(w)*exp(-gamma*(w-wL)/(1-wL)+k*1/(1-wL)*integralrEf(w))*1/(1-wL)
      return(res)
    }
    
    # Average B
    averBf <- function(wL, cPDF){
      Bf1 <- Vectorize(function(w)Bf(w, gswf1(w), ca))
      f1 <- function(w)Bf1(w)*PDFf(w, wL, cPDF)
      
      res <- integrate(f1, wL, 1, rel.tol=.Machine$double.eps^0.25)
      return(res)
    }
    
    mu <- optimize(muf, c(-20, 0), tol=.Machine$double.eps, wL=wL)$minimum
    gswf1 <- Vectorize(function(w)gswf(w, mu, wL)*mul)
    integralfnoc <- integralfnocf(wL)$value
    cPDF <- 1/(integralfnoc+1/k)
    averB <- averBf(wL, cPDF)$value
    return(averB)
  })
}

# Optimize wL
wLf <- Vectorize(function(wL){
  averB <- averBf(wL, parms, 1)
  message(wL, " ", averB)
  return(averB)
})

# mu
muf <- function(wL, parms){
  with(as.list(c(wL, parms)), {
    h <- l*a*LAI/nZ*p
    h2 <- l*LAI/nZ*p/1000
    gamma <- 1/((MAP/365/k)/1000)*nZ
    
    # Functions
    # Boundary conditions
    gsBCf <- function(w, mu, wL){
      f1 <- function(gs){
        px <- pxf(w, gs)
        res <- (1/2)*gs*h*VPD*((-gs)*LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))+2*((-1+exp(-(-(px/d))^c))*h3+mu-(1/2)*LAI*((-ca)*gs-gs*Km+Rd-Vcmax+sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))))-(2*c*gs*h*h3*(-(px/d))^c*VPD*w^b)/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))))
        return(res^2)
      }
      
      gsmax <- gsmaxf(w)
      res <- optimize(f1, c(gsmax*0.1, gsmax), tol=.Machine$double.eps)$minimum
      return(res)
    }
    
    # dgs/dw
    dgsdwf <- function(w, y, parms, mu, wL){
      with(as.list(c(y, parms, mu, wL)), {
        px <- pxf(w, gs)
        res <- -((-((2*b*(-1+c)*c*gs*h*h3*pe*(-(px/d))^c*VPD*w^(-1+b))/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2))+(2*b*c^2*gs*h*h3*pe*(-(px/d))^(2*c)*VPD*w^(-1+b))/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)-(2*b*c*h3*pe*(-(px/d))^c)/(exp((-(px/d))^c)*(w*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))))+(2*b*c^2*gs*h*h3*pe*(-(px/d))^(2*c)*VPD*w^(-1+b)*(pe-c*pe+c*px*w^b))/(h2*kxmax*((-px)*w^b+c*(-(px/d))^c*(-pe+px*w^b))^3)+(2*gamma*((-1+exp(-(-(px/d))^c))*h3+mu-(1/2)*LAI*((-ca)*gs-gs*Km+Rd-Vcmax+sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))))/(1-wL)+(gs*(-gamma+k/(gs*h*VPD))*(LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))+(2*c*h*h3*(-(px/d))^c*VPD*w^b)/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))))/(1-wL))/(gs*h^2*VPD^2*((LAI*((2*(ca+Km))/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)^2/((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))^(3/2)+(-3*ca^2*gs-3*gs*Km^2-2*ca*(3*gs*Km+Rd-Vcmax)-2*(2*cp*Vcmax+Km*(Rd+Vcmax)))/(gs*sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))))/(h^2*VPD^2)-(2*(-1+c)*c*exp((-(px/d))^c)*h3*(-(px/d))^c*w^(2*b))/(h2^2*kxmax^2*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(2*c^2*exp((-(px/d))^c)*h3*(-(px/d))^(2*c)*w^(2*b))/(h2^2*kxmax^2*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)-(2*(LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))+(2*c*h*h3*(-(px/d))^c*VPD*w^b)/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))))/(gs*h^2*VPD^2)+(2*c*h3*(-(px/d))^c*w^b*(2*h2*kxmax*px^2*w^(2*b)+c^2*(-(px/d))^c*(pe-px*w^b)*((-exp((-(px/d))^c))*gs*h*(-1+(-(px/d))^c)*VPD*w^b+2*h2*kxmax*(-(px/d))^c*(pe-px*w^b))-c*(-(px/d))^c*w^b*(4*h2*kxmax*px*(-pe+px*w^b)+exp((-(px/d))^c)*gs*h*VPD*(pe+px*w^b))))/(gs*h*h2^2*kxmax^2*VPD*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^3))))
        return(list(res))
      })
    }
    
    # Optimize mu
    muf <- function(mu, wL){
      gswL <- gsBCf(wL, mu, wL)
      gs11 <- as.numeric(ode(y=c(gs=gswL), times=c(wL, 1), func=dgsdwf, parms=NULL, mu=mu, wL=wL, maxsteps=100000)[2, 2])#, method=c("ode45"), maxsteps=100000
      gs12 <- gsBCf(1, mu, wL)
      res <- gs11-gs12
      return(res^2)
    }
    
    mu <- optimize(muf, c(-20, 0), tol=.Machine$double.eps, wL=wL)$minimum
    return(mu)
  })
}

# gs(w)
gswf <- function(w, wL, mu, parms){
  with(as.list(c(w, wL, mu, parms)), {
    h <- l*a*LAI/nZ*p
    h2=l*LAI/nZ*p/1000
    gamma <- 1/((MAP/365/k)/1000)*nZ
    
    # Functions
    # Boundary conditions
    gsBCf <- function(w, mu, wL){
      f1 <- function(gs){
        px <- pxf(w, gs)
        res <- (1/2)*gs*h*VPD*((-gs)*LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))+2*((-1+exp(-(-(px/d))^c))*h3+mu-(1/2)*LAI*((-ca)*gs-gs*Km+Rd-Vcmax+sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))))-(2*c*gs*h*h3*(-(px/d))^c*VPD*w^b)/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))))
        return(res^2)
      }
      
      gsmax <- gsmaxf(w)
      res <- optimize(f1, c(gsmax*0.1, gsmax), tol=.Machine$double.eps)$minimum
      return(res)
    }
    
    # dgs/dw
    dgsdwf <- function(w, y, parms, mu, wL){
      with(as.list(c(y, parms, mu, wL)), {
        px <- pxf(w, gs)
        res <- -((-((2*b*(-1+c)*c*gs*h*h3*pe*(-(px/d))^c*VPD*w^(-1+b))/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2))+(2*b*c^2*gs*h*h3*pe*(-(px/d))^(2*c)*VPD*w^(-1+b))/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)-(2*b*c*h3*pe*(-(px/d))^c)/(exp((-(px/d))^c)*(w*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))))+(2*b*c^2*gs*h*h3*pe*(-(px/d))^(2*c)*VPD*w^(-1+b)*(pe-c*pe+c*px*w^b))/(h2*kxmax*((-px)*w^b+c*(-(px/d))^c*(-pe+px*w^b))^3)+(2*gamma*((-1+exp(-(-(px/d))^c))*h3+mu-(1/2)*LAI*((-ca)*gs-gs*Km+Rd-Vcmax+sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))))/(1-wL)+(gs*(-gamma+k/(gs*h*VPD))*(LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))+(2*c*h*h3*(-(px/d))^c*VPD*w^b)/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))))/(1-wL))/(gs*h^2*VPD^2*((LAI*((2*(ca+Km))/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)^2/((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))^(3/2)+(-3*ca^2*gs-3*gs*Km^2-2*ca*(3*gs*Km+Rd-Vcmax)-2*(2*cp*Vcmax+Km*(Rd+Vcmax)))/(gs*sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))))/(h^2*VPD^2)-(2*(-1+c)*c*exp((-(px/d))^c)*h3*(-(px/d))^c*w^(2*b))/(h2^2*kxmax^2*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(2*c^2*exp((-(px/d))^c)*h3*(-(px/d))^(2*c)*w^(2*b))/(h2^2*kxmax^2*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)-(2*(LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))+(2*c*h*h3*(-(px/d))^c*VPD*w^b)/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))))/(gs*h^2*VPD^2)+(2*c*h3*(-(px/d))^c*w^b*(2*h2*kxmax*px^2*w^(2*b)+c^2*(-(px/d))^c*(pe-px*w^b)*((-exp((-(px/d))^c))*gs*h*(-1+(-(px/d))^c)*VPD*w^b+2*h2*kxmax*(-(px/d))^c*(pe-px*w^b))-c*(-(px/d))^c*w^b*(4*h2*kxmax*px*(-pe+px*w^b)+exp((-(px/d))^c)*gs*h*VPD*(pe+px*w^b))))/(gs*h*h2^2*kxmax^2*VPD*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^3))))
        return(list(res))
      })
    }
    
    # gs(w)
    gswf <- function(w, mu, wL){
      gswL <- gsBCf(wL, mu, wL)
      res <- ifelse(w==wL, gswL, as.numeric(ode(y=c(gs=gswL), times=c(wL, w), func=dgsdwf, parms=NULL, mu=mu, wL=wL)[2, 2], maxsteps=100000))#, method=c("ode45"), maxsteps=100000
      return(res)
    }
    
    res <- gswf(w, mu, wL)
    return(res)
  })
}

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
  fnoc <- function(w)1/Ef(w)*exp(-gamma*(w-wL)/(1-wL)+k*1/(1-wL)*integralrEf(w))*1/(1-wL)
  
  f <- function(w)cPDF*fnoc(w)
  Bf1 <- function(w)Bf(w, gswf(w), ca)
  fB <- Vectorize(function(w)f(w)*Bf1(w))
  
  integralfnoc <- integrate(fnoc, wL, 1, rel.tol=.Machine$double.eps^0.25)
  cPDF <- 1/(integralfnoc$value+1/k)
  averB <- integrate(fB, wL, 1, rel.tol=.Machine$double.eps^0.25)
  return(averB$value)
}
