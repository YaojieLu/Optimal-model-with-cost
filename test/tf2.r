
# Maximum averB
averBf <- function(wL, parms){
  with(as.list(c(wL, parms)), {
    h <- l*a*LAI/nZ*p
    h2 <- l*LAI/nZ*p/1000
    gamma <- 1/((MAP/365/k)/1000)*nZ
    
    # Functions
    # Boundary conditions
    gsBCf <- function(w, mu, wL){
      f1 <- function(gs){
        px <- pxf(w, gs)
        res <- (-gs)*LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))+2*((-1+exp(-(-(px/d))^c))*h3+mu-(1/2)*LAI*((-ca)*gs-gs*Km+Rd-Vcmax+sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))))-(2*c*gs*h*h3*(-(px/d))^c*VPD*w^b)/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))
        return(res^2)
      }
      
      gsmax <- gsmaxf(w)
      res <- optimize(f1, c(0, gsmax), tol=.Machine$double.eps)$minimum
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
      #integralrEf <- Vectorize(function(w)integrate(rEf, wL, w, rel.tol=.Machine$double.eps^0.25)$value)
      integralrEf <- Vectorize(function(w)integral(rEf, wL, w, reltol=1e-12))
      fnoc <- function(w)1/Ef(w)*exp(-gamma*(w-wL)/(1-wL)+k*1/(1-wL)*integralrEf(w))*1/(1-wL)
      
      #res <- integrate(fnoc, wL, 1, rel.tol=.Machine$double.eps^0.25)
      res <- integral(fnoc, wL, 1, reltol=1e-12)
      return(res)
    }
    
    # PDF of w
    PDFf <- function(w, wL, cPDF){
      Ef <- function(w)h*VPD*gswf1(w)
      rEf <- function(w)1/Ef(w)
      #integralrEf <- Vectorize(function(w)integrate(rEf, wL, w, rel.tol=.Machine$double.eps^0.25)$value)
      integralrEf <- Vectorize(function(w)integral(rEf, wL, w, reltol=1e-12))
      
      res <- cPDF/Ef(w)*exp(-gamma*(w-wL)/(1-wL)+k*1/(1-wL)*integralrEf(w))*1/(1-wL)
      return(res)
    }
    
    # Average B
    averBf <- function(wL, cPDF){
      Bf1 <- Vectorize(function(w)Bf(w, gswf1(w), ca))
      f1 <- function(w)Bf1(w)*PDFf(w, wL, cPDF)
      
      #res <- integrate(f1, wL, 1, rel.tol=.Machine$double.eps^0.25)
      res <- integral(f1, wL, 1, reltol=1e-12)
      return(res)
    }
    browser()
    mu <- optimize(muf, c(-30, 0), tol=.Machine$double.eps, wL=wL)
    gswf1 <- Vectorize(function(w)gswf(w, mu$minimum, wL))
    integralfnoc <- integralfnocf(wL)
    #cPDF <- 1/(integralfnoc$value+1/k)
    cPDF <- 1/(integralfnoc+1/k)
    averB <- averBf(wL, cPDF)
    browser()
    #return(averB$value)
    return(averB)
  })
}
