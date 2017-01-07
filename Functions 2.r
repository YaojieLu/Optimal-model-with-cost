
# averB given wL
averBf <- function(wL, parms){
  with(as.list(c(wL, parms)), {
    h <- l*a*LAI/nZ*p
    h2 <- l*LAI/nZ*p/1000
    gamma <- 1/((MAP/365/k)/1000)*nZ
    
    # Functions
    # Boundary conditions at wL
    gsBCwLf <- function(wL){
      f1 <- function(gs){
        w <- wL
        px <- pxf(w, gs)
        res <- (1/2)*(2*(-1+exp(-(-(px/d))^c))*h3+LAI*(gs*(ca+Km)-Rd+Vcmax-sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))-2*gs*((1/2)*LAI*(ca+Km+((-ca^2)*gs-gs*Km^2-Km*Rd-2*cp*Vcmax-Km*Vcmax+ca*(-2*gs*Km-Rd+Vcmax))/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))+(c*h*h3*(-(px/d))^c*VPD*w^b)/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))))
        return(res)
      }
      
      gsmax <- gsmaxf(wL)
      res <- uniroot(f1, c(0, gsmax), tol=.Machine$double.eps)
      return(res)
    }
    
    # Boundary conditions at 1
    gsBC1f <- function(mu, wL){
      f1 <- function(gs){
        w <- 1
        px <- pxf(w, gs)
        res <- (1/2)*(2*(-1+exp(-(-(px/d))^c))*h3+2*mu+LAI*(gs*(ca+Km)-Rd+Vcmax-sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))-2*gs*((1/2)*LAI*(ca+Km+((-ca^2)*gs-gs*Km^2-Km*Rd-2*cp*Vcmax-Km*Vcmax+ca*(-2*gs*Km-Rd+Vcmax))/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))+(c*h*h3*(-(px/d))^c*VPD*w^b)/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))))
        return(res)
      }
      
      gsmax <- gsmaxf(1)
      res <- uniroot(f1, c(0, gsmax), tol=.Machine$double.eps)
      return(res)
    }
    
    # dgs/dw
    dgsdwf <- function(w, y, parms, mu, wL){
      with(as.list(c(y, parms, mu, wL)), {
        px <- pxf(w, gs)
        res <- -((2*gamma*((-1+exp(-(-(px/d))^c))*h3+mu-(1/2)*LAI*((-ca)*gs-gs*Km+Rd-Vcmax+sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))))-(2*b*(-1+c)*c*gs*h*h3*pe*(-(px/d))^c*VPD*w^(-1+b))/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(2*b*c^2*gs*h*h3*pe*(-(px/d))^(2*c)*VPD*w^(-1+b))/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)-(2*b*c*h3*pe*(-(px/d))^c)/(exp((-(px/d))^c)*(w*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))))+(2*b*c^2*gs*h*h3*pe*(-(px/d))^(2*c)*VPD*w^(-1+b)*(pe-c*pe+c*px*w^b))/(h2*kxmax*((-px)*w^b+c*(-(px/d))^c*(-pe+px*w^b))^3)+gs*(-gamma+k/(gs*h*VPD))*(LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))+(2*c*h*h3*(-(px/d))^c*VPD*w^b)/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))))/(gs*h^2*VPD^2*((LAI*((2*(ca+Km))/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)^2/((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))^(3/2)+(-3*ca^2*gs-3*gs*Km^2-2*ca*(3*gs*Km+Rd-Vcmax)-2*(2*cp*Vcmax+Km*(Rd+Vcmax)))/(gs*sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))))/(h^2*VPD^2)-(2*(-1+c)*c*exp((-(px/d))^c)*h3*(-(px/d))^c*w^(2*b))/(h2^2*kxmax^2*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(2*c^2*exp((-(px/d))^c)*h3*(-(px/d))^(2*c)*w^(2*b))/(h2^2*kxmax^2*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)-(2*(LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))+(2*c*h*h3*(-(px/d))^c*VPD*w^b)/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))))/(gs*h^2*VPD^2)+(2*c*h3*(-(px/d))^c*w^b*(2*h2*kxmax*px^2*w^(2*b)+c^2*(-(px/d))^c*(pe-px*w^b)*((-exp((-(px/d))^c))*gs*h*(-1+(-(px/d))^c)*VPD*w^b+2*h2*kxmax*(-(px/d))^c*(pe-px*w^b))-c*(-(px/d))^c*w^b*(4*h2*kxmax*px*(-pe+px*w^b)+exp((-(px/d))^c)*gs*h*VPD*(pe+px*w^b))))/(gs*h*h2^2*kxmax^2*VPD*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^3))))
        return(list(res))
      })
    }
    
    # Optimize mu
    muf <- function(mu, wL){
      gs11 <- as.numeric(ode(y=c(gs=gswL), times=c(wL, 1), func=dgsdwf, parms=NULL, mu=mu, wL=wL, maxsteps=100000)[2, 2])
      gs12 <- gsBC1f(mu, wL)$root
      res <- gs11-gs12
      return(res)
    }
    
    # gs(w)
    gswf <- function(w, mu, wL){
      res <- as.numeric(ode(y=c(gs=gswL), times=c(wL, w), func=dgsdwf, parms=NULL, mu=mu, wL=wL, maxsteps=100000)[2, 2])
      return(res)
    }
    
    # integralfnoc of PDF
    integralfnocf <- function(wL){
      Ef <- function(w)h*VPD*gswf1(w)
      rEf <- function(w)1/Ef(w)
      integralrEf <- Vectorize(function(w)integrate(rEf, wL, w, rel.tol=.Machine$double.eps^0.3)$value)
      fnoc <- function(w)1/Ef(w)*exp(-gamma*w+k*integralrEf(w))
      
      res <- integrate(fnoc, wL, 1, rel.tol=.Machine$double.eps^0.3)
      return(res)
    }
    
    # PDF of w
    PDFf <- function(w, wL, cPDF){
      Ef <- function(w)h*VPD*gswf1(w)
      rEf <- function(w)1/Ef(w)
      integralrEf <- Vectorize(function(w)integrate(rEf, wL, w, rel.tol=.Machine$double.eps^0.3)$value)
      
      res <- cPDF/Ef(w)*exp(-gamma*w+k*integralrEf(w))
      return(res)
    }
    
    # Average B
    averBf <- function(wL, cPDF){
      Bf1 <- Vectorize(function(w)Bf(w, gswf1(w)))
      f1 <- function(w)Bf1(w)*PDFf(w, wL, cPDF)
      
      res <- integrate(f1, wL, 1, rel.tol=.Machine$double.eps^0.3)
      return(res)
    }

    browser()
    gswL <- gsBCwLf(wL)$root
    df <- Vectorize(function(mu)dgsdwf(wL, y=c(gs=gswL), parms=NULL, mu, wL))
    muf1 <- Vectorize(function(mu)muf(mu, wL))
    mu <- uniroot(muf, c(-20, 0), tol=.Machine$double.eps, wL=wL)
    gswf1 <- Vectorize(function(w)gswf(w, mu$root, wL))
    #Bf1 <- Vectorize(function(w)Bf(w, gswf1(w)))
    #curve(Bf1, wL, 1, xlim=c(0, 1))
    #pxf1 <- Vectorize(function(w)pxf(w, gswf1(w)))
    #curve(pxf1, wL, 1, xlim=c(0, 1))
    integralfnoc <- integralfnocf(wL)
    cPDF <- 1/(integralfnoc$value+1/k*exp(-gamma*wL))
    averB <- averBf(wL, cPDF)
    return(averB$value)
  })
}
