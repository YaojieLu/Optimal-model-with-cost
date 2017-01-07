
options(digits=20)
source("Functions 1.r")

ca <- 400
k <- 0.05
MAP <- 1000
h3 <- 10

parms <- c(LAI=1,
           Vcmax=50, cp=30, Km=703, Rd=1,
           a=1.6, nZ=0.5, p=43200, l=1.8e-5, VPD=0.02,
           pe=-1.58*10^-3, b=4.38, kxmax=5, c=2.64, d=3.54)

gsBCf <- function(w, parms){
  with(as.list(c(parms)), {
    h <- l*a*LAI/nZ*p
    h2 <- l*LAI/nZ*p/1000
    gamma <- 1/((MAP/365/k)/1000)*nZ
    gsmax <- gsmaxf(w)
    
    f1 <- Vectorize(function(gs){
      px <- pxf(w, gs)
      res <- (1/2)*(2*(-1+exp(-(-(px/d))^c))*h3+LAI*(gs*(ca+Km)-Rd+Vcmax-sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))-2*gs*((1/2)*LAI*(ca+Km+((-ca^2)*gs-gs*Km^2-Km*Rd-2*cp*Vcmax-Km*Vcmax+ca*(-2*gs*Km-Rd+Vcmax))/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))+(c*h*h3*(-(px/d))^c*VPD*w^b)/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))))
      #res <- (1/2)*(2*(-1+exp(-(-(px/d))^c))*h3+2*mu+LAI*(gs*(ca+Km)-Rd+Vcmax-sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))-2*gs*((1/2)*LAI*(ca+Km+((-ca^2)*gs-gs*Km^2-Km*Rd-2*cp*Vcmax-Km*Vcmax+ca*(-2*gs*Km-Rd+Vcmax))/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))+(c*h*h3*(-(px/d))^c*VPD*w^b)/(h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))))
      return(res)
    })
    
    browser()
    curve(f1, 0, gsmax*0.8)
    res <- uniroot(f1, c(0, gsmax), tol=.Machine$double.eps)
    return(res)
  })
}

gsBCf(0.187923801472611, parms)