
f1 <- function(gs,
               w=0.02,
               ca=ca1, k=k1, MAP=MAP1,
               mu=mu2,
               pe=-1.58, w0=2/3, b=4.38, kmax=5, c=5.71, d=10.05, h3=10,
               a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
               h2=h/1000,
               Vcmax=50, cp=30, Km=703, Rd=1,
               MDP=MAP/365, gamma=1/(MDP/k/1000)*nZ){
  
  pxf1 <- function(pxw, pxgs){
    ps <- pe*(w0+(1-w0)*pxw)^(-b)
    kf <- function(x)kmax*exp(-(-x/d)^c)
    f1 <- function(x)-((ps-x)*h2*kf(x))
    minimum <- optimize(f1, c(-20,0))$minimum
    f2 <- function(x)((ps-x)*h2*kf(x)-h*VPD*pxgs)^2
    px <- optimize(f2, c(minimum, ps))$minimum
    return(px)
  }
  pxf2 <- Vectorize(pxf1)
  
  px <- pxf2(pxw=w, pxgs=gs)
  integrand <- (gs^3*h*LAI*((2*ca)/gs+(2*Km)/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)^2/(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))^(3/2)+(-3*ca^2*gs-3*gs*Km^2-2*ca*(3*gs*Km+Rd-Vcmax)-2*(2*cp*Vcmax+Km*(Rd+Vcmax)))/(gs*sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))*VPD-(2*(-1+c)*c*exp((-(px/d))^c)*gs^3*h^3*h3*(-(px/d))^c*VPD^3*(w+w0-w*w0)^(2*b))/(h2^2*kmax^2*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)+(2*c^2*exp((-(px/d))^c)*gs^3*h^3*h3*(-(px/d))^(2*c)*VPD^3*(w+w0-w*w0)^(2*b))/(h2^2*kmax^2*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)+(4*c*gs^2*h^2*h3*(-(px/d))^c*VPD^2*(w+w0-w*w0)^b)/(h2*kmax*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b)))+2*gs^2*h*VPD*(gs*LAI*(-(ca/gs)-Km/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/(gs*sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))-(2*c*h*h3*(-(px/d))^c*VPD*(w+w0-w*w0)^b)/(h2*kmax*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b)))))/(gs^2*h*VPD*(-2*gamma*((-1+exp(-(-(px/d))^c))*h3+mu-(1/2)*LAI*((-ca)*gs-gs*Km+Rd-Vcmax+sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))-(2*b*(-1+c)*c*gs*h*h3*pe*(-(px/d))^c*VPD*(-1+w0)*(w+w0-w*w0)^(-1+b))/(h2*kmax*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)+(2*b*c^2*gs*h*h3*pe*(-(px/d))^(2*c)*VPD*(-1+w0)*(w+w0-w*w0)^(-1+b))/(h2*kmax*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)-(2*b*c^2*gs*h*h3*pe*(-(px/d))^(2*c)*VPD*(-1+w0)*(w+w0-w*w0)^b)/(h2*kmax*(w*(-1+w0)-w0)*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)-(2*b*c*h3*pe*(-(px/d))^c*(-1+w0))/(exp((-(px/d))^c)*((w+w0-w*w0)*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))))+gs*(-gamma+k/(gs*h*VPD))*((-LAI)*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax))))-(2*c*h*h3*(-(px/d))^c*VPD*(w+w0-w*w0)^b)/(h2*kmax*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))))))
  return(integrand)
}

curve(f1,0,0.3)
abline(h=0, col="red")
