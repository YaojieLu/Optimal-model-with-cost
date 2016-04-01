
source("test (functions).r")

ca1 <- 400
k1 <- 0.025
MAP1 <- 365

w0opt <- w0optf(ca1)
mu <- muf(ca1, k1, MAP1)

# gsmax
gsmaxf1 <- function(w,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                    pe=-1.58*10^-3, b=4.38, h2=h/1000, xkmax=5, c=5.71, d=10.05, h3=10){
  
  # xylem conductance function
  xkf <- function(x)xkmax*exp(-(-x/d)^c)
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)-(ps-x)*h2*xkf(x)
    res <- optimize(f1, c(-20,0), tol=.Machine$double.eps)$minimum
    return(res)
  }
  ps <- pe*w^(-b)
  pxmin <- pxminf(w)
  res <- (ps-pxmin)*h2*xkf(pxmin)/(h*VPD)
  return(res)
}
gsmaxf2 <- Vectorize(gsmaxf1)

# disf is used to identify whether the integrand function has any discontinuity between 0 and ESS(w). It returns the lowest singularity 
disf <- function(w, ca, k, MAP,
                 Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
                 a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                 pe=-1.58*10^-3, b=4.38, h2=h/1000, xkmax=5, c=5.71, d=10.05, h3=10,
                 MDP=MAP/365, gamma=1/(MDP/k/1000)*nZ){
  
  # xylem conductance function
  xkf <- function(x)xkmax*exp(-(-x/d)^c)
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)-((ps-x)*h2*xkf(x))
    res <- optimize(f1, c(-20,0), tol=.Machine$double.eps)$minimum
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
  
  f1 <- function(gs){
    w <- w
    px <- pxf(w, gs)
    res <- -((2*gamma*((-1+exp(-(-(px/d))^c))*h3+mu-(1/2)*LAI*((-ca)*gs-gs*Km+Rd-Vcmax+sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax))))))/gs^2)+(2*b*(-1+c)*c*h*h3*pe*(-(px/d))^c*VPD*w^(-1+b))/(gs*h2*xkmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)-(4*b*c^2*h*h3*pe*(-(px/d))^(2*c)*VPD*w^(-1+b))/(gs*h2*xkmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(2*b*c*h3*pe*(-(px/d))^c)/(exp((-(px/d))^c)*(gs^2*w*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))))+(-gamma+k/(gs*h*VPD))*(-((LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))/gs)-(2*c*h*h3*(-(px/d))^c*VPD*w^b)/(gs*h2*xkmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))))
    return(res)
  }
  
  x0 <- ESSf1(w, ca)
  x1 <- optimize(f1, c(0, x0), tol=.Machine$double.eps)$minimum
  x2 <- optimize(f1, c(0, x1), tol=.Machine$double.eps, maximum=T)$maximum
  x3 <- try(uniroot(f1, c(x2/1e5, x2), tol=.Machine$double.eps)$root, silent=T)
  res <- ifelse(x3<x0, x3, x0)
  return(res)
}
disf(0.18, ca1, k1, MAP1)
ESSf1(0.18, ca1)
