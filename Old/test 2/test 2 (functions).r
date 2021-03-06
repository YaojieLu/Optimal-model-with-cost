
# B(w, gs)
Bf <- function(w, gs,
               ca,
               Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
               a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
               pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=5.71, d=10.05, h3=10){
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)-((ps-x)*h2*kxf(x))
    res <- optimize(f1, c(-20,0), tol=.Machine$double.eps)$minimum
    return(res)
  }
  # xylem water potential function
  pxf <- function(w, gs){
    ps <- pe*w^(-b)
    pxmin <- pxminf(w)
    f1 <- function(x)((ps-x)*h2*kxf(x)-h*VPD*gs)^2
    res <- optimize(f1, c(pxmin, ps), tol=.Machine$double.eps)$minimum
    return(res)
  }
  # maximum gs function for given w
  gsmaxf <- function(w){
    ps <- pe*w^(-b)
    pxmin <- pxminf(w)
    res <- (ps-pxmin)*h2*kxf(pxmin)/(h*VPD)
    return(res)
  }
  # m(w, gs)
  mf <- function(w, gs){
    px <- pxf(w, gs)
    kx <- kxf(px)
    PLC <- 1-kx/kxmax
    res <- h3*PLC
    return(res)
  }
  # A(gs)
  Af <- function(gs)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-sqrt(((ca-Km)*gs+Rd-Vcmax)^2+4*gs*((ca*gs+Rd)*Km+cp*Vcmax)))
  
  res <- Af(gs)-mf(w, gs)
  return(res)
}
# ESS stomatal behaviour function
ESSf1 <- function(w,
                  ca,
                  Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
                  a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                  pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=5.71, d=10.05, h3=10){
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)-((ps-x)*h2*kxf(x))
    res <- optimize(f1, c(-20,0), tol=.Machine$double.eps)$minimum
    return(res)
  }
  # maximum gs function for given w
  gsmaxf <- function(w){
    ps <- pe*w^(-b)
    pxmin <- pxminf(w)
    gsmax <- (ps-pxmin)*h2*kxf(pxmin)/(h*VPD)
    return(gsmax)
  }
  gsmax <- gsmaxf(w)
  
  f1 <- function(gs)Bf(w, gs, ca)
  res <- optimize(f1, c(0, gsmax), tol=.Machine$double.eps, maximum=T)
  return(res$maximum)
}
# ESS B(w)
ESSBf1 <- function(w, ca)Bf(w, ESSf1(w, ca), ca)
# function for wL and gsL (where optimal stomatal behaviour stops transpiring)
LHSoptf <- function(ca){
  f1 <- function(w)ESSBf1(w, ca)
  wL <- uniroot(f1, c(0.1, 1), tol=.Machine$double.eps)
  gsL <- ESSf1(wL$root, ca)
  res <- c(wL$root, gsL)
  return(res)
}
# function for identifying mu
muf <- function(ca, k, MAP,
                Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
                a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=5.71, d=10.05, h3=10,
                MDP=MAP/365, gamma=1/(MDP/k/1000)*nZ){
  
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)-((ps-x)*h2*kxf(x))
    res <- optimize(f1, c(-20,0), tol=.Machine$double.eps)$minimum
    return(res)
  }
  # xylem water potential function
  pxf <- function(w, gs){
    ps <- pe*w^(-b)
    pxmin <- pxminf(w)
    f1 <- function(x)((ps-x)*h2*kxf(x)-h*VPD*gs)^2
    res <- optimize(f1, c(pxmin, ps), tol=.Machine$double.eps)$minimum
    return(res)
  }
  # function for identifying mu
  muf <- function(){
    f1 <- function(mu){
      f1 <- function(w){
        w <- w
        O <- function(gs){
          integrand <- function(gs){
            gs <- gs
            pxf1 <- Vectorize(pxf)
            px <- pxf1(w, gs)
            res <- (h^2*VPD^2*((LAI*((2*ca)/gs+(2*Km)/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)^2/(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))^(3/2)+(-3*ca^2*gs-3*gs*Km^2-2*ca*(3*gs*Km+Rd-Vcmax)-2*(2*cp*Vcmax+Km*(Rd+Vcmax)))/(gs*sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax))))))/(h^2*VPD^2)-(2*(-1+c)*c*exp((-(px/d))^c)*h3*(-(px/d))^c*w^(2*b))/(h2^2*kxmax^2*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(2*c^2*exp((-(px/d))^c)*h3*(-(px/d))^(2*c)*w^(2*b))/(h2^2*kxmax^2*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(4*c*h3*(-(px/d))^c*w^b)/(gs*h*h2*kxmax*VPD*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))+(2*(LAI*(-(ca/gs)-Km/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/(gs*sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))-(2*c*h*h3*(-(px/d))^c*VPD*w^b)/(gs*h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))))/(h^2*VPD^2)))/(gs*(-((2*gamma*((-1+exp(-(-(px/d))^c))*h3+mu-(1/2)*LAI*((-ca)*gs-gs*Km+Rd-Vcmax+sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax))))))/gs^2)+(2*b*(-1+c)*c*h*h3*pe*(-(px/d))^c*VPD*w^(-1+b))/(gs*h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)-(4*b*c^2*h*h3*pe*(-(px/d))^(2*c)*VPD*w^(-1+b))/(gs*h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(2*b*c*h3*pe*(-(px/d))^c)/(exp((-(px/d))^c)*(gs^2*w*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))))+(-gamma+k/(gs*h*VPD))*(-((LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))/gs)-(2*c*h*h3*(-(px/d))^c*VPD*w^b)/(gs*h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))))))
            return(res)
          }
          res <- integrate(integrand, gsL, gs, rel.tol=.Machine$double.eps^0.5)$value-(w-wLopt)
          return(res^2)
        }
        u <- optimize(O, c(gsL, ESSf1(w, ca)), tol=.Machine$double.eps)
        return(u$minimum)
      }
      
      gs1 <- f1(1)
      px1 <- pxf(1, gs1)
      test <- (1/2)*(-2*(-1+exp(-(-(px1/d))^c))*h3-gs1^2*LAI*(-(ca/gs1)-Km/gs1+(ca^2*gs1+gs1*Km^2+Km*Rd+ca*(2*gs1*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/(gs1*sqrt(ca^2*gs1^2+gs1^2*Km^2+(Rd-Vcmax)^2+2*ca*gs1*(gs1*Km+Rd-Vcmax)+2*gs1*(2*cp*Vcmax+Km*(Rd+Vcmax)))))+LAI*((-ca)*gs1-gs1*Km+Rd-Vcmax+sqrt(ca^2*gs1^2+gs1^2*Km^2+(Rd-Vcmax)^2+2*ca*gs1*(gs1*Km+Rd-Vcmax)+2*gs1*(2*cp*Vcmax+Km*(Rd+Vcmax))))+(2*c*gs1*h*h3*(-(px1/d))^c*VPD)/(h2*kxmax*px1+c*h2*kxmax*pe*(-(px1/d))^c-c*h2*kxmax*px1*(-(px1/d))^c))
      res <- (test-mu)^2
      return(res)
    }
    
    resf <- function(x)optimize(f1, c(x, 1e9), tol=.Machine$double.eps)
    estint <- function(p1, p2, p3){
      est <- try(resf(p1), silent=TRUE)
      i <- 0
      if(inherits(est, "try-error")){
        while(inherits(est, "try-error")){
          i <- i + 1
          est <- try(resf(p1+p2*i/(10^p3)), silent=TRUE)
          message(p1+p2*i/(10^p3))
        }
      }
      return(p1+p2*i/(10^p3))
    }
    intlowerm1 <- estint(-20, 1, -1)
    intlower0 <- estint(intlowerm1-1e1, 1, 0)
    intlower1 <- estint(intlower0-1e0, 1, 1)
    intlower2 <- estint(intlower1-1e-1, 1, 2)
    intlower3 <- estint(intlower2-1e-2, 1, 3)
    res <- resf(intlower3)
    return(res)
  }
  res <- muf()
  return(res$minimum)
}
# disf is used to identify whether the integrand function has any discontinuity between 0 and ESS(w). It returns the lowest singularity 
####################################################################
###### Important!!! may fail with different parameter values  ######
####################################################################
disf <- function(w, ca, k, MAP,
                 Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
                 a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                 pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=5.71, d=10.05, h3=10,
                 MDP=MAP/365, gamma=1/(MDP/k/1000)*nZ){
  
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)-((ps-x)*h2*kxf(x))
    res <- optimize(f1, c(-20,0), tol=.Machine$double.eps)$minimum
    return(res)
  }
  # xylem water potential function
  pxf <- function(w, gs){
    ps <- pe*w^(-b)
    pxmin <- pxminf(w)
    f1 <- function(x)((ps-x)*h2*kxf(x)-h*VPD*gs)^2
    res <- optimize(f1, c(pxmin, ps), tol=.Machine$double.eps)$minimum
    return(res)
  }
  
  f1 <- function(gs){
    w <- w
    px <- pxf(w, gs)
    res <- -((2*gamma*((-1+exp(-(-(px/d))^c))*h3+mu-(1/2)*LAI*((-ca)*gs-gs*Km+Rd-Vcmax+sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax))))))/gs^2)+(2*b*(-1+c)*c*h*h3*pe*(-(px/d))^c*VPD*w^(-1+b))/(gs*h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)-(4*b*c^2*h*h3*pe*(-(px/d))^(2*c)*VPD*w^(-1+b))/(gs*h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(2*b*c*h3*pe*(-(px/d))^c)/(exp((-(px/d))^c)*(gs^2*w*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))))+(-gamma+k/(gs*h*VPD))*(-((LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))/gs)-(2*c*h*h3*(-(px/d))^c*VPD*w^b)/(gs*h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))))
    return(res)
  }
  
  gsESS <- ESSf1(w, ca)
  x1 <- optimize(f1, c(0, gsESS), tol=.Machine$double.eps)$minimum
  x2 <- optimize(f1, c(0, x1), maximum=T, tol=.Machine$double.eps)$maximum
  res1 <- uniroot(f1, c(1e-10, x2), tol=.Machine$double.eps)
  res2 <- uniroot(f1, c(x2, gsESS), tol=.Machine$double.eps)
  res <- c(res1$root, res2$root)
  return(res)
}
# optimal stomatal behaviour function
optf1 <- function(w, ca, k, MAP,
                  Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
                  a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                  pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=5.71, d=10.05, h3=10,
                  MDP=MAP/365, gamma=1/(MDP/k/1000)*nZ){
  
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)-((ps-x)*h2*kxf(x))
    res <- optimize(f1, c(-20,0), tol=.Machine$double.eps)$minimum
    return(res)
  }
  # xylem water potential function
  pxf <- function(w, gs){
    ps <- pe*w^(-b)
    pxmin <- pxminf(w)
    f1 <- function(x)((ps-x)*h2*kxf(x)-h*VPD*gs)^2
    res <- optimize(f1, c(pxmin, ps), tol=.Machine$double.eps)$minimum
    return(res)
  }
  
  # optimal gs for given w
  w <- w
  O <- function(gs){
    integrand <- function(gs){
      gs <- gs
      pxf1 <- Vectorize(pxf)
      px <- pxf1(w, gs)
      res <- (h^2*VPD^2*((LAI*((2*ca)/gs+(2*Km)/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)^2/(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))^(3/2)+(-3*ca^2*gs-3*gs*Km^2-2*ca*(3*gs*Km+Rd-Vcmax)-2*(2*cp*Vcmax+Km*(Rd+Vcmax)))/(gs*sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax))))))/(h^2*VPD^2)-(2*(-1+c)*c*exp((-(px/d))^c)*h3*(-(px/d))^c*w^(2*b))/(h2^2*kxmax^2*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(2*c^2*exp((-(px/d))^c)*h3*(-(px/d))^(2*c)*w^(2*b))/(h2^2*kxmax^2*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(4*c*h3*(-(px/d))^c*w^b)/(gs*h*h2*kxmax*VPD*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))+(2*(LAI*(-(ca/gs)-Km/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/(gs*sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))-(2*c*h*h3*(-(px/d))^c*VPD*w^b)/(gs*h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))))/(h^2*VPD^2)))/(gs*(-((2*gamma*((-1+exp(-(-(px/d))^c))*h3+mu-(1/2)*LAI*((-ca)*gs-gs*Km+Rd-Vcmax+sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax))))))/gs^2)+(2*b*(-1+c)*c*h*h3*pe*(-(px/d))^c*VPD*w^(-1+b))/(gs*h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)-(4*b*c^2*h*h3*pe*(-(px/d))^(2*c)*VPD*w^(-1+b))/(gs*h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(2*b*c*h3*pe*(-(px/d))^c)/(exp((-(px/d))^c)*(gs^2*w*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))))+(-gamma+k/(gs*h*VPD))*(-((LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))/gs)-(2*c*h*h3*(-(px/d))^c*VPD*w^b)/(gs*h2*kxmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))))))
      return(res)
    }
    res <- integrate(integrand, gsL, gs, rel.tol=.Machine$double.eps^0.5)$value-(w-wLopt)
    return(res^2)
  }
  Bf1 <- function(gs)Bf(w, gs, ca)
  x1 <- uniroot(Bf1, c(0, 1), tol=.Machine$double.eps)
  res1 <- try(optimize(O, c(x1$root, ESSf1(w, ca)), tol=.Machine$double.eps), silent=T)
  #browser()
  if(is.numeric(try(res1$minimum, silent=T))){
    res <- res1$minimum
  }else{
    res2 <- optimize(O, disf(w, ca, k, MAP), tol=.Machine$double.eps)
    res <- res2$minimum
  }
  return(res)
}
# optimal B(w)
optBf1 <- function(w, ca, k, MAP)Bf(w, optf1(w, ca, k, MAP), ca)
# optimal WUE(w)
optWUEf1 <- function(w, ca, k, MAP)optBf1(w, ca, k, MAP)/optf1(w, ca, k, MAP)
