
# function for w0opt (where optimal stomatal behaviour stops transpiring)
w0optf <- function(ca,
                   Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
                   a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                   pe=-1.58*10^-3, b=4.38, h2=h/1000, xkmax=5, c=5.71, d=10.05, h3=10){
  
  # photosynthesis rate function
  Af <- function(gs){LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-sqrt(((ca-Km)*gs+Rd-Vcmax)^2+4*gs*((ca*gs+Rd)*Km+cp*Vcmax)))}
  # xylem conductance function
  xkf <- function(x)xkmax*exp(-(-x/d)^c)
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)-(ps-x)*h2*xkf(x)
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
  # PLC cost function
  mf <- function(w, gs){
    px <- pxf(w, gs)
    xk <- xkf(px)
    PLC <- 1-xk/xkmax
    res <- h3*PLC
    return(res)
  }
  # maximum gs function for given w
  gsmaxf <- function(w){
    ps <- pe*w^(-b)
    pxmin <- pxminf(w)
    gsmax <- (ps-pxmin)*h2*xkf(pxmin)/(h*VPD)
    return(gsmax)
  }
  # net photosynthesis rate function
  Bf <- function(w, gs)Af(gs)-mf(w, gs)
  # maximum B function for given w
  Bmaxf <- function(w){
    f1 <- function(gs)-Bf(w, gs)
    gsmax <- gsmaxf(w)
    res <- optimize(f1, c(0, gsmax), tol=.Machine$double.eps)$objective
    return(res^2)
  }
  res <- optimize(Bmaxf, c(0.1, 1), tol=.Machine$double.eps)$minimum
  return(res)
}
# ESS stomatal behaviour function
ESSf1 <- function(w,
                  ca,
                  Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
                  a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                  pe=-1.58*10^-3, b=4.38, h2=h/1000, xkmax=5, c=5.71, d=10.05, h3=10){
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
  # maximum gs function for given w
  gsmaxf <- function(w){
    ps <- pe*w^(-b)
    pxmin <- pxminf(w)
    gsmax <- (ps-pxmin)*h2*xkf(pxmin)/(h*VPD)
    return(gsmax)
  }
  gsmax <- gsmaxf(w)
  
  f1 <- function(gs){
    px <- pxf(w, gs)
    test <- ((-c)*h*h2*h3*(ca+Km)*xkmax*LAI*px*(-(px/d))^c*(ca*(Rd-Vcmax)+2*cp*Vcmax+Km*(Rd+Vcmax))*VPD*w^(2*b)-c^2*h*h3*(-(px/d))^(2*c)*(ca*(Rd-Vcmax)+2*cp*Vcmax+Km*(Rd+Vcmax))*VPD*w^b*(h*h3*VPD*w^b+ca*h2*xkmax*LAI*(pe-px*w^b)+h2*Km*xkmax*LAI*(pe-px*w^b))-sqrt(c*h*h3*(cp+Km)*(-(px/d))^c*Vcmax*(Km*Rd+ca*(Rd-Vcmax)+cp*Vcmax)*VPD*w^b*(h2*(ca+Km)*xkmax*LAI*px*w^b+c*(-(px/d))^c*(h*h3*VPD*w^b+ca*h2*xkmax*LAI*(pe-px*w^b)+h2*Km*xkmax*LAI*(pe-px*w^b)))*(h2*(ca+Km)*xkmax*LAI*px*w^b+c*(-(px/d))^c*(2*h*h3*VPD*w^b+ca*h2*xkmax*LAI*(pe-px*w^b)+h2*Km*xkmax*LAI*(pe-px*w^b)))^2))/((-(px/d))^c*w^b)/(c*h*h3*(ca+Km)^2*VPD*(h2*(ca+Km)*xkmax*LAI*px*w^b+c*(-(px/d))^c*(h*h3*VPD*w^b+ca*h2*xkmax*LAI*(pe-px*w^b)+h2*Km*xkmax*LAI*(pe-px*w^b))))
    res <- (gs-test)^2
    return(res)
  }
  res <- optimize(f1, c(0, gsmax), tol=.Machine$double.eps)$minimum
  return(res)
}
ESSf2 <- Vectorize(function(w)ESSf1(w, ca=ca1))
# ESS B(w)
ESSBf1 <- function(w, ca, k, MAP,
                   Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
                   a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                   pe=-1.58*10^-3, b=4.38, h2=h/1000, xkmax=5, c=5.71, d=10.05, h3=10){
  
  # photosynthesis rate function
  Af <- function(gs){LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-sqrt(((ca-Km)*gs+Rd-Vcmax)^2+4*gs*((ca*gs+Rd)*Km+cp*Vcmax)))}
  # xylem conductance function
  xkf <- function(x)xkmax*exp(-(-x/d)^c)
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)-(ps-x)*h2*xkf(x)
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
  # PLC cost function
  mf <- function(w, gs){
    px <- pxf(w, gs)
    xk <- xkf(px)
    PLC <- 1-xk/xkmax
    res <- h3*PLC
    return(res)
  }
  # net photosynthesis rate function
  Bf <- function(w, gs)Af(gs)-mf(w, gs)
  f1 <- function(w, ca)Bf(w, ESSf1(w, ca))
  res <- f1(w, ca)
  return(res)
}

# function for w0ESS (where ESS stomatal behaviour stops transpiring)
w0ESSf <- function(w0, ca,
                   gs=0,
                   Vcmax=50, cp=30, Km=703, Rd=1,
                   a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                   pe=-1.58*10^-3, b=4.38, h2=h/1000, kmax=5, c=5.71, d=10.05, h3=10){
  res <- ((-c)*h*h2*h3*(ca+Km)*kmax*LAI*pe*(ca*(Rd-Vcmax)+2*cp*Vcmax+Km*(Rd+Vcmax))*VPD*w0^b*(-(pe/(w0^b*d)))^c-c^2*h^2*h3^2*(ca*(Rd-Vcmax)+2*cp*Vcmax+Km*(Rd+Vcmax))*VPD^2*w0^(2*b)*(-(pe/(w0^b*d)))^(2*c)-sqrt(c*h*h3*(cp+Km)*Vcmax*(Km*Rd+ca*(Rd-Vcmax)+cp*Vcmax)*VPD*w0^b*(-(pe/(w0^b*d)))^c*(h2*(ca+Km)*kmax*LAI*pe+c*h*h3*VPD*w0^b*(-(pe/(w0^b*d)))^c)*(h2*(ca+Km)*kmax*LAI*pe+2*c*h*h3*VPD*w0^b*(-(pe/(w0^b*d)))^c)^2))/(w0^b*(-(pe/(w0^b*d)))^c)/(c*h*h3*(ca+Km)^2*VPD*(h2*(ca+Km)*kmax*LAI*pe+c*h*h3*VPD*w0^b*(-(pe/(w0^b*d)))^c))
  return(abs(res))
}
# function for identifying mu
muf <- function(ca, k, MAP,
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
  # maximum gs function for given w
  gsmaxf <- function(w){
    ps <- pe*w^(-b)
    pxmin <- pxminf(w)
    gsmax <- (ps-pxmin)*h2*xkf(pxmin)/(h*VPD)
    return(gsmax)
  }
  # function for identifying mu
  muf <- function(){
    f1 <- function(mu){
      f1 <- function(w){
        w <- w
        O <- function(gs){
          px <- pxf(w, gs)
          integrand <- function(gs)(h^2*VPD^2*((LAI*((2*ca)/gs+(2*Km)/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)^2/(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))^(3/2)+(-3*ca^2*gs-3*gs*Km^2-2*ca*(3*gs*Km+Rd-Vcmax)-2*(2*cp*Vcmax+Km*(Rd+Vcmax)))/(gs*sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax))))))/(h^2*VPD^2)-(2*(-1+c)*c*exp((-(px/d))^c)*h3*(-(px/d))^c*w^(2*b))/(h2^2*xkmax^2*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(2*c^2*exp((-(px/d))^c)*h3*(-(px/d))^(2*c)*w^(2*b))/(h2^2*xkmax^2*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(4*c*h3*(-(px/d))^c*w^b)/(gs*h*h2*xkmax*VPD*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))+(2*(LAI*(-(ca/gs)-Km/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/(gs*sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))-(2*c*h*h3*(-(px/d))^c*VPD*w^b)/(gs*h2*xkmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))))/(h^2*VPD^2)))/(gs*(-((2*gamma*((-1+exp(-(-(px/d))^c))*h3+mu-(1/2)*LAI*((-ca)*gs-gs*Km+Rd-Vcmax+sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax))))))/gs^2)+(2*b*(-1+c)*c*h*h3*pe*(-(px/d))^c*VPD*w^(-1+b))/(gs*h2*xkmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)-(4*b*c^2*h*h3*pe*(-(px/d))^(2*c)*VPD*w^(-1+b))/(gs*h2*xkmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(2*b*c*h3*pe*(-(px/d))^c)/(exp((-(px/d))^c)*(gs^2*w*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))))+(-gamma+k/(gs*h*VPD))*(-((LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))/gs)-(2*c*h*h3*(-(px/d))^c*VPD*w^b)/(gs*h2*xkmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))))))
          res <- integrate(integrand, 0, gs, rel.tol=.Machine$double.eps^0.5)$value-(w-w0opt)
          return(res^2)
        }
        u <- optimize(O, c(0, ESSf1(w, ca)), tol=.Machine$double.eps)$minimum
        return(u)
      }
      
      gs1 <- f1(1)
      px1 <- pxf(1, gs1)
      test <- (1/2)*(-2*(-1+exp(-(-(px1/d))^c))*h3-gs1^2*LAI*(-(ca/gs1)-Km/gs1+(ca^2*gs1+gs1*Km^2+Km*Rd+ca*(2*gs1*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/(gs1*sqrt(ca^2*gs1^2+gs1^2*Km^2+(Rd-Vcmax)^2+2*ca*gs1*(gs1*Km+Rd-Vcmax)+2*gs1*(2*cp*Vcmax+Km*(Rd+Vcmax)))))+LAI*((-ca)*gs1-gs1*Km+Rd-Vcmax+sqrt(ca^2*gs1^2+gs1^2*Km^2+(Rd-Vcmax)^2+2*ca*gs1*(gs1*Km+Rd-Vcmax)+2*gs1*(2*cp*Vcmax+Km*(Rd+Vcmax))))+(2*c*gs1*h*h3*(-(px1/d))^c*VPD)/(h2*xkmax*px1+c*h2*xkmax*pe*(-(px1/d))^c-c*h2*xkmax*px1*(-(px1/d))^c))
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
  res <- muf()$minimum
  return(res)
}
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
  #browser()
  res <- ifelse(x3<x2, x3, x2)
  #res <- ifelse(x3<x0, x3, x0)
  return(res)
}
# optimal stomatal behaviour function
optf1 <- function(w, ca, k, MAP,
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
  
  # optimal gs for given w
  w <- w
  O <- function(gs){
    px <- pxf(w, gs)
    integrand <- function(gs)(h^2*VPD^2*((LAI*((2*ca)/gs+(2*Km)/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)^2/(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))^(3/2)+(-3*ca^2*gs-3*gs*Km^2-2*ca*(3*gs*Km+Rd-Vcmax)-2*(2*cp*Vcmax+Km*(Rd+Vcmax)))/(gs*sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax))))))/(h^2*VPD^2)-(2*(-1+c)*c*exp((-(px/d))^c)*h3*(-(px/d))^c*w^(2*b))/(h2^2*xkmax^2*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(2*c^2*exp((-(px/d))^c)*h3*(-(px/d))^(2*c)*w^(2*b))/(h2^2*xkmax^2*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(4*c*h3*(-(px/d))^c*w^b)/(gs*h*h2*xkmax*VPD*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))+(2*(LAI*(-(ca/gs)-Km/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/(gs*sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))-(2*c*h*h3*(-(px/d))^c*VPD*w^b)/(gs*h2*xkmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))))/(h^2*VPD^2)))/(gs*(-((2*gamma*((-1+exp(-(-(px/d))^c))*h3+mu-(1/2)*LAI*((-ca)*gs-gs*Km+Rd-Vcmax+sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax))))))/gs^2)+(2*b*(-1+c)*c*h*h3*pe*(-(px/d))^c*VPD*w^(-1+b))/(gs*h2*xkmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)-(4*b*c^2*h*h3*pe*(-(px/d))^(2*c)*VPD*w^(-1+b))/(gs*h2*xkmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(2*b*c*h3*pe*(-(px/d))^c)/(exp((-(px/d))^c)*(gs^2*w*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))))+(-gamma+k/(gs*h*VPD))*(-((LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))/gs)-(2*c*h*h3*(-(px/d))^c*VPD*w^b)/(gs*h2*xkmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))))))
    res <- integrate(integrand, 0, gs, rel.tol=.Machine$double.eps^0.5)$value-(w-w0opt)
    return(res^2)
  }
  res1 <- try(optimize(O, c(0, ESSf1(w, ca)), tol=.Machine$double.eps)$minimum, silent=T)
  res <- ifelse(is.numeric(res1), res1, optimize(O, c(0, disf(w, ca, k, MAP)), tol=.Machine$double.eps)$minimum)
  #browser()
  return(res)
}
# optimal B(w)
optBf1 <- function(w, ca, k, MAP,
                   Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
                   a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                   pe=-1.58*10^-3, b=4.38, h2=h/1000, xkmax=5, c=5.71, d=10.05, h3=10){
  
  # photosynthesis rate function
  Af <- function(gs){LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-sqrt(((ca-Km)*gs+Rd-Vcmax)^2+4*gs*((ca*gs+Rd)*Km+cp*Vcmax)))}
  # xylem conductance function
  xkf <- function(x)xkmax*exp(-(-x/d)^c)
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)-(ps-x)*h2*xkf(x)
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
  # PLC cost function
  mf <- function(w, gs){
    px <- pxf(w, gs)
    xk <- xkf(px)
    PLC <- 1-xk/xkmax
    res <- h3*PLC
    return(res)
  }
  # net photosynthesis rate function
  Bf <- function(w, gs)Af(gs)-mf(w, gs)
  f1 <- function(w)Bf(w, optf1(w, ca, k, MAP))
  res <- f1(w)
  return(res)
}
