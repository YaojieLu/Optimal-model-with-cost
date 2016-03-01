
options(digits=20)

ca1 <- 400
k1 <- 0.1
MAP1 <- 3650
gs1 <- 0.1
# Function
# w0 & ESS
w0f <- function(w,
                ca=ca1, Vcmax=50, cp=30, Km=703, Rd=1,
                a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                pe=-1.58*10^-3, b=4.38, h2=h/1000, kmax=5, c=5.71, d=10.05, h3=10){
  res <- ((-c)*h*h2*h3*(ca+Km)*kmax*LAI*pe*(ca*(Rd-Vcmax)+2*cp*Vcmax+Km*(Rd+Vcmax))*VPD*w^b*(-(pe/(w^b*d)))^c-c^2*h^2*h3^2*(ca*(Rd-Vcmax)+2*cp*Vcmax+Km*(Rd+Vcmax))*VPD^2*w^(2*b)*(-(pe/(w^b*d)))^(2*c)-sqrt(c*h*h3*(cp+Km)*Vcmax*(Km*Rd+ca*(Rd-Vcmax)+cp*Vcmax)*VPD*w^b*(-(pe/(w^b*d)))^c*(h2*(ca+Km)*kmax*LAI*pe+c*h*h3*VPD*w^b*(-(pe/(w^b*d)))^c)*(h2*(ca+Km)*kmax*LAI*pe+2*c*h*h3*VPD*w^b*(-(pe/(w^b*d)))^c)^2))/(w^b*(-(pe/(w^b*d)))^c)/(c*h*h3*(ca+Km)^2*VPD*(h2*(ca+Km)*kmax*LAI*pe+c*h*h3*VPD*w^b*(-(pe/(w^b*d)))^c))
  return(res^2)
}

w0 <- optimize(w0f, c(0.5, 0.9), tol=.Machine$double.eps)$minimum

ESS1 <- function(w,
                 ca=400, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
                 a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                 pe=-1.58*10^-3, b=4.38, h2=h/1000, kmax=5, c=5.71, d=10.05, h3=10){
  
  pxf1 <- function(gs){
    ps <- pe*(w0+(1-w0)*w)^(-b)
    kf <- function(x)kmax*exp(-(-x/d)^c)
    f1 <- function(x)-((ps-x)*h2*kf(x))
    minimum <- optimize(f1, c(-20,0), tol = .Machine$double.eps)$minimum
    f2 <- function(x)((ps-x)*h2*kf(x)-h*VPD*gs)^2
    px <- optimize(f2, c(minimum, ps), tol = .Machine$double.eps)$minimum
    return(px)
  }
  pxf2 <- Vectorize(pxf1)
  
  gsmax1 <- function(w){
    ps <- pe*(w0+(1-w0)*w)^(-b)
    kf <- function(x)kmax*exp(-(-x/d)^c)
    f1 <- function(x)-((ps-x)*h2*kf(x))
    pxmin <- optimize(f1, c(-20,0), tol = .Machine$double.eps)$minimum
    gsmax <- (ps-pxmin)*h2*kf(pxmin)/(h*VPD)
    return(gsmax)
  }
  gsmax2 <- Vectorize(gsmax1)
  gsmax <- gsmax2(w)
  
  f <- function(gs){
    px <- pxf2(gs)
    test <- ((-c)*h*h2*h3*(ca+Km)*kmax*LAI*px*(-(px/d))^c*(ca*(Rd-Vcmax)+2*cp*Vcmax+Km*(Rd+Vcmax))*VPD*(w+w0-w*w0)^(2*b)-c^2*h*h3*(-(px/d))^(2*c)*(ca*(Rd-Vcmax)+2*cp*Vcmax+Km*(Rd+Vcmax))*VPD*(w+w0-w*w0)^b*(h*h3*VPD*(w+w0-w*w0)^b+ca*h2*kmax*LAI*(pe-px*(w+w0-w*w0)^b)+h2*Km*kmax*LAI*(pe-px*(w+w0-w*w0)^b))-sqrt(c*h*h3*(cp+Km)*(-(px/d))^c*Vcmax*(Km*Rd+ca*(Rd-Vcmax)+cp*Vcmax)*VPD*(w+w0-w*w0)^b*(h2*(ca+Km)*kmax*LAI*px*(w+w0-w*w0)^b+c*(-(px/d))^c*(h*h3*VPD*(w+w0-w*w0)^b+ca*h2*kmax*LAI*(pe-px*(w+w0-w*w0)^b)+h2*Km*kmax*LAI*(pe-px*(w+w0-w*w0)^b)))*(h2*(ca+Km)*kmax*LAI*px*(w+w0-w*w0)^b+c*(-(px/d))^c*(2*h*h3*VPD*(w+w0-w*w0)^b+ca*h2*kmax*LAI*(pe-px*(w+w0-w*w0)^b)+h2*Km*kmax*LAI*(pe-px*(w+w0-w*w0)^b)))^2))/((-(px/d))^c*(w+w0-w*w0)^b)/(c*h*h3*(ca+Km)^2*VPD*(h2*(ca+Km)*kmax*LAI*px*(w+w0-w*w0)^b+c*(-(px/d))^c*(h*h3*VPD*(w+w0-w*w0)^b+ca*h2*kmax*LAI*(pe-px*(w+w0-w*w0)^b)+h2*Km*kmax*LAI*(pe-px*(w+w0-w*w0)^b))))
    res <- (gs-test)^2
    return(res)
  }
  
  res <- optimize(f, c(0, gsmax), tol = .Machine$double.eps)$minimum
  return(res)
}
ESS2 <- Vectorize(ESS1)
ESS3 <- Vectorize(function(w)ifelse(w<w0, 0, ESS2((w-w0)/(1-w0))))

# gs1
gsf1 <- function(x,
                 ca=ca1, k=k1, MAP=MAP1,
                 a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                 h2=h/1000,
                 pe=-1.58*10^-3, b=4.38, kmax=5, c=5.71, d=10.05, h3=10,
                 Vcmax=50, cp=30, Km=703, Rd=1,
                 MDP=MAP/365, gamma=1/(MDP/k/1000)*nZ){
  
  pxf1 <- function(w, gs){
    ps <- pe*w^(-b)
    kf <- function(x)kmax*exp(-(-x/d)^c)
    f1 <- function(x)-((ps-x)*h2*kf(x))
    minimum <- optimize(f1, c(-20,0), tol=.Machine$double.eps)$minimum
    f2 <- function(x)((ps-x)*h2*kf(x)-h*VPD*gs)^2
    px <- optimize(f2, c(minimum, ps), tol=.Machine$double.eps)$minimum
    return(px)
  }
  pxf2 <- Vectorize(pxf1)
  muf <- function(gs1){
    px1 <- pxf2(1, gs1)
    mu <- (1/2)*(-2*(-1+exp(-(-(px1/d))^c))*h3-gs1^2*LAI*(-(ca/gs1)-Km/gs1+(ca^2*gs1+gs1*Km^2+Km*Rd+ca*(2*gs1*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/(gs1*sqrt(ca^2*gs1^2+gs1^2*Km^2+(Rd-Vcmax)^2+2*ca*gs1*(gs1*Km+Rd-Vcmax)+2*gs1*(2*cp*Vcmax+Km*(Rd+Vcmax)))))+LAI*((-ca)*gs1-gs1*Km+Rd-Vcmax+sqrt(ca^2*gs1^2+gs1^2*Km^2+(Rd-Vcmax)^2+2*ca*gs1*(gs1*Km+Rd-Vcmax)+2*gs1*(2*cp*Vcmax+Km*(Rd+Vcmax))))+(2*c*gs1*h*h3*(-(px1/d))^c*VPD)/(h2*kmax*px1+c*h2*kmax*pe*(-(px1/d))^c-c*h2*kmax*px1*(-(px1/d))^c))
    return(mu)
  }
  mu <- muf(gs1)
  
  gsf <- function(w){
    w <- w
    O <- function(gs){
      px <- pxf2(w, gs)
      integrand <- function(gs)(h^2*VPD^2*((LAI*((2*ca)/gs+(2*Km)/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)^2/(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))^(3/2)+(-3*ca^2*gs-3*gs*Km^2-2*ca*(3*gs*Km+Rd-Vcmax)-2*(2*cp*Vcmax+Km*(Rd+Vcmax)))/(gs*sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax))))))/(h^2*VPD^2)-(2*(-1+c)*c*exp((-(px/d))^c)*h3*(-(px/d))^c*w^(2*b))/(h2^2*kmax^2*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(2*c^2*exp((-(px/d))^c)*h3*(-(px/d))^(2*c)*w^(2*b))/(h2^2*kmax^2*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(4*c*h3*(-(px/d))^c*w^b)/(gs*h*h2*kmax*VPD*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))+(2*(LAI*(-(ca/gs)-Km/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/(gs*sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))-(2*c*h*h3*(-(px/d))^c*VPD*w^b)/(gs*h2*kmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b)))))/(h^2*VPD^2)))/(gs*(-((2*gamma*((-1+exp(-(-(px/d))^c))*h3+mu-(1/2)*LAI*((-ca)*gs-gs*Km+Rd-Vcmax+sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax))))))/gs^2)+(2*b*(-1+c)*c*h*h3*pe*(-(px/d))^c*VPD*w^(-1+b))/(gs*h2*kmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)-(4*b*c^2*h*h3*pe*(-(px/d))^(2*c)*VPD*w^(-1+b))/(gs*h2*kmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))^2)+(2*b*c*h3*pe*(-(px/d))^c)/(exp((-(px/d))^c)*(gs^2*w*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))))+(-gamma+k/(gs*h*VPD))*(-((LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))/gs)-(2*c*h*h3*(-(px/d))^c*VPD*w^b)/(gs*h2*kmax*(px*w^b+c*(-(px/d))^c*(pe-px*w^b))))))
      res <- integrate(integrand, gs, gs1, rel.tol=.Machine$double.eps^0.25)$value-(1-w)
      return(res^2)
    }
    u <- optimize(O, c(0, ESS3(w)), tol=.Machine$double.eps)$minimum
    return(u)
  }
  
  res <- gsf(x)
  return(res)
}
gsf2 <- Vectorize(gsf1)

wleft1 <- optimize(gsf2, c(0.6, 1), tol=.Machine$double.eps)
wleft2 <- wleft1$minimum
#curve(gsf2, wleft2, 1)
wleft1
# PDF
averBf1 <- function(ca, k, MAP,
                    Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                    pe=-1.58*10^-3, b=4.38, h2=h/1000, kmax=5, c=5.71, d=10.05, h3=10,
                    gamma=1/((MAP/365/k)/1000)*nZ){
  
  mf <- function(w, gs){
    ps <- pe*w^(-b)
    kf <- function(x)kmax*exp(-(-x/d)^c)
    f1 <- function(x)-((ps-x)*h2*kf(x))
    minimum <- optimize(f1, c(-20,0), tol = .Machine$double.eps)$minimum
    f2 <- function(x)((ps-x)*h2*kf(x)-h*VPD*gs)^2
    px <- optimize(f2, c(minimum, ps), tol = .Machine$double.eps)$minimum
    
    k <- kf(px)
    PLC <- 1-k/kmax
    res <- h3*PLC
    return(res)
  }
  Af <- function(gs)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))
  Bf1 <- function(w, gs)Af(gs)-mf(w, gs)
  Bf2 <- Vectorize(Bf1)
  Bf3 <- function(w)Bf2(w, gsf2(w))
  
  Ev <- function(w){h*VPD*gsf2(w)}
  rEv <- function(w){1/Ev(w)}
  integralrEv <- Vectorize(function(w){integrate(rEv, w, 1)$value})
  fnoc <- function(w){1/Ev(w)*exp(-gamma*w-k*integralrEv(w))}
  integralfnoc <- integrate(fnoc, wleft2, 1)$value
  cpdf <- 1/(integralfnoc+1/k)
  f <- function(w){cpdf*fnoc(w)}
  Bfff <- function(w)Bf3(w)*f(w)
  averB <- integrate(Bfff, wleft2, 1)$value
  return(averB)
}

averBf2 <- Vectorize(averBf1)
averB <- averBf2(ca1, k1, MAP1)
