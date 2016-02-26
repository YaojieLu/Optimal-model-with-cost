
options(digits=20)

ca1 <- 400
k1 <- 0.025
MAP1 <- 365

# w0
w0f <- function(w,
                ca=400, Vcmax=50, cp=30, Km=703, Rd=1,
                a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                pe=-1.58, b=4.38, h2=h/1000, kmax=5, c=5.71, d=10.05, h3=10){
  res <- ((-c)*h*h2*h3*(ca+Km)*kmax*LAI*pe*(ca*(Rd-Vcmax)+2*cp*Vcmax+Km*(Rd+Vcmax))*VPD*w^b*(-(pe/(w^b*d)))^c-c^2*h^2*h3^2*(ca*(Rd-Vcmax)+2*cp*Vcmax+Km*(Rd+Vcmax))*VPD^2*w^(2*b)*(-(pe/(w^b*d)))^(2*c)-sqrt(c*h*h3*(cp+Km)*Vcmax*(Km*Rd+ca*(Rd-Vcmax)+cp*Vcmax)*VPD*w^b*(-(pe/(w^b*d)))^c*(h2*(ca+Km)*kmax*LAI*pe+c*h*h3*VPD*w^b*(-(pe/(w^b*d)))^c)*(h2*(ca+Km)*kmax*LAI*pe+2*c*h*h3*VPD*w^b*(-(pe/(w^b*d)))^c)^2))/(w^b*(-(pe/(w^b*d)))^c)/(c*h*h3*(ca+Km)^2*VPD*(h2*(ca+Km)*kmax*LAI*pe+c*h*h3*VPD*w^b*(-(pe/(w^b*d)))^c))
  return(res^2)
}

w0 <- optimize(w0f, c(0.1, 0.9), tol=.Machine$double.eps)$minimum

## mu
muf <- function(ca, k, MAP){
  
  muf1 <- function(mu,
                   pe=-1.58, b=4.38, kmax=5, c=5.71, d=10.05, h3=10,
                   a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                   h2=h/1000,
                   Vcmax=50, cp=30, Km=703, Rd=1,
                   MDP=MAP/365, gamma=1/(MDP/k/1000)*nZ){
    
    pxf1 <- function(w, gs){
      ps <- pe*(w0+(1-w0)*w)^(-b)
      kf <- function(x)kmax*exp(-(-x/d)^c)
      f1 <- function(x)-((ps-x)*h2*kf(x))
      minimum <- optimize(f1, c(-20,0), tol=.Machine$double.eps)$minimum
      f2 <- function(x)((ps-x)*h2*kf(x)-h*VPD*gs)^2
      px <- optimize(f2, c(minimum, ps), tol=.Machine$double.eps)$minimum
      return(px)
    }
    pxf2 <- Vectorize(pxf1)
    
    f1 <- function(w){
      w <- w
      O <- function(gs){
        px <- pxf2(w, gs)
        integrand <- function(gs)(h^2*VPD^2*((LAI*((2*ca)/gs+(2*Km)/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)^2/(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))^(3/2)+(-3*ca^2*gs-3*gs*Km^2-2*ca*(3*gs*Km+Rd-Vcmax)-2*(2*cp*Vcmax+Km*(Rd+Vcmax)))/(gs*sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax))))))/(h^2*VPD^2)-(2*(-1+c)*c*exp((-(px/d))^c)*h3*(-(px/d))^c*(w+w0-w*w0)^(2*b))/(h2^2*kmax^2*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)+(2*c^2*exp((-(px/d))^c)*h3*(-(px/d))^(2*c)*(w+w0-w*w0)^(2*b))/(h2^2*kmax^2*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)+(4*c*h3*(-(px/d))^c*(w+w0-w*w0)^b)/(gs*h*h2*kmax*VPD*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b)))+(2*(LAI*(-(ca/gs)-Km/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/(gs*sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))-(2*c*h*h3*(-(px/d))^c*VPD*(w+w0-w*w0)^b)/(gs*h2*kmax*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b)))))/(h^2*VPD^2)))/(gs*(-((2*gamma*((-1+exp(-(-(px/d))^c))*h3+mu-(1/2)*LAI*((-ca)*gs-gs*Km+Rd-Vcmax+sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax))))))/gs^2)-(2*b*(-1+c)*c*h*h3*pe*(-(px/d))^c*VPD*(-1+w0)*(w+w0-w*w0)^(-1+b))/(gs*h2*kmax*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)+(2*b*c^2*h*h3*pe*(-(px/d))^(2*c)*VPD*(-1+w0)*(w+w0-w*w0)^(-1+b))/(gs*h2*kmax*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)-(2*b*c^2*h*h3*pe*(-(px/d))^(2*c)*VPD*(-1+w0)*(w+w0-w*w0)^b)/(gs*h2*kmax*(w*(-1+w0)-w0)*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)-(2*b*c*h3*pe*(-(px/d))^c*(-1+w0))/(exp((-(px/d))^c)*(gs^2*(w+w0-w*w0)*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))))+(-gamma+k/(gs*h*VPD))*(-((LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))/gs)-(2*c*h*h3*(-(px/d))^c*VPD*(w+w0-w*w0)^b)/(gs*h2*kmax*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))))))
        res <- integrate(integrand, 0, gs, rel.tol=.Machine$double.eps^0.75)$value-(w-w0)
        return(res^2)
      }
      u <- optimize(O, c(0,1), tol=.Machine$double.eps)$minimum
      return(u)
    }
    
    f2 <- Vectorize(f1)
    gs1 <- f2(1)
    px1 <- pxf2(1, gs1)
    res <- (1/2)*(-2*(-1+exp(-(-(px1/d))^c))*h3-gs1^2*LAI*(-(ca/gs1)-Km/gs1+(ca^2*gs1+gs1*Km^2+Km*Rd+ca*(2*gs1*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/(gs1*sqrt(ca^2*gs1^2+gs1^2*Km^2+(Rd-Vcmax)^2+2*ca*gs1*(gs1*Km+Rd-Vcmax)+2*gs1*(2*cp*Vcmax+Km*(Rd+Vcmax)))))+LAI*((-ca)*gs1-gs1*Km+Rd-Vcmax+sqrt(ca^2*gs1^2+gs1^2*Km^2+(Rd-Vcmax)^2+2*ca*gs1*(gs1*Km+Rd-Vcmax)+2*gs1*(2*cp*Vcmax+Km*(Rd+Vcmax))))+(2*c*gs1*h*h3*(-(px1/d))^c*VPD)/(h2*kmax*px1+c*h2*kmax*pe*(-(px1/d))^c-c*h2*kmax*px1*(-(px1/d))^c))
    return(res)
  }
  
  muf2 <- function(mu){(muf1(mu)-mu)^2}
  resf <- function(x)optimize(muf2, c(x, 1e9), tol=.Machine$double.eps)
  
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

mu1 <- muf(ca=ca1, k=k1, MAP=MAP1)
mu2 <- mu1$minimum

# gsw
gswf1 <- function(w,
                  mu=mu2,
                  ca=ca1, k=k1, MAP=MAP1,
                  pe=-1.58, b=4.38, kmax=5, c=5.71, d=10.05, h3=10,
                  a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                  h2=h/1000,
                  Vcmax=50, cp=30, Km=703, Rd=1,
                  MDP=MAP/365, gamma=1/(MDP/k/1000)*nZ){
  
  pxf1 <- function(w, gs){
    ps <- pe*(w0+(1-w0)*w)^(-b)
    kf <- function(x)kmax*exp(-(-x/d)^c)
    f1 <- function(x)-((ps-x)*h2*kf(x))
    minimum <- optimize(f1, c(-20,0), tol=.Machine$double.eps)$minimum
    f2 <- function(x)((ps-x)*h2*kf(x)-h*VPD*gs)^2
    px <- optimize(f2, c(minimum, ps), tol=.Machine$double.eps)$minimum
    return(px)
  }
  pxf2 <- Vectorize(pxf1)
  
  f1 <- function(w){
    w <- w
    O <- function(gs){
      px <- pxf2(w, gs)
      integrand <- function(gs)(h^2*VPD^2*((LAI*((2*ca)/gs+(2*Km)/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)^2/(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))^(3/2)+(-3*ca^2*gs-3*gs*Km^2-2*ca*(3*gs*Km+Rd-Vcmax)-2*(2*cp*Vcmax+Km*(Rd+Vcmax)))/(gs*sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax))))))/(h^2*VPD^2)-(2*(-1+c)*c*exp((-(px/d))^c)*h3*(-(px/d))^c*(w+w0-w*w0)^(2*b))/(h2^2*kmax^2*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)+(2*c^2*exp((-(px/d))^c)*h3*(-(px/d))^(2*c)*(w+w0-w*w0)^(2*b))/(h2^2*kmax^2*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)+(4*c*h3*(-(px/d))^c*(w+w0-w*w0)^b)/(gs*h*h2*kmax*VPD*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b)))+(2*(LAI*(-(ca/gs)-Km/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/(gs*sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))-(2*c*h*h3*(-(px/d))^c*VPD*(w+w0-w*w0)^b)/(gs*h2*kmax*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b)))))/(h^2*VPD^2)))/(gs*(-((2*gamma*((-1+exp(-(-(px/d))^c))*h3+mu-(1/2)*LAI*((-ca)*gs-gs*Km+Rd-Vcmax+sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax))))))/gs^2)-(2*b*(-1+c)*c*h*h3*pe*(-(px/d))^c*VPD*(-1+w0)*(w+w0-w*w0)^(-1+b))/(gs*h2*kmax*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)+(2*b*c^2*h*h3*pe*(-(px/d))^(2*c)*VPD*(-1+w0)*(w+w0-w*w0)^(-1+b))/(gs*h2*kmax*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)-(2*b*c^2*h*h3*pe*(-(px/d))^(2*c)*VPD*(-1+w0)*(w+w0-w*w0)^b)/(gs*h2*kmax*(w*(-1+w0)-w0)*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)-(2*b*c*h3*pe*(-(px/d))^c*(-1+w0))/(exp((-(px/d))^c)*(gs^2*(w+w0-w*w0)*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))))+(-gamma+k/(gs*h*VPD))*(-((LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))/gs)-(2*c*h*h3*(-(px/d))^c*VPD*(w+w0-w*w0)^b)/(gs*h2*kmax*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))))))
      res <- integrate(integrand, 0, gs)$value-(w-w0)
      return(res^2)
    }
    u <- optimize(O, c(0,1), tol=.Machine$double.eps)$minimum
    return(u)
  }
  
  f2 <- Vectorize(f1)
  res <- f2(w)
  return(res)
}
gswf2 <- Vectorize(gswf1)

# Figure
windows(8, 6)
par(mgp = c(2.2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(4, 4, 2.5, 4), mfrow=c(1,1))
curve(gswf2, w0+0.032, 1, type="l",
     xlab = expression(italic(w)~"(%)"),
     ylab = expression(italic(g[s])~(mol~m^-2~s^-1)),
     xlim = c(0, 1), ylim = c(0, 1),
     cex.lab = 1.3)
abline(v=w0)