
ca1 <- 400
k1 <- 0.1
MAP1 <- 3650

# Functions
## mu
muf <- function(ca, k, MAP){
  
  muf1 <- function(mu,
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
    
    f1 <- function(x){
      w <- x
      O <- function(gs){
        integrand <- function(gs){
          px <- pxf2(pxw=w, pxgs=gs)
          res <- (gs^3*h*LAI*((2*ca)/gs+(2*Km)/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)^2/(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))^(3/2)+(-3*ca^2*gs-3*gs*Km^2-2*ca*(3*gs*Km+Rd-Vcmax)-2*(2*cp*Vcmax+Km*(Rd+Vcmax)))/(gs*sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))*VPD-(2*(-1+c)*c*exp((-(px/d))^c)*gs^3*h^3*h3*(-(px/d))^c*VPD^3*(w+w0-w*w0)^(2*b))/(h2^2*kmax^2*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)+(2*c^2*exp((-(px/d))^c)*gs^3*h^3*h3*(-(px/d))^(2*c)*VPD^3*(w+w0-w*w0)^(2*b))/(h2^2*kmax^2*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)+(4*c*gs^2*h^2*h3*(-(px/d))^c*VPD^2*(w+w0-w*w0)^b)/(h2*kmax*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b)))+2*gs^2*h*VPD*(gs*LAI*(-(ca/gs)-Km/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/(gs*sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))-(2*c*h*h3*(-(px/d))^c*VPD*(w+w0-w*w0)^b)/(h2*kmax*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b)))))/(gs^2*h*VPD*(-2*gamma*((-1+exp(-(-(px/d))^c))*h3+mu-(1/2)*LAI*((-ca)*gs-gs*Km+Rd-Vcmax+sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))-(2*b*(-1+c)*c*gs*h*h3*pe*(-(px/d))^c*VPD*(-1+w0)*(w+w0-w*w0)^(-1+b))/(h2*kmax*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)+(2*b*c^2*gs*h*h3*pe*(-(px/d))^(2*c)*VPD*(-1+w0)*(w+w0-w*w0)^(-1+b))/(h2*kmax*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)-(2*b*c^2*gs*h*h3*pe*(-(px/d))^(2*c)*VPD*(-1+w0)*(w+w0-w*w0)^b)/(h2*kmax*(w*(-1+w0)-w0)*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)-(2*b*c*h3*pe*(-(px/d))^c*(-1+w0))/(exp((-(px/d))^c)*((w+w0-w*w0)*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))))+gs*(-gamma+k/(gs*h*VPD))*((-LAI)*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax))))-(2*c*h*h3*(-(px/d))^c*VPD*(w+w0-w*w0)^b)/(h2*kmax*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))))))
          return(res)
        }
        res <- (integrate(integrand, 0, gs)$value-w)^2
        return(res)
      }
      u <- optimize(O, c(0,1))$minimum
      return(u)
    }
    
    f2 <- Vectorize(f1)
    gs1 <- f2(1)
    px1 <- pxf2(pxw=1, pxgs=gs1)
    z1 <- 1/(h*VPD*gs1)
    resf <- function(w)-((1/(2*z1^2))*(-((2*c*h3*(-(px1/d))^c*(w+w0-w*w0)^b*z1)/(h2*kmax*(px1*(w+w0-w*w0)^b+c*(-(px1/d))^c*(pe-px1*(w+w0-w*w0)^b))))-2*(1-exp(-(-(px1/d))^c))*h3*z1^2+(LAI*((-ca)*h*VPD*z1-h*Km*VPD*z1+(ca^2+Km^2+2*cp*h*Vcmax*VPD*z1+h*Km*(Rd+Vcmax)*VPD*z1+ca*(2*Km+h*(Rd-Vcmax)*VPD*z1))/sqrt((ca^2+Km^2+2*h*Km*(Rd+Vcmax)*VPD*z1+2*ca*(Km+h*(Rd-Vcmax)*VPD*z1)+h*VPD*z1*(4*cp*Vcmax+h*(Rd-Vcmax)^2*VPD*z1))/(h^2*VPD^2*z1^2))))/(h^2*VPD^2)+(LAI*z1*(ca+Km-h*Rd*VPD*z1+h*Vcmax*VPD*z1-h*VPD*z1*sqrt((ca^2+Km^2+2*h*Km*(Rd+Vcmax)*VPD*z1+2*ca*(Km+h*(Rd-Vcmax)*VPD*z1)+h*VPD*z1*(4*cp*Vcmax+h*(Rd-Vcmax)^2*VPD*z1))/(h^2*VPD^2*z1^2))))/(h*VPD)))
    res <- resf(1)
    return(res)
  }
  
  muf2 <- function(mu){(muf1(mu)-mu)^2}
  resf <- function(x)optimize(muf2, c(x, 1e9))
  
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

## gs(w)
gsw1 <- function(w,
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
  
  gsmax1 <- function(w){
    
    ps <- pe*(w0+(1-w0)*w)^(-b)
    kf <- function(x)kmax*exp(-(-x/d)^c)
    f1 <- function(x)-((ps-x)*h2*kf(x))
    pxmin <- optimize(f1, c(-20,0))$minimum
    f2 <- function(gs)((ps-pxmin)*h2*kf(pxmin)-h*VPD*gs)^2
    gsmax <- optimize(f2, c(0, 2))$minimum
    return(gsmax)
  }
  gsmax2 <- Vectorize(gsmax1)
  
  O <- function(gs){
    integrand <- function(gs){
      px <- pxf2(pxw=w, pxgs=gs)
      res <- (gs^3*h*LAI*((2*ca)/gs+(2*Km)/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)^2/(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))^(3/2)+(-3*ca^2*gs-3*gs*Km^2-2*ca*(3*gs*Km+Rd-Vcmax)-2*(2*cp*Vcmax+Km*(Rd+Vcmax)))/(gs*sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))*VPD-(2*(-1+c)*c*exp((-(px/d))^c)*gs^3*h^3*h3*(-(px/d))^c*VPD^3*(w+w0-w*w0)^(2*b))/(h2^2*kmax^2*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)+(2*c^2*exp((-(px/d))^c)*gs^3*h^3*h3*(-(px/d))^(2*c)*VPD^3*(w+w0-w*w0)^(2*b))/(h2^2*kmax^2*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)+(4*c*gs^2*h^2*h3*(-(px/d))^c*VPD^2*(w+w0-w*w0)^b)/(h2*kmax*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b)))+2*gs^2*h*VPD*(gs*LAI*(-(ca/gs)-Km/gs+(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/(gs*sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))-(2*c*h*h3*(-(px/d))^c*VPD*(w+w0-w*w0)^b)/(h2*kmax*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b)))))/(gs^2*h*VPD*(-2*gamma*((-1+exp(-(-(px/d))^c))*h3+mu-(1/2)*LAI*((-ca)*gs-gs*Km+Rd-Vcmax+sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax)))))-(2*b*(-1+c)*c*gs*h*h3*pe*(-(px/d))^c*VPD*(-1+w0)*(w+w0-w*w0)^(-1+b))/(h2*kmax*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)+(2*b*c^2*gs*h*h3*pe*(-(px/d))^(2*c)*VPD*(-1+w0)*(w+w0-w*w0)^(-1+b))/(h2*kmax*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)-(2*b*c^2*gs*h*h3*pe*(-(px/d))^(2*c)*VPD*(-1+w0)*(w+w0-w*w0)^b)/(h2*kmax*(w*(-1+w0)-w0)*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))^2)-(2*b*c*h3*pe*(-(px/d))^c*(-1+w0))/(exp((-(px/d))^c)*((w+w0-w*w0)*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))))+gs*(-gamma+k/(gs*h*VPD))*((-LAI)*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax))))-(2*c*h*h3*(-(px/d))^c*VPD*(w+w0-w*w0)^b)/(h2*kmax*(px*(w+w0-w*w0)^b+c*(-(px/d))^c*(pe-px*(w+w0-w*w0)^b))))))
      return(res)
    }
    res <- (integrate(integrand, 0, gs)$value-w)^2
    return(res)
  }
  
  u <- function(x)optimize(O, c(1e-10,x))$minimum
  estint <- function(p1, p2, p3){
    est <- try(u(p1), silent=TRUE)
    i <- 0
    if(inherits(est, "try-error")){
      while(inherits(est, "try-error")){
        i <- i + 1
        est <- try(u(p1+p2*i/(10^p3)), silent=TRUE)
      }
    }
    return(p1+p2*i/(10^p3))
  }
  
  intupper0 <- estint(gsmax2(w), -1, 1)
  intupper1 <- estint(max(0,intupper0)+1e-1, -1, 2)
  intupper2 <- estint(max(0,intupper1)+1e-2, -1, 3)
  intupper3 <- estint(max(0,intupper2)+1e-3, -1, 4)
  intupper4 <- estint(max(0,intupper3)+1e-4, -1, 5)
  res <- u(intupper4)
  message(w)
  return(res)
}

gsw2 <- Vectorize(gsw1)

w0 <- optimize(gsw2, c(0, 1))

# Results
wdata <- c(w0$minimum, seq(ceiling(w0$minimum*500)/500, 1, by=1/500))
gsdata <- sapply(wdata, gsw2)

res <- cbind(wdata, gsdata)
colnames(res) <- c("w", "gs") 

namef1 <- function(ca, k, MAP){ 
  filename <- paste("gs(w) for ca=",bquote(.(ca1)),"k=",bquote(.(k1)),"MAP=",bquote(.(MAP1)), ".csv", sep="") 
  write.csv(res, file = filename) 
}
namef1(ca1,k1,MAP1)

# Plots
gsw <- read.csv("gs(w).csv")
windows(8, 6)
par(mgp = c(2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(3.5, 4, 2, 2))
plot(gsw$w, gsw$gs, type="l",
     main = c(as.expression(bquote(italic(c[a])==.(ca1)~~italic(k)==.(k1)~~MAP==.(MAP1)))),
     xlab = expression(italic(w)~"(%)"),
     ylab = expression(italic(g[s])~(mol~m^-2~s^-1)),
     xlim = c(0, 1), ylim = c(0, 0.6),
     cex.lab = 1.3)

## Optimal
gw1 <- function(w,
                averA,
                ca=ca1, k=k1, MAP=MAP1,
                Vcmax=50, cp=30, Km=703, Rd=1, VPD=0.02, p=43200, LAI=1, nZ=0.5,
                a=1.6, l=1.8e-5, MDP=MAP/365, gamma=1/((MDP/k)/1000)*nZ, h=l*a*LAI/nZ*p){
  w <- w
  O <- function(g){
    integrand <- function(g){(4*h*(cp+Km)*LAI*Vcmax*(Km*Rd+ca*(Rd-Vcmax)+cp*Vcmax)*VPD*g)/(((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))*(ca*k*LAI*Rd+k*Km*LAI*Rd-ca*k*LAI*Vcmax+2*cp*k*LAI*Vcmax+k*Km*LAI*Vcmax+gamma*h*LAI*Rd^2*VPD-2*gamma*h*LAI*Rd*Vcmax*VPD+gamma*h*LAI*Vcmax^2*VPD+ca^2*k*LAI*g+2*ca*k*Km*LAI*g+k*Km^2*LAI*g+ca*gamma*h*LAI*Rd*VPD*g+gamma*h*Km*LAI*Rd*VPD*g-ca*gamma*h*LAI*Vcmax*VPD*g+2*cp*gamma*h*LAI*Vcmax*VPD*g+gamma*h*Km*LAI*Vcmax*VPD*g-ca*k*LAI*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))-k*Km*LAI*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))+2*gamma*h*averA*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))+gamma*h*LAI*Rd*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))-gamma*h*LAI*Vcmax*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))))}
    integrate(integrand, 0, g)$value-w
  }
  u <- uniroot(O, c(0,1), tol=1e-10, extendInt = "yes")
  return(u$root)
}

gw2 <- Vectorize(gw1)
optdata <- read.csv("Derived variables.csv")
gw3 <- function(w)gw2(w, averA=optdata[optdata$ca==ca1 & optdata$k==k1 & optdata$MAP==MAP1, "A"])

curve(gw3, 0, 1, add=T, col="red")
Legend = c("With PLC cost", "Without PLC cost")
legend("topleft", Legend, lty = c(1, 1), col = c("black", "red"))

namef2 <- function(ca, k, MAP){ 
  filename <- paste("ca=",bquote(.(ca1)),"k=",bquote(.(k1)),"MAP=",bquote(.(MAP1)), ".pdf", sep="") 
  dev.copy2pdf(file = filename) 
}
namef2(ca1,k1,MAP1)
