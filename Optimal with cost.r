
# Functions
# With cost
bf <- function(k, MAP,
               hk){
  
  bf1 <- function(b,
                  Amax=15.85545680, Kh=0.03530976, VPD=0.02, p=43200, LAI=1, nZ=0.5,
                  a=1.6, l=1.8e-5, MDP=MAP/365, gamma=1/((MDP/k)/1000)*nZ, h=l*a*LAI/nZ*p){
    
    f1 <- function(x){
      w <- x
      O <- function(gs){
        integrand <- function(gs){-((2*Amax*gs*h*Kh*VPD)/((gs+Kh)*(2*gs*Kh*(hk*k-b*gamma*h*VPD)+gs^2*(hk*k-(b+Amax)*gamma*h*VPD)+Kh*((-Amax)*k+hk*k*Kh-b*gamma*h*Kh*VPD))))}
        integrate(integrand, 0, gs)$value-w
      }
      u <- uniroot(O, c(0,1), tol=1e-10, extendInt = "yes")
      return(u$root)
    }
    
    f2 <- Vectorize(f1)
    gs1 <- f2(1)
    res <- -((Amax*gs1^2)/(gs1+Kh)^2)
    return(res)
  }
  
  estint <- function(a, b, c){
    estb <- try(bf1(a), silent=TRUE)
    i <- 0
    if(inherits(estb, "try-error")){
      while(inherits(estb, "try-error")){
        i <- i + 1
        estb <- try(bf1(a+b*i/(10^c)), silent=TRUE)
      }
    }
    return(a+b*i/(10^c))
  }
  
  intlower1 <- estint(-20, 1, 3)#!!! sensitive to environmental conditions
  intupper1 <- estint(0, -1, 3)#!!! sensitive to environmental conditions
  int <- c(intlower1, intupper1)
  
  bf2 <- function(b){(bf1(b)-b)^2}
  res <- optimize(bf2, int)$minimum
  return(res)
}

ca1 <- 400
k1 <- 0.1
MAP1 <- 365
hk1 <- 0
b1 <- bf(k=k1, MAP=MAP1, hk=hk1)

gsw1 <- function(w,
                 k, MAP,
                 hk, b,
                 Amax=15.85545680, Kh=0.03530976, VPD=0.02, p=43200, LAI=1, nZ=0.5,
                 a=1.6, l=1.8e-5, MDP=MAP/365, gamma=1/((MDP/k)/1000)*nZ, h=l*a*LAI/nZ*p){
  w <- w
  O <- function(gs){
    integrand <- function(gs){-((2*Amax*gs*h*Kh*VPD)/((gs+Kh)*(2*gs*Kh*(hk*k-b*gamma*h*VPD)+gs^2*(hk*k-(b+Amax)*gamma*h*VPD)+Kh*((-Amax)*k+hk*k*Kh-b*gamma*h*Kh*VPD))))}
    integrate(integrand, 0, gs)$value-w
  }
  u <- uniroot(O, c(0,1), tol=1e-10, extendInt = "yes")
  return(u$root)
}

gsw2 <- function(w)gsw1(w, k=k1, MAP=MAP1, hk=hk1, b=b1)
gsw3 <- Vectorize(gsw2)

# Without cost
gwopt <- function(averA,
                  ca, k, MAP,
                  Vcmax=50, cp=30, Km=703, Rd=1, VPD=0.02, p=43200, LAI=1, nZ=0.5,
                  a=1.6, l=1.8e-5, MDP=MAP/365, gamma=1/((MDP/k)/1000)*nZ, h=l*a*LAI/nZ*p){
  
  f1 <- function(w){
    w <- w
    O <- function(g){
      integrand <- function(g){(4*h*(cp+Km)*LAI*Vcmax*(Km*Rd+ca*(Rd-Vcmax)+cp*Vcmax)*VPD*g)/(((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))*(ca*k*LAI*Rd+k*Km*LAI*Rd-ca*k*LAI*Vcmax+2*cp*k*LAI*Vcmax+k*Km*LAI*Vcmax+gamma*h*LAI*Rd^2*VPD-2*gamma*h*LAI*Rd*Vcmax*VPD+gamma*h*LAI*Vcmax^2*VPD+ca^2*k*LAI*g+2*ca*k*Km*LAI*g+k*Km^2*LAI*g+ca*gamma*h*LAI*Rd*VPD*g+gamma*h*Km*LAI*Rd*VPD*g-ca*gamma*h*LAI*Vcmax*VPD*g+2*cp*gamma*h*LAI*Vcmax*VPD*g+gamma*h*Km*LAI*Vcmax*VPD*g-ca*k*LAI*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))-k*Km*LAI*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))+2*gamma*h*averA*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))+gamma*h*LAI*Rd*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))-gamma*h*LAI*Vcmax*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))))}
      integrate(integrand, 0, g)$value-w
    }
    u <- uniroot(O, c(0,1), tol=1e-10, extendInt = "yes")
    return(u$root)
  }
  
  f2 <- Vectorize(f1)
  
  g <- vector("numeric", length=length(w))
  g <- f2(w=w)
  
  return(g)
}

dvs <- read.csv("Derived variables.csv")
averAold <- dvs[dvs$ca==ca1 & dvs$k==k1 & dvs$MAP==MAP1, "A"]
w <- seq(0, 1, length=1e3)
gopt <- gwopt(averAold, ca1, k1, MAP1)

# ESS
gX <- (sqrt(15.85545680*0.03530976)-sqrt(hk1)*0.03530976)/sqrt(hk1)

# Figure
windows(8, 6)
par(mgp = c(2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(3.5, 4, 1, 2))
# With cost
curve(gsw3, 0, 1,
      xlab = expression(italic(w)~"(%)"), 
      ylab = expression(italic(g[s])~(mol~m^-2~s^-1)),
      xlim = c(0, 1), ylim = c(0, 0.1),
      cex.lab = 1.3)
# Without cost
points(w, gopt, type = "l", col = "red")
# ESS
abline(h=gX, col="blue")
legend("topleft", legend=c(as.expression(bquote(cost==.(hk1)%*%italic(g[s]))), "cost = 0", "ESS"), col = c("black", "red", "blue"), lty=c(1, 1, 1), lwd = c(2, 2, 2))

# averA
averAf1 <- function(k, MAP, 
                   b, hk,
                   Amax=15.85545680, Kh=0.03530976, VPD=0.02, p=43200, LAI=1, nZ=0.5,
                   a=1.6, l=1.8e-5, MDP=MAP/365, gamma=1/((MDP/k)/1000)*nZ, h=l*a*LAI/nZ*p){
  
  gsw <- function(w)gsw3(w)
  A <- function(w){Amax*gsw(w)/(gsw(w)+Kh)-hk*gsw(w)}
  Ev <- function(w){h*VPD*gsw(w)}
  rEv <- function(w){1/Ev(w)}
  integralEv <- Vectorize(function(w){integrate(rEv, 0, w)$value})
  fnoc <- function(w){1/Ev(w)*exp(-gamma*w+k*integralEv(w))}
  integralfnoc <- integrate(fnoc, 0, 1)$value
  c <- 1/(integralfnoc+1/k)
  f <- function(w){c*fnoc(w)}
  fA <- function(w){f(w)*A(w)}
  averA <- integrate(fA, 0, 1)$value
  return(averA)
}

averAf2 <- function(k, MAP, 
                    b, hk,
                    Amax=15.85545680, Kh=0.03530976, VPD=0.02, p=43200, LAI=1, nZ=0.5,
                    a=1.6, l=1.8e-5, MDP=MAP/365, gamma=1/((MDP/k)/1000)*nZ, h=l*a*LAI/nZ*p){
  
  gsw <- function(w)ifelse(gsw3(w)<gX, gsw3(w), gX)
  A <- function(w){Amax*gsw(w)/(gsw(w)+Kh)-hk*gsw(w)}
  Ev <- function(w){h*VPD*gsw(w)}
  rEv <- function(w){1/Ev(w)}
  integralEv <- Vectorize(function(w){integrate(rEv, 0, w)$value})
  fnoc <- function(w){1/Ev(w)*exp(-gamma*w+k*integralEv(w))}
  integralfnoc <- integrate(fnoc, 0, 1)$value
  c <- 1/(integralfnoc+1/k)
  f <- function(w){c*fnoc(w)}
  fA <- function(w){f(w)*A(w)}
  averA <- integrate(fA, 0, 1)$value
  return(averA)
}

averA21 <- averAf1(k1, MAP1, b1, hk1)
averA22 <- averAf2(k1, MAP1, b1, hk1)
averA21
averA22

# A, h, B
#A <- function(gs,Amax=15.85545680, Kh=0.03530976)Amax*gs/(Kh+gs)
#h <- function(gs, s=hk)s*gs
#B <- function(gs)A(gs) - h(gs)

#windows(8, 6)
#par(mgp = c(2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(3.5, 4, 1, 2))
#curve(B, 0, gsw3(1), ylim=c(0, 10))
#abline(v=gX)