
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
        integrand <- function(gs){
          (gs*h*((2*Amax*gs)/(gs+Kh)^3-(2*Amax)/(gs+Kh)^2)*VPD)/((-gamma)*h*(b-gs*hk+(Amax*gs)/(gs+Kh))*VPD+(-hk-(Amax*gs)/(gs+Kh)^2+Amax/(gs+Kh))*(-k+gamma*gs*h*VPD))
          }
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
  res <- optimize(bf2, int)
  return(res)
}

ca1 <- 400
k1 <- 0.05
MAP1 <- 1825
hk1 <- 60
b1 <- bf(k=k1, MAP=MAP1, hk=hk1)
b1
