
slope <- 65

bf1 <- function(b,
                k=0.05, MDP=5,
                Amax=15.85545680, Kh=0.03530976, hk=slope, VPD=0.02, p=43200, LAI=1, nZ=0.5,
                a=1.6, l=1.8e-5, gamma=1/((MDP/k)/1000)*nZ, h=l*a*LAI/nZ*p){
  
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
  #res <- gs1*hk-(Amax*gs1)/(gs1+Kh)+gs1^2*h*(-hk-(Amax*gs1)/(gs1+Kh)^2+Amax/(gs1+Kh))*VPD
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
      message(a+b*i/(10^c))
    }
  }
  return(a+b*i/(10^c))
}

intlower1 <- estint((slope-65)*1.5/5-1, 1, 3)#!!! sensitive to environmental conditions
intupper1 <- estint((slope-65)*0.5/5, -1, 3)#!!! sensitive to environmental conditions
int <- c(intlower1, intupper1)

bf2 <- function(b){(bf1(b)-b)^2}
optimize(bf2, int)$minimum