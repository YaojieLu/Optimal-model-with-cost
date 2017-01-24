
######################### No cost #########################

# Af(gs)
Af <- function(gs, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))

# dgs/dw
dgsdwf <- function(w, y, parms, mu){
  with(as.list(c(y, parms, mu)), {
    h <- l*a*LAI/nZ*p
    gamma <- 1/((MAP/365/k)/1000)*nZ

    res <- (((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))*(ca^2*gs*k*LAI+2*ca*gs*k*Km*LAI+gs*k*Km^2*LAI+ca*k*LAI*Rd+k*Km*LAI*Rd-ca*k*LAI*Vcmax+2*cp*k*LAI*Vcmax+k*Km*LAI*Vcmax-ca*k*LAI*sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))-k*Km*LAI*sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))+ca*gamma*gs*h*LAI*Rd*VPD+gamma*gs*h*Km*LAI*Rd*VPD+gamma*h*LAI*Rd^2*VPD-ca*gamma*gs*h*LAI*Vcmax*VPD+2*cp*gamma*gs*h*LAI*Vcmax*VPD+gamma*gs*h*Km*LAI*Vcmax*VPD-2*gamma*h*LAI*Rd*Vcmax*VPD+gamma*h*LAI*Vcmax^2*VPD-2*gamma*h*mu*sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))*VPD+gamma*h*LAI*Rd*sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))*VPD-gamma*h*LAI*Vcmax*sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))*VPD))/(4*gs*h*(cp+Km)*LAI*Vcmax*(Km*Rd+ca*(Rd-Vcmax)+cp*Vcmax)*VPD)
    return(list(c(res)))
  })
}

# Boundary condition at w=wL
gsBCwLf <- function(w, parms){
  with(as.list(c(w, parms)), {
    h <- l*a*LAI/nZ*p
    
    f1 <- function(gs){
      res <- (LAI*(ca*gs*(Rd-Vcmax)+gs*(2*cp*Vcmax+Km*(Rd+Vcmax))+(Rd-Vcmax)*(Rd-Vcmax+sqrt(ca^2*gs^2+gs^2*Km^2+(Rd-Vcmax)^2+2*ca*gs*(gs*Km+Rd-Vcmax)+2*gs*(2*cp*Vcmax+Km*(Rd+Vcmax))))))/(2*sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))
      return(res)
    }
    
    res <- uniroot(f1, c(0, 1), tol=.Machine$double.eps)
    return(res)
  })
}

# Boundary condition at w=1
gsBC1f <- function(parms, mu){
  with(as.list(c(parms, mu)), {
    
    w <- 1
    h <- l*a*LAI/nZ*p

    f1 <- function(gs){
      res <- (1/2)*(-2*mu+gs*LAI*(ca+Km+((-ca^2)*gs-gs*Km^2-Km*Rd-2*cp*Vcmax-Km*Vcmax+ca*(-2*gs*Km-Rd+Vcmax))/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))+LAI*((-ca)*gs-gs*Km+Rd-Vcmax+sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))))
      return(res)
    }
    
    res <- uniroot(f1, c(0, 1), tol=.Machine$double.eps)
    return(res)
  })
}

# optimize mu
muf <- function(mu){
  gs1 <- gsBC1f(parms, mu)
  out <- ode(y=c(gs=gs1$root), times=c(1, wL), func=dgsdwf, parms=parms, mu=mu)#, maxsteps=100000, method=c("lsodes")
  res <- out[2, 2]-gswL
  return(res)
}

# gs(w)
gswf <- function(w, mu){
  res <- ode(y=c(gs=gswL), times=c(wL, w), func=dgsdwf, parms=parms, mu=mu)[2, 2]
  return(res)
}
