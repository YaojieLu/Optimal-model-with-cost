
######################### No cost #########################

# optimize mu
muf <- function(mu, parms){
  with(as.list(c(parms, mu)), {
    h <- l*a*LAI/nZ*p
    gamma <- 1/((MAP/365/k)/1000)*nZ
    
    # Boundary condition at w=1
    gs1f <- function(mu){
      f1 <- function(gs)ca*gs*LAI+gs*Km*LAI+2*mu-LAI*Rd+LAI*Vcmax-LAI*sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))-gs*LAI*(ca+Km-(ca^2*gs+gs*Km^2+Km*Rd+ca*(2*gs*Km+Rd-Vcmax)+2*cp*Vcmax+Km*Vcmax)/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))
      res <- uniroot(f1, c(0, 1), tol=.Machine$double.eps)
      return(res)
    }
    
    # integrand
    f1 <- function(gs)1/((((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))*(ca^2*gs*k*LAI+2*ca*gs*k*Km*LAI+gs*k*Km^2*LAI+ca*k*LAI*Rd+k*Km*LAI*Rd-ca*k*LAI*Vcmax+2*cp*k*LAI*Vcmax+k*Km*LAI*Vcmax-ca*k*LAI*sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))-k*Km*LAI*sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))+ca*gamma*gs*h*LAI*Rd*VPD+gamma*gs*h*Km*LAI*Rd*VPD+gamma*h*LAI*Rd^2*VPD-ca*gamma*gs*h*LAI*Vcmax*VPD+2*cp*gamma*gs*h*LAI*Vcmax*VPD+gamma*gs*h*Km*LAI*Vcmax*VPD-2*gamma*h*LAI*Rd*Vcmax*VPD+gamma*h*LAI*Vcmax^2*VPD-2*gamma*h*mu*sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))*VPD+gamma*h*LAI*Rd*sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))*VPD-gamma*h*LAI*Vcmax*sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax))*VPD))/(4*gs*h*(cp+Km)*LAI*Vcmax*(Km*Rd+ca*(Rd-Vcmax)+cp*Vcmax)*VPD))
    
    gs1 <- gs1f(mu)$root
    #browser()
    res <- integrate(f1, gsw0, gs1, rel.tol=.Machine$double.eps^0.25)$value-(1-w0)
    return(res^2)
  })
}
