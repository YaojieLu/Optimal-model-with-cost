
gw <- function(q){
  f1 <- function(w,
                 ca=400, lambda=0.05, averR=50,
                 Vcmax=50,cp=30, Km=703, Rd=1, VPD=0.02, p=43200,LAI=1, nZ=0.5,
                 a=1.6, l=1.8e-5, gamma=1/(averR/1000)*nZ,h=l*a*LAI/nZ*p){
    w <- w
    O <- function(g){
      integrand <- function(g){
        1/(((-2*Rd*Vcmax+(Vcmax)^2+2*(-ca+2*cp+Km)*Vcmax*g+(Rd+(ca+Km)*g)^2)*(ca*lambda*Rd+Km*lambda*Rd-ca*lambda*Vcmax+2*cp*lambda*Vcmax+Km*lambda*Vcmax+gamma*h*Rd^2*VPD-2*gamma*h*Rd*Vcmax*VPD+gamma*h*(Vcmax)^2*VPD+ca^2*lambda*g+2*ca*Km*lambda*g+Km^2*lambda*g+ca*gamma*h*Rd*VPD*g+gamma*h*Km*Rd*VPD*g-ca*gamma*h*Vcmax*VPD*g+2*cp*gamma*h*Vcmax*VPD*g+gamma*h*Km*Vcmax*VPD*g-ca*lambda*sqrt(-2*Rd*Vcmax+(Vcmax)^2+2*(-ca+2*cp+Km)*Vcmax*g+(Rd+(ca+Km)*g)^2)-Km*lambda*sqrt(-2*Rd*Vcmax+(Vcmax)^2+2*(-ca+2*cp+Km)*Vcmax*g+(Rd+(ca+Km)*g)^2)+2*gamma*h*q*VPD*sqrt(-2*Rd*Vcmax+(Vcmax)^2+2*(-ca+2*cp+Km)*Vcmax*g+(Rd+(ca+Km)*g)^2)+gamma*h*Rd*VPD*sqrt(-2*Rd*Vcmax+(Vcmax)^2+2*(-ca+2*cp+Km)*Vcmax*g+(Rd+(ca+Km)*g)^2)-gamma*h*Vcmax*VPD*sqrt(-2*Rd*Vcmax+(Vcmax)^2+2*(-ca+2*cp+Km)*Vcmax*g+(Rd+(ca+Km)*g)^2)))/(4*h*(cp+Km)*Vcmax*(Km*Rd+ca*(Rd-Vcmax)+cp*Vcmax)*VPD*g))
      }
      res <- integrate(integrand, 0, g)$value-w
    }
    u <- uniroot(O, c(0,1), tol=1e-10, extendInt = "yes")
    return(u$root)
  }
  
  # vectorized function for curve()
  f2 <- Vectorize(f1)
  
  # soil water content
  w <- seq(0,1,length=1e4)
  
  # solution
  g <- vector("numeric", length=length(w))
  g <- sapply(w, f2)
  
  res <- data.frame(w=w, g=g)
  
  return(res)
}

t0 <- proc.time()[3]
data <- gw(11.5875204255557)
t1 <- proc.time()[3]
message("time=", round((t1-t0)/60 ,2), " min")

with(data, points(w, g, type = "l", col = "red"))

legend("bottomright", c("Optimal without cost", "Optimal with cost", "ESS"), col=c("red", "black", "blue"), lty=c(1, 1, 1))
