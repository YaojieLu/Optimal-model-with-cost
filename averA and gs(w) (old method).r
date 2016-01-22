
hk1 <- 0

qopt <- function(ca, k, MAP, hk=hk1){
  
  # qopt1: calculate average photosynthesis rate (A) for given q (the assumed Abar)
  qopt1 <- function(q,
                    Vcmax=50, cp=30, Km=703, Rd=1, VPD=0.02, p=43200, LAI=1, nZ=0.5,
                    a=1.6, l=1.8e-5, MDP=MAP/365, gamma=1/((MDP/k)/1000)*nZ, h=l*a*LAI/nZ*p){
    # f1(w): the solution to Eq. 14
    f1 <- function(w){
      w <- w
      O <- function(g){
        # Following Eq. 14: integrand(g)dg = dw
        integrand <- function(g){
          (4*h*(cp+Km)*LAI*Vcmax*(Km*Rd+ca*(Rd-Vcmax)+cp*Vcmax)*VPD*g)/(((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))*(ca*k*LAI*Rd+k*Km*LAI*Rd-ca*k*LAI*Vcmax+2*cp*k*LAI*Vcmax+k*Km*LAI*Vcmax+gamma*h*LAI*Rd^2*VPD-2*gamma*h*LAI*Rd*Vcmax*VPD+gamma*h*LAI*Vcmax^2*VPD+ca^2*k*LAI*g+2*ca*k*Km*LAI*g+k*Km^2*LAI*g+ca*gamma*h*LAI*Rd*VPD*g+gamma*h*Km*LAI*Rd*VPD*g-ca*gamma*h*LAI*Vcmax*VPD*g+2*cp*gamma*h*LAI*Vcmax*VPD*g+gamma*h*Km*LAI*Vcmax*VPD*g-ca*k*LAI*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))+2*hk*k*LAI*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))-k*Km*LAI*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))+2*gamma*h*q*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))+gamma*h*LAI*Rd*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))-gamma*h*LAI*Vcmax*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))))
        }
        integrate(integrand, 0, g)$value-w
      }
      u <- uniroot(O, c(0,1), tol=1e-10, extendInt = "yes")
      return(u$root)
    }
    
    f2 <- Vectorize(f1)
    
    # Eq. 7
    A <- function(w){LAI*(1/2*(Vcmax+(Km+ca)*f2(w)-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*f2(w)+((ca+Km)*f2(w)+Rd)^2-2*Rd*Vcmax)^(1/2))+hk*f2(w))}
    # Eq. 1
    Ev <- function(w){h*VPD*f2(w)}
    rEv <- function(w){1/Ev(w)}
    integralEv <- Vectorize(function(w){integrate(rEv, 0, w)$value})
    fnoc <- function(w){1/Ev(w)*exp(-gamma*w+k*integralEv(w))}
    integralfnoc <- integrate(fnoc, 0, 1)$value
    # c: the integral constant in Eqs. 3 & 4
    c <- 1/(integralfnoc+1/k)
    # Eq. 4
    f <- function(w){c*fnoc(w)}
    
    fA <- function(w){f(w)*A(w)}
    # averA: average A
    averA <- integrate(fA, 0, 1)$value
    return(averA)
  }
  
  #estint: identify the interval of q, within which "optim" gives meaningful result
  estint <- function(a, b, c){
    estA <- try(qopt1(a), silent=TRUE)
    if(inherits(estA, "try-error")){
      i <- 0
      while(inherits(estA, "try-error")){
        i <- i + 1
        estA <- try(qopt1(a+b*i/(10^c)), silent=TRUE)
      }
    }
    return(a+b*i/(10^c))
  }
  
  int <- c()
  intlower1 <- estint(0, 1, 3) # the lower bound
  # We guess that the real Abar should be smaller than 30 under
  # all the combinations of environmental conditions considered here.
  # So if you change the environmental conditions into those more
  # favourable to photosynthesis, you may need to increase this limit
  intupper1 <- estint(15, -1, 3) # the upper bound
  
  # Estimating interval with higher accuracy if necessary
  ifelse(
    round(intlower1, 3)==round(intupper1, 3),
    {
      intlower2 <- estint(intlower1-10^-3, 1, 6)
      intupper2 <- estint(intupper1+10^-3, -1, 6)
      int[1] <- intlower2
      int[2] <- intupper2
    },
    {
      int[1] <- intlower1
      int[2] <- intupper1
    }
  )
  
  #qopt2: square of the difference between q and the resulting averA
  qopt2 <- function(q){(qopt1(q)-q)^2}
  estA <- optimize(qopt2, int)$minimum
  
  return(estA)
}

#environmental conditions
ca <- c(400)  # Atmospheric CO2 concentration (ppm)
k <- c(0.1) # Rainfall frequency (per day)
MAP <- seq(1, 10, by=1)*365 # MAP=MDP*365; MAP: mean annual precipitation; MDP: mean daily precipitation
env <- as.vector(expand.grid(ca, k, MAP))

# Initialize
dvs <- matrix(nrow=nrow(env), ncol=1)

# Run every parameter combination
for(i in 1:nrow(env)){
  # Identify the correct value for Abar under given environmental conditions
  dvs[i, 1] <- qopt(env[i, 1], env[i, 2], env[i, 3])
  message(env[i, 1], " ", env[i, 2], " ", env[i, 3], " ", dvs[i, 1])
}

# Collect results
res <- cbind(env, dvs)
colnames(res) <- c("ca", "k", "MAP", "A") 

write.csv(res, "Derived variables with cost.csv", row.names = FALSE)

# gw: the optimal gs(w), following Eq. 14
gw1 <- function(w,
               ca, k, MAP,
               averA,
               hk=hk1,
               Vcmax=50, cp=30, Km=703, Rd=1, VPD=0.02, p=43200, LAI=1, nZ=0.5,
               a=1.6, l=1.8e-5, MDP=MAP/365, gamma=1/((MDP/k)/1000)*nZ, h=l*a*LAI/nZ*p){
  w <- w
  O <- function(g){
    integrand <- function(g){(4*h*(cp+Km)*LAI*Vcmax*(Km*Rd+ca*(Rd-Vcmax)+cp*Vcmax)*VPD*g)/(((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))*(ca*k*LAI*Rd+k*Km*LAI*Rd-ca*k*LAI*Vcmax+2*cp*k*LAI*Vcmax+k*Km*LAI*Vcmax+gamma*h*LAI*Rd^2*VPD-2*gamma*h*LAI*Rd*Vcmax*VPD+gamma*h*LAI*Vcmax^2*VPD+ca^2*k*LAI*g+2*ca*k*Km*LAI*g+k*Km^2*LAI*g+ca*gamma*h*LAI*Rd*VPD*g+gamma*h*Km*LAI*Rd*VPD*g-ca*gamma*h*LAI*Vcmax*VPD*g+2*cp*gamma*h*LAI*Vcmax*VPD*g+gamma*h*Km*LAI*Vcmax*VPD*g-ca*k*LAI*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))+2*hk*k*LAI*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))-k*Km*LAI*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))+2*gamma*h*averA*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))+gamma*h*LAI*Rd*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))-gamma*h*LAI*Vcmax*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))))}
    integrate(integrand, 0, g)$value-w
  }
  u <- uniroot(O, c(0,1), tol=1e-10, extendInt = "yes")
  return(u$root)
}

gw2 <- Vectorize(gw1)
dvs <- read.csv("Derived variables with cost.csv")
MAP <- 365
gw3 <- function(w)gw2(w, ca=400, k=0.1, MAP=MAP, averA=dvs[dvs$ca==400 & dvs$k==0.1 & dvs$MAP==MAP, "A"])

# Figure
windows(8, 6)
par(mgp = c(2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(3.5, 4, 1, 2))
curve(gw3, 0, 1,
      xlab = expression(italic(w)~"(%)"), 
      ylab = expression(italic(g[s])~(mol~m^-2~s^-1)),
      xlim = c(0, 1), ylim = c(0, 0.05),
      cex.lab = 1.3)
