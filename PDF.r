
dff <- function(k, MAP,
                 Amax=15.85545680, Kh=0.03530976, VPD=0.02, p=43200, LAI=1, nZ=0.5,
                 a=1.6, l=1.8e-5, MDP=MAP/365, gamma=1/((MDP/k)/1000)*nZ, h=l*a*LAI/nZ*p){
  
  Ev <- function(w){h*VPD*gsw3(w)}
  rEv <- function(w){1/Ev(w)}
  integralEv <- Vectorize(function(w){integrate(rEv, 0, w)$value})
  fnoc <- function(w){1/Ev(w)*exp(-gamma*w+k*integralEv(w))}
  integralfnoc <- integrate(fnoc, 0, 1)$value
  c <- 1/(integralfnoc+1/k)
  f <- function(w){c*fnoc(w)}
  pdf <- f(w2)
  #cdf <- cumsum(pdf)/50
  return(pdf)
}

w2 <- seq(1/50, 1, length=50)
df <- dff(k=k, MAP=MAP)

windows(8, 6)
par(mgp=c(2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(2, 4, 1.5, 2), mfrow=c(1,1))

plot(w2, df, 
     ylab="Probability Density",
     xlim=c(0, 1),
     ylim=c(0, 2),
     cex.lab=1.3)