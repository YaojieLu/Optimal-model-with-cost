
#options(digits=20)
Af <- function(gs,
               ca=400, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))

mf <- function(gs,
               a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
               pe=-1.58*10^-3, b=4.38, h2=h/1000, kmax=5, c=5.71, d=10.05, h3=10){
  
  ps <- pe*w^(-b)
  kf <- function(x)kmax*exp(-(-x/d)^c)
  f1 <- function(x)-((ps-x)*h2*kf(x))
  minimum <- optimize(f1, c(-20,0))$minimum
  f2 <- function(x)((ps-x)*h2*kf(x)-h*VPD*gs)^2
  px <- optimize(f2, c(minimum, ps))$minimum
  
  k <- kf(px)
  PLC <- 1-k/kmax
  res <- h3*PLC
  return(res)
}

Bf <- Vectorize(function(gs)Af(gs)-mf(gs))

B0f <- function(ca=400){
  f1 <- function(gs)Bf(gs)
  res <- uniroot(f1, c(0, 1), tol=.Machine$double.eps)
  return(res$root)
}

WUEf <- function(gs)Bf(gs)/gs

gsmaxf1 <- function(w,
                   a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                   pe=-1.58*10^-3, b=4.38, h2=h/1000, kmax=5, c=5.71, d=10.05){
  ps <- pe*w^(-b)
  kf <- function(x)kmax*exp(-(-x/d)^c)
  f1 <- function(x)-((ps-x)*h2*kf(x))
  pxmin <- optimize(f1, c(-20,0), tol = .Machine$double.eps)$minimum
  gsmax <- (ps-pxmin)*h2*kf(pxmin)/(h*VPD)
  return(gsmax)
}
gsmaxf2 <- Vectorize(gsmaxf1)

# Figure
windows(8, 6)
par(mgp = c(2.2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(4, 4, 2.5, 1), mfrow=c(1,1))
w <- 0.15
curve(WUEf, B0f(), gsmaxf2(w))
#curve(WUEf, B0f(), gsmaxf2(w),
#      main=expression(italic(m(w,g[s]))),
#      xlab=expression(italic(g[s])), 
#      ylab="WUE",
#      xlim=c(0, 1.6), ylim=c(0, 300),
#      cex.lab=1.3)
#abline(v=B0f())
#w <- 0.2
#curve(WUEf, B0f(), gsmaxf2(w), col="red", add=T)
#w <- 1
#curve(WUEf, B0f(), gsmaxf2(w), col="blue", add=T)
#legend("topright", c("0.15", "0.2", "1"), lty=c(1, 1, 1), col=c("black", "red", "blue"), title=expression(italic(w)))
