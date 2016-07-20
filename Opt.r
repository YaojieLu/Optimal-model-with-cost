
options(digits=20)
library(deSolve)
source("Functions.r")

ca <- 400
k <- 0.05
MAP <- 1000
parms <- c(LAI=1,
           Vcmax=50, cp=30, Km=703, Rd=1,
           a=1.6, nZ=0.5, p=43200, l=1.8e-5, VPD=0.02,
           pe=-1.58*10^-3, b=4.38, c=2.64, d=3.54, kxmax=5, h3=10)

wL <- 0.2

# integralfnoc of PDF
integralfnocf <- function(wL, parms){
  with(as.list(c(wL, parms)), {
    
    h <- l*a*LAI/nZ*p
    gamma <- 1/((MAP/365/k)/1000)*nZ
    
    muopt <- optimize(muf, c(-20, 0), tol=.Machine$double.eps, wL=wL)

    gswf1 <- Vectorize(function(w)gswf(w, muopt$minimum, wL))
    Ef <- function(w){h*VPD*gswf1(w)}
    rEf <- function(w){1/Ef(w)}
    integralrEf <- Vectorize(function(w){integrate(rEf, wL, w, rel.tol=.Machine$double.eps^0.5)$value})
    fnoc <- function(w){1/Ef(w)*exp(-gamma*(w-wL)/(1-wL)+k*integralrEf(w)*1/(1-wL))*1/(1-wL)}
    
    res <- integrate(fnoc, wL, 1, rel.tol=.Machine$double.eps^0.5)
    
    return(res)
  })
}

integralfnoc <- integralfnocf(wL, parms)
cPDF <- 1/(integralfnoc$value+1/k)

# PDF of w
PDFf <- function(w, wL, cPDF, parms){
  with(as.list(c(wL, cPDF, parms)), {
    
    h <- l*a*LAI/nZ*p
    gamma <- 1/((MAP/365/k)/1000)*nZ
    
    muopt <- optimize(muf, c(-20, 0), tol=.Machine$double.eps, wL=wL)
    
    gswf1 <- Vectorize(function(w)gswf(w, muopt$minimum, wL))
    Ef <- function(w){h*VPD*gswf1(w)}
    rEf <- function(w){1/Ef(w)}
    integralrEf <- Vectorize(function(w){integrate(rEf, wL, w, rel.tol=.Machine$double.eps^0.5)$value})
    res <- cPDF/Ef(w)*exp(-gamma*(w-wL)/(1-wL)+k*integralrEf(w)*1/(1-wL))*1/(1-wL)
    return(res)
  })
}

# Average A
averAf <- function(wL, cPDF, parms){
  with(as.list(c(wL, cPDF, parms)), {
    
    h <- l*a*LAI/nZ*p
    gamma <- 1/((MAP/365/k)/1000)*nZ
    
    muopt <- optimize(muf, c(-20, 0), tol=.Machine$double.eps, wL=wL)
    
    gswf1 <- Vectorize(function(w)gswf(w, muopt$minimum, wL))
    Af1 <- Vectorize(function(w)Af(gswf1(w), ca))
    f1 <- function(w)Af1(w)*PDFf(w, wL, cPDF, parms)
    res <- integrate(f1, wL, 1, rel.tol=.Machine$double.eps^0.5)
    return(res)
  })
}

averAf(wL, cPDF, parms)
