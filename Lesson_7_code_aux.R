
# Spectral Distribution of an ARMA Process --------------------------------

## spectral density

armapq.spec <- function(ar.coef,ma.coef,sigma,mesh)
{
  p <- length(ar.coef)
  q <- length(ma.coef)
  lambda <- pi*seq(0,mesh)/mesh
  spec.ar <- rep(1,mesh+1) 
  if(p > 0)
  {
    for(k in 1:p)
    {
      spec.ar <- spec.ar - ar.coef[k]*exp(-1i*lambda*k)
    }
  }
  spec.ma <- rep(1,mesh+1) 
  if(q > 0)
  {
    for(k in 1:q)
    {
      spec.ma <- spec.ma + ma.coef[k]*exp(-1i*lambda*k)
    }
  }
  spec <- sigma^2*Mod(spec.ma)^2/Mod(spec.ar)^2
  return(spec)
}

# ok 

# we get the spectral distribution by numerical integration, via trapezoidal rule.

armapq.distr <- function(ar.coef,ma.coef,sigma,mesh)
{
  
  f.spec <- armapq.spec(ar.coef,ma.coef,sigma,mesh)
  f.spec <- c(rev(f.spec),f.spec[-1]) # pourquoi ça avant d'integrer
  f.distr <- rep(0,2*mesh+1)
  for(i in 1:(2*mesh))
  {	
    f.distr[i+1] <- f.distr[i] + (f.spec[i]+f.spec[i+1])/2
  }
  f.distr <- pi*f.distr/mesh
  return(f.distr)
}


#Then we apply this in the case of the cycle ARMA(2,1). We plot the spectral density
#and the spectral distribution.
sigma<-1
mesh <- 1000
rho <- .9
omega <- pi/6
ar.coef <- c(2*rho*cos(omega),-1*rho^2)
ma.coef <- -1*rho*cos(omega)
spec <- armapq.spec(ar.coef,ma.coef,1,mesh)
plot(ts(c(rev(spec),spec[-1]),start=-1,frequency=mesh),xlab="Cycles",ylab="Spectral Density",main="")
spec <- armapq.distr(ar.coef,ma.coef,1,mesh)
plot(ts(spec,start=-1,frequency=mesh),xlab="Cycles",ylab="Spectral Distribution",main="")


# DFT of an AR(1) ---------------------------------------------------------

# We simulate an AR(1).

arp.sim <- function(n,burn,ar.coefs,innovar) # more general for AR(p)
{
  p <- length(ar.coefs)
  z <- rnorm(n+burn+p,sd=sqrt(innovar))
  x <- z[1:p] # init x, here no constraint on variance to have a stationnary process (because of burn in ?)
  for(t in (p+1):(p+n+burn))
  {
    next.x <- sum(ar.coefs*x[(t-1):(t-p)]) + z[t] # mult elmem by elem 
    x <- c(x,next.x)
  }	
  x <- x[(p+burn+1):(p+burn+n)]
  return(x)
}



phi <- .8
innovar <- 1
n <- 200
x.sim <- arp.sim(n,500,phi,innovar)
plot(ts(x.sim))

# We compute the DFT. We begin with code to compute $Q$.
mesh<-10

get.qmat <- function(mesh)
{
  mesh2 <- floor(mesh/2)
  inds <- seq(mesh2-mesh+1,mesh2)
  Q.mat <- exp(1i*2*pi*mesh^{-1}*t(t(seq(1,mesh)) %x% inds))*mesh^{-1/2} # rangé dans la matrice direct dans l'ordre 
  # cf indices, rremplyssage par ligne par defaut 
  return(Q.mat)
}

# exp(1i*2*pi*mesh^{-1}*t(t(seq(1,mesh)) %x% inds))
# 1i*2*pi*mesh^{-1}
# seq(1,mesh) # colonne par defaut
# t(seq(1,mesh)) # ligne 
# t(t(seq(1,mesh))) # vraie colonne, necessaire pour multiplicaion matricelle
# t(t(seq(1,mesh))) %*% inds# mult matricielle classque
# t(t(seq(1,mesh))) %x% inds # produit cartesien
# str(inds)
# inds
# n<-10
# Q.mat <- get.qmat(n)
# Q.mat

phi <- .8
innovar <- 1
n <- 200
x.sim <- arp.sim(n,500,phi,innovar)
plot(ts(x.sim))
Q.mat <- get.qmat(n)



# Then we use Proposition 7.2.7 to get the DFT vector. We plot the real part, imaginary part,
# and the modulus.

x.dft <- Conj(t(Q.mat)) %*% x.sim

plot(ts(Re(x.dft),start=-floor(n/2)+1,frequency=1),
     xlab="Frequency Index",ylab="Real Part")

plot(ts(Im(x.dft),start=-floor(n/2)+1,frequency=1),
     xlab="Frequency Index",ylab="Imaginary Part")

plot(ts(Mod(x.dft)^2,start=-floor(n/2)+1,frequency=1),
     xlab="Frequency Index",ylab="Squared Modulus")


# We compute the DFT of the Wolfer sunspot data, and plot.

wolfer <- read.table("Data/wolfer.dat")
wolfer <- ts(wolfer,start=1749,frequency=12)
plot(wolfer)
n <- length(wolfer)
Q.mat <- get.qmat(n) # ne depend que de n
x.dft <- Conj(t(Q.mat)) %*% wolfer
plot(ts(Re(x.dft),start=-floor(n/2)+1,frequency=1),
     xlab="Frequency Index",ylab="Real Part")
plot(ts(Im(x.dft),start=-floor(n/2)+1,frequency=1),
     xlab="Frequency Index",ylab="Imaginary Part")
plot(ts(log(Mod(x.dft)^2),start=-floor(n/2)+1,frequency=1),
     xlab="Frequency Index",ylab="Log Squared Modulus")

# We see higher values of the squared modulus (the uncentered periodogram) near
# frequency zero.


# We compute the DFT of the Mauna Loa growth rate, and plot.

mau <- read.table("Data/mauna.dat",header=TRUE,sep="")
mau <- ts(mau,start=1958,frequency=12)
mau.gr <- diff(log(mau))
plot(mau)

n <- length(mau.gr)
Q.mat <- get.qmat(n)
x.dft <- Conj(t(Q.mat)) %*% mau.gr
plot(ts(Re(x.dft),start=-floor(n/2)+1,frequency=1),
     xlab="Frequency Index",ylab="Real Part")
plot(ts(Im(x.dft),start=-floor(n/2)+1,frequency=1),
     xlab="Frequency Index",ylab="Imaginary Part")
plot(ts(log(Mod(x.dft)^2),start=-floor(n/2)+1,frequency=1),
     xlab="Frequency Index",ylab="Log Squared Modulus")

# We see higher values of the uncentered periodogram in the shape of peaks, at four
# non-zero frequencies. These correspond to cyclical seasonal effects.

# Phase and Gain for Simple Moving Average. 
lambda <- pi*seq(-1000,1000)/1000
p <- 3
simplema.frf <- sin((p+1/2)*lambda)/((2*p+1)*sin(lambda/2))
simplema.gain <- Mod(simplema.frf)
simplema.phase <- atan(Im(simplema.frf)/Re(simplema.frf))
simplema.delay <- -1*simplema.phase/lambda
plot(ts(Re(simplema.frf),start=-1,frequency=1000),
     xlab="Frequency",ylab="Real Part Frf")
plot(ts(Im(simplema.frf),start=-1,frequency=1000),
     xlab="Frequency",ylab="Imaginary Part Frf")
plot(ts(simplema.gain,start=-1,frequency=1000),
     xlab="Frequency",ylab="Gain")
plot(ts(simplema.phase,start=-1,frequency=1000),
     xlab="Frequency",ylab="Phase")
plot(ts(simplema.delay,start=-1,frequency=1000),
     xlab="Frequency",ylab="Phase Delay")

