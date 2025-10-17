mesh <- 1000
seq(-mesh,mesh)/mesh
lambda <- pi*seq(-mesh,mesh)/mesh # values for plotting

rho.1 <- .4
spec <- 1 + 2*rho.1*cos(lambda) # values of spectral d / fact gamma(0) pres 

ts(spec,start=-1,frequency=mesh) # key here : frequency can be anything ? integrer ? autre 

plot(ts(spec,start=-1,frequency=mesh),xlab="Cycles",ylab="") 
abline(h=0,col=2) # no effect



maq.spec <- function(ma.acf,mesh) # args: macf = vect (gamma(0), rho(1) à rho(q)) / pas tracage
{
  q <- length(ma.acf)-1
  lambda <- pi*seq(0,mesh)/mesh # seq de 0 à 1, avec un pas de mesh fois pi
  spec <- ma.acf[1]*cos(0*lambda) #
  if(q > 0)
  {
    for(k in 1:q)
    {
      spec <- spec + 2*ma.acf[k+1]*cos(k*lambda)
    }
  }
  return(spec)
}

# $q=2$, $\gamma(0)=1$, $\rho(1)= 0$, and $\rho(2) =.6$.

spec <- maq.spec(c(1,0,.6),mesh)
plot(ts(spec,start=0,frequency=mesh),xlab="Cycles",ylab="")
abline(h=0,col=2)


# spectral density arma (p,q) ---------------------------------------------

spec <- NULL
mesh <- 1000
rho <- .8
omega <- pi/6
ar.coef <- c(2*rho*cos(omega),-1*rho^2)
ar.coef
ma.coef <- -1*rho*cos(omega)
ma.coef


armapq.spec <- function(ar.coef,ma.coef,sigma,mesh)
{
  p <- length(ar.coef)
  q <- length(ma.coef)
  lambda <- pi*seq(0,mesh)/mesh
  lambda
  spec.ar <- rep(1,mesh+1) 
  spec.ar
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


spec <- NULL
mesh <- 1000
rho <- .8
omega <- pi/6
ar.coef <- c(2*rho*cos(omega),-1*rho^2)
ma.coef <- -1*rho*cos(omega)
spec <- armapq.spec(ar.coef,ma.coef,1,mesh) #sigma = 1
length(spec)
plot(ts(spec,start=0,frequency=mesh),xlab="Cycles",ylab="",main="") # time series, freq = mesh




# Ideal Low pass ----------------------------------------------------------

mu <- pi/5
mesh <- 1000
lambda <- pi*seq(0,mesh)/mesh
psi.frf <- rep(0,mesh+1) # filter coefficients values
psi.frf[lambda <= mu] <- 1 # 1, 0 sinon
plot(ts(psi.frf,start=0,frequency=mesh),ylab="",xlab="Cycles")


## Eigenvalues of an MA(1) Toeplitz Matrix.

theta <- .8
# sigma <- 1

n <- 10
lambda <- 2*pi*seq(0,n-1)/n # fourier freq direct, idem que entre -pi et pi
Gamma <- toeplitz(c(1+theta^2, theta, rep(0,n-2))) # on declare toutes les valeurs de la 1ere ligne
Gamma
eigen(Gamma)$values
# avec calcul densité spectrale, pkoi ordre inverse pile ?
rev(sort(1+theta^2 + 2*theta*cos(lambda)))

n <- 20
lambda <- 2*pi*seq(0,n-1)/n
Gamma <- toeplitz(c(1+theta^2, theta, rep(0,n-2)))
eigen(Gamma)$values
rev(sort(1+theta^2 + 2*theta*cos(lambda)))

n <- 30
lambda <- 2*pi*seq(0,n-1)/n
Gamma <- toeplitz(c(1+theta^2, theta, rep(0,n-2)))
eigen(Gamma)$values
rev(sort(1+theta^2 + 2*theta*cos(lambda)))

## Eigenvalues of an MA(1) Inverse Toeplitz Matrix.

theta <- .8

n <- 10
lambda <- 2*pi*seq(0,n-1)/n
Gamma <- toeplitz(c(1+theta^2, theta, rep(0,n-2)))
eigen(solve(Gamma))$values
1/sort(1+theta^2 + 2*theta*cos(lambda))

n <- 20
lambda <- 2*pi*seq(0,n-1)/n
Gamma <- toeplitz(c(1+theta^2, theta, rep(0,n-2)))
eigen(solve(Gamma))$values
1/sort(1+theta^2 + 2*theta*cos(lambda))

n <- 30
lambda <- 2*pi*seq(0,n-1)/n
Gamma <- toeplitz(c(1+theta^2, theta, rep(0,n-2)))
eigen(solve(Gamma))$values
1/sort(1+theta^2 + 2*theta*cos(lambda))

