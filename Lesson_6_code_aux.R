mesh <- 1000
lambda <- pi*seq(-mesh,mesh)/mesh # values for plotting

rho.1 <- .4
spec <- 1 + 2*rho.1*cos(lambda) # values of spectral d / fact gamma(0) pres 

ts(spec,start=-1,frequency=mesh) # key here : frequancy can be anything ? integrer ? autre 

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
