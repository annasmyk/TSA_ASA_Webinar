# multp polinomials 

# see polynomials in R functions

# coeffs pol 1, coeffs pol 2

a<- c(1,2,3)
b<- c(2,3,4)

polymul <- function(a,b) 
{
  bb <- c(b,rep(0,length(a)-1))
  B <- toeplitz(bb)
  B[lower.tri(B)] <- 0
  aa <- rev(c(a,rep(0,length(b)-1)))
  prod <- B %*% matrix(aa,length(aa),1)
  return(rev(prod[,1]))
}


theta <- c(1,.4,.2,-.3)
gamma <- polymul(theta,rev(theta))
print(gamma)
