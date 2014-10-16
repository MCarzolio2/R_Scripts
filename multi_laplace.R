# Multivariate Laplace

rMultiLaplace <- function(n,lambda,psi,mu=0){
  
  # Determinant of psi must be 1
  if(abs(det(psi)-1)>sqrt(.Machine$double.eps)) cat('Error: det(psi)!=1 \n')
  
  # Length of mu must agree with dim of psi
  if(!length(mu)%in%c(dim(psi),1)) cat('Error: !length(mu)%in%c(dim(psi),1)')
  
  z <- rexp(n,rate=1/lambda)
  choleski_decomp <- t(chol(psi))
  
  result <- sapply(z,FUN=function(x) mu + sqrt(x)*choleski_decomp%*%rnorm(p))
  
  return(result)
  
}

logMultiLaplacePDF <- function(X,mu,psi,lambda){
  
  # x is an MxN matrix of N observations of a multivariate Laplace
  # see: http://eo.uit.no/publications/te-spl-06.pdf
  
  N <- ncol(X)
  
  term1 <- N*(-p/2*log(2*pi))
  term2 <- N*(log(2)-log(lambda))
  term3 <- 0
  
  p <- length(X[,i])
  
  for(i in 1:N){
    
    q <- t(X[,i]-mu)%*%solve(psi)%*%(X[,i]-mu)
    
    bessel <- besselK(sqrt(2/lambda*q),nu=p/2-1,expon.scaled=FALSE)
    term3 <- term3 + log(bessel)-(p/2-1)/2*log(lambda/2*q)
  }
  
  
  return(term1+term2+term3)
  
}