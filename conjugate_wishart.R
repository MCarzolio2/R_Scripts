# rWishart draw for Gibbs sampling in conjugate inference on precision matrices

conj_wishart <- function(y,mu,nu,D_inv){
  
  # y is an N x M matrix of observations, 
  #   where each column corresponds to a 
  #   multivariate normal draw of size N
  # mu is an N mean vector
  # nu is the prior df
  # D_inv is the inverse of the prior scale matrix
  
  # Check inputs:
  D_inv <- as.matrix(D_inv)
  if(diff(dim(D_inv))!=0) cat("Error: Asymmetric scale matrix \n")
  
  if(class(y)!="matrix") cat("Error: class(y)!=matrix \n")
  N <- nrow(y)
  M <- ncol(y)
  
  if(nu+1<N) cat("Error: df + 1 < N \n")
  if(!length(mu)%in%c(1,N)) cat("Error: !length(mu)%in%c(1,N) \n")
  
  # Compute posterior draw:
  errs <- matrix(apply(y,2,FUN=function(x) x-mu),ncol=M)
  
  dev_mat <- matrix(0,ncol=N,nrow=N)
  
  for(i in 1:M) dev_mat <- dev_mat + sum(errs[,i]*errs[,i],na.rm=TRUE)
  
  posterior_D <- solve(D_inv+dev_mat)
  
  return(drop(rWishart(n = 1,df = nu + M,posterior_D)))
  
}