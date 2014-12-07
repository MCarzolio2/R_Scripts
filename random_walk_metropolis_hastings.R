###########################################
##### Random Walk Metropolis-Hastings #####
###########################################

rwmh <- function(N,log_target,x_0,prop_sd,burnin=0){
  
  # Returns N samples from log_target
  
  d <- length(x_0)
  samples <- matrix(ncol=N,nrow=d)
  current_x <- x_0
  
  for(i in 1:(N-1+burnin)){
    
    if(i>burnin){
      samples[,i-burnin] <- current_x
    }
    
    prop_x <- rnorm(d,mean=current_x,sd=prop_sd)
    log_mh_ratio <- log_target(prop_x)-log_target(current_x)
    if(log(runif(1))<log_mh_ratio) current_x <- prop_x
    
  }
  samples[,N] <- current_x
  return(samples)
}