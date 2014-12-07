###########################################
####### Weighted Particle Tempering #######
###########################################

source('./random_walk_metropolis_hastings.R')

if(!'parallel'%in%installed.packages()[,1]) install.packages('parallel')

particle_weights <- function(log_target,particles,nu){
  
  log_evals <- (1-nu)*unlist(lapply(particles,log_target))
  max_log_evals <- max(log_evals)
  adjusted_log_evals <- log_evals-max_log_evals
  weights <- exp(adjusted_log_evals)/sum(exp(adjusted_log_evals))
  return(weights)
  
}

wpt <- function(N,log_target,x_0,p,nu,prop_sd,burnin=0,parallel=TRUE,cl=NULL,keep.open=FALSE){
  
  if(parallel){
    require(parallel)
    if(is.null(cl)){
      cl <- makeCluster(p)
      varlist <- c('N','log_target','x_0','prop_sd','rwmh')
      clusterExport(cl,varlist)
    } 
    
    particles <- clusterCall(cl,function() matrix(rwmh(N+burnin,log_target = function(x) nu*log_target(x),x_0,prop_sd)[,sample(1:N)],ncol=N+burnin))
    if(!keep.open) {
      stopCluster(cl)
      rm(cl)
    }
  } else {
    particles <- list()
    for(j in 1:p){
      particles <- c(particles,list(matrix(rwmh(N+burnin,log_target = function(x) nu*log_target(x),x_0,prop_sd)[,sample(1:N)],ncol=N+burnin)))
    }
  }

  d <- length(x_0)
  
  samples <- matrix(ncol=N,nrow=d)
  weights <- sapply(1:(N+burnin),FUN=function(i) particle_weights(log_target,lapply(particles,FUN=function(x) x[,i]),nu))
  
  sample_ind <- sample(1:p,1,prob=weights[,1])
  current_sample <- particles[[sample_ind]][,1]
  
  for(i in 1:(N-1+burnin)){
    if(i>burnin){
      samples[,i-burnin] <- current_sample
    }
    
    sample_ind <- sample(1:p,1,prob=weights[,i])
    prop_sample <- particles[[sample_ind]][,i]
    
    pseudo_samples <- unlist(lapply(particles,FUN=function(x) x[,sample(1:i,1)])[-sample(1:p,1)])
    pseudo_weights <- particle_weights(log_target,list(current_sample,pseudo_samples),nu)
    
    log_mh_ratio <- (1-nu)*(log_target(prop_sample)-log_target(current_sample))+log(pseudo_weights[1])-log(weights[sample_ind,i])
      
    if(log(runif(1))<log_mh_ratio) current_sample <- prop_sample
  }
  samples[,N] <- current_sample
  return(samples)
}