###########################################
####### Weighted Particle Tempering #######
###########################################

source('./random_walk_metropolis_hastings.R')

if(!'parallel'%in%installed.packages()[,1]) install.packages('parallel')

particle_weights <- function(log_target,particles,delta){
  
  log_evals <- (1-delta)*unlist(lapply(particles,log_target))
  weights <- exp(log_evals)
  return(weights)
  
}

wpt <- function(N,log_target,x_0,p,nu,prop_sd,delta=NULL,burnin=0,parallel=TRUE,cl=NULL,keep.open=FALSE,high_dim_exp=FALSE){
  d <- length(x_0)
  if(is.null(delta)) delta <- nu
  if(parallel){
    require(parallel)
    if(is.null(cl)){
      cl <- makeCluster(p)
    } 
    varlist <- c('N','log_target','x_0','prop_sd','rwmh','nu','burnin','d')
    clusterExport(cl,varlist,envir=environment())
    if(high_dim_exp) {
      particles <- clusterCall(cl,function() matrix(rwmh(N+burnin,log_target = function(x) nu*log_target(x),x_0=runif(d,0,10),prop_sd)[,sample(1:(N+burnin))],ncol=N+burnin,nrow=d))
      } else{
        particles <- clusterCall(cl,function() matrix(rwmh(N+burnin,log_target = function(x) nu*log_target(x),x_0,prop_sd)[,sample(1:(N+burnin))],ncol=N+burnin,nrow=d))
      }
    if(!keep.open) {
      stopCluster(cl)
      rm(cl)
    }
  } else {
    particles <- list()
    for(j in 1:p){
      if(high_dim_exp){
        particles <- c(particles,list(matrix(rwmh(N+burnin,log_target = function(x) nu*log_target(x),x_0=runif(d,0,10),prop_sd)[,sample(1:(N+burnin))],ncol=N+burnin,nrow=d)))
      } else {
        particles <- c(particles,list(matrix(rwmh(N+burnin,log_target = function(x) nu*log_target(x),x_0,prop_sd)[,sample(1:(N+burnin))],ncol=N+burnin,nrow=d))) 
      }
      
    }
  }
  samples <- matrix(ncol=N,nrow=d)
  weights <- sapply(1:(N+burnin),FUN=function(i) particle_weights(log_target,lapply(particles,FUN=function(x) x[,i]),delta))
  pos_weights <- sum(weights[,1])>0
  sample_ind <- if(pos_weights) sample(1:p,1,prob=weights[,1]) else sample(1:p,1)
  
  current_sample <- particles[[sample_ind]][,1]
  log_mh_ratio_vector <- vector('numeric',N)
  accept_rate <- 0
  
  for(i in 1:(N-1+burnin)){
    if(i>burnin){
      samples[,i-burnin] <- current_sample
    }
    pos_weights <- sum(weights[,i])>0
    sample_ind <- if(pos_weights) sample(1:p,1,prob=weights[,i]) else sample(1:p,1)
    prop_sample <- particles[[sample_ind]][,i]
    if(!pos_weights) print(i)
    dub <- weights[,i]/sum(weights[,i])
    
    log_mh_ratio <- if(pos_weights) (delta-nu)*(log_target(prop_sample)-log_target(current_sample))+log(sum(weights[,i]))-log(sum(particle_weights(log_target,list(current_sample),delta),weights[-sample_ind,i])) else (delta-nu)*(log_target(prop_sample)-log_target(current_sample))+log(sum(weights[,i]))
    
    log_mh_ratio_vector[i] <- log_mh_ratio
    
    if(log(runif(1))<log_mh_ratio){
      current_sample <- prop_sample 
      accept_rate <- 1+accept_rate
    } else{
      current_sample <- rwmh(1,log_target = function(x) log_target(x),current_sample,prop_sd*nu)
    }
  }
  samples[,N] <- current_sample
  if(!high_dim_exp) return(samples) else return(list(samples,log_mh_ratio_vector,accept_rate/(N+burnin)))
}